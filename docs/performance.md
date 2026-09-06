# Performance record

Timings from `benchmarks/bench_flow.py` on the committed meshes: the
steady stabilised cylinder (10.8k nodes, 25k unknowns, 6 Newton
iterations), 20 theta steps and 20 fixed RK45 steps on the same mesh, and
two outer iterations of the RANS coupling on the laminar plate mesh.  The
columns attribute the wall time to element assembly, matrix factorisation,
triangular solves and everything else (Python overhead, sparse conversions,
output).  All runs on the same 4-core machine; the machine was shared with
other jobs during some of them, so differences below 10 % mean nothing.

Each section is one optimisation phase (see `docs/parallel.md`).

### Baseline (before any optimisation)

x86_64, 2.4.6 numpy

```
case      wall [s]  assembly  factorise   solve   other  note
steady        9.10      2.35       6.08    0.00    0.67  6 Newton it.
theta        73.16     43.44      22.40    0.00    7.32  20 steps
rk45         43.77     37.48       1.82    0.00    4.46  20 steps
rans         33.40      5.31      19.62    0.00    8.47  2 outer it.
```

### Phase 1: residual-only evaluations, fixed sparsity pattern (NumPy, SuperLU)

x86_64, 2.4.6 numpy

```
case      wall [s]  assembly  factorise   solve   other  note
steady        9.07      1.73       6.23    0.00    1.11  6 Newton it.
theta        32.38     10.57      20.24    0.00    1.57  20 steps
rk45          8.45      4.55       1.63    0.00    2.27  20 steps
rans         27.38      3.62      16.07    0.00    7.69  2 outer it.
```

### Phase 3: numba assembly kernels, 4 threads (SuperLU)

x86_64, 2.4.6 numpy

```
case      wall [s]  assembly  factorise   solve   other  note
steady        9.56      0.77       5.70    0.00    3.09  6 Newton it.
theta        29.24      4.37      20.88    0.00    4.00  20 steps
rk45          7.16      0.98       1.72    0.00    4.46  20 steps
rans         26.02      1.11      15.58    0.00    9.33  2 outer it.
```


### Phases 2-3: numba kernels (4 threads) + PETSc/MUMPS factorisation, one process

x86_64, 2.4.6 numpy, 100 % of the machine: 4 threads, 12.0 GB budget, assembly numba, direct solver petsc, linear algebra scipy

```
case      wall [s]  assembly  factorise   solve   other  note
steady        4.81      1.09       2.55    0.13    1.04  6 Newton it.
theta        18.31      7.37       9.90    1.89    0.00  20 steps
rk45          3.24      0.89       0.48    1.52    0.34  20 steps
rans         17.12      1.74      12.80    0.64    1.95  2 outer it.
```

### Phase 5: four MPI ranks, replicated assembly, distributed MUMPS (one thread per rank)

`mpirun -n 4 python benchmarks/bench_flow.py --cases steady theta --threads 1`
on the same 4-core machine (rank 0 shown; every rank holds the full system):

```
case      wall [s]  assembly  factorise   solve   other  note
steady        4.69      1.55       2.09    0.08    0.96  6 Newton it.
theta        16.98      7.01       9.69    1.09    0.00  20 steps
```

On one machine the distributed factorisation across four ranks is no faster
than MUMPS in one process with four threads of assembly: at 24k unknowns the
2D factorisation takes two seconds either way, and the replicated assembly
runs on one thread per rank.  The MPI path is there for problems that do not
fit one process, not to speed this size up.

With the distributed assembly of `parallel/flow.py` (every rank assembles
its own elements and inserts its own rows; the raw table is at the end of
this file) the same four ranks solve the steady cylinder in 1.2 s and 20
theta steps in 1.7 s, against 3.6 s and 5.2 s for the replicated path and
4.8 s / 6.5 s for one process with four threads: on one machine the
domain decomposition is now the fastest way to run the solver, and it is
the path that scales beyond one machine.

### Summary (4-core machine)

| case | baseline | phase 1 (NumPy, SuperLU) | phase 3 (numba, SuperLU) | phases 2-3 (numba, MUMPS) |
|---|---|---|---|---|
| steady Newton, 6 iterations | 9.1 s | 9.1 s | 9.6 s | 4.8 s |
| theta, 20 steps | 73 s | 32 s | 29 s | 18 s |
| RK45, 20 steps | 44 s | 8.5 s | 7.2 s | 3.2 s |
| RANS, 2 outer iterations | 33 s | 27 s | 26 s | 17 s |

What remains after these phases is the factorisation itself (MUMPS, one
process) and, in the theta scheme, the per-step refactorisation; the
Krylov path (`solve_flow(method="fieldsplit")`) removes the factorisation
at the price of iterations that grow with the Reynolds number (see the
Krylov section below).

### Validation cases with the final defaults (numba, MUMPS, one process)

| case | before | after | result |
|---|---|---|---|
| turbulent flat plate, Re 1e6 (`turbulent_flat_plate.py`) | 9 min | 4 min | Cf 0.00327 vs Coles-Fernholz 0.00323, unchanged |
| cylinder shedding, 1600 theta steps (`cylinder_shedding.py`) | 77 min | 19 min | St 0.3013, C_D,max 3.2303, C_L,max 0.9958, unchanged |
| cylinder shedding, twice-finer mesh (`cylinder_bl_fine.msh`, 41k nodes, 92k unknowns) | 74 min (one factorisation per step) | 39 min (factorisation carried across steps: 699 for 1600 steps) | St 0.3014, C_D,max 3.2320, C_L,max 0.9972, identical |
| cylinder shedding, dt 0.0025 (3200 steps, coarse mesh) | - | 29 min | St 0.3016, C_D,max 3.2271, C_L,max 0.9885 |

The refined shedding run is the data point for the scaling of the direct
solver: 3.8x the unknowns cost 3.9x the time per step (2.8 s against
0.71 s).  Measured on the two Jacobians (numba assembly, MUMPS, 4 threads):

| mesh | unknowns | assembly | factorise | solve | factor entries |
|---|---|---|---|---|---|
| `cylinder_bl.msh` | 24281 | 0.04 s | 0.36 s | 0.008 s | 4.3e6 |
| `cylinder_bl_fine.msh` | 92330 | 0.11 s | 2.66 s | 0.040 s | 2.0e7 |

The factorisation grew 7.4x, as the n^1.5 of a 2D sparse LU predicts
(3.8^1.5 = 7.4), the assembly linearly, so with one factorisation per step
a theta step on the fine mesh was 95 % factorisation; the factors take
0.25 GB.  The theta scheme therefore now carries the factorised Jacobian
from one step into the next (`reuse_jacobian=True`, the default) and
refreshes it only when a modified-Newton iteration contracts the residual
by less than a factor 3: the fine-mesh run needs 699 factorisations for
1600 steps and 39 minutes instead of 74, with the converged states, and
hence the forces, unchanged (the residual is always the true one).  The
coarse mesh, where the factorisation was 60 % of a step, gains
proportionally less (its 19 minutes were measured with one factorisation
per step).  Beyond this the Krylov path is the remaining lever.

### Krylov solver hardening (fieldsplit, `linalg/krylov.py`)

FGMRES iterations to a relative residual of 1e-8 on the eliminated
Jacobians at the converged states of the validation cases, 4 threads;
"mass" is the viscosity-scaled pressure mass matrix as Schur-complement
preconditioner (the previous default), "LSC" the least-squares commutator;
the velocity block is solved exactly (LU, to isolate the Schur
approximation), by ILU(1) (the CPU default) or by two GAMG V-cycles (the
device choice).  "-" means no convergence within 300 iterations.

| Jacobian | unknowns | direct solve | mass + LU | LSC + LU | LSC + ILU(1) | LSC + GAMG |
|---|---|---|---|---|---|---|
| cavity Re 100, Q2/Q1 | 5427 | 0.05 s | 75 | 36 | 34 (0.06 s) | 35 |
| cavity Re 100, stabilised | 5427 | 0.05 s | 77 | 33 | 35 (0.07 s) | 39 |
| cavity Re 1000, stabilised | 9539 | 0.10 s | - | 58 | 75 (0.19 s) | 98 |
| cylinder Re 20, stabilised, steady | 24281 | 0.36 s | - | 76 | 133 (0.95 s) | 142 |
| cylinder Re 100, theta step dt 0.005 | 24281 | 0.36 s | - | 25 | 28 (0.43 s) | 28 |

Before this the solver used GMRES restart 50, ILU(0) with a single sweep,
and the mass matrix: 300 to 500 iterations without convergence on every
case above.  Three things mattered: the Schur approximation (the mass
matrix ignores the grad-div term, which is hundreds of times the viscosity
on the cylinder, and the convection), the restart length (cavity Re 100:
156 to 75 iterations from restart 50 to 200) and ILU(0), which stalls on
the steady cylinder where ILU(1) converges.  Whole solves through the
Krylov path reproduce the direct ones (cavity Re 1000 to 2e-12, cylinder
Re 20 to 7e-13, 20 theta steps of the cylinder to 3e-9) with the same
Newton counts; in 2D they take 1.1 to 2.7 times the direct time, as
expected at these sizes.


### Distributed assembly, 4 ranks x 1 thread (bench_mpi.py)

100 % of the machine: 1 threads x 4 ranks, 12.0 GB budget, assembly numba, direct solver petsc, linear algebra scipy

```
case                           ranks  wall [s]  note
steady replicated                  4      3.58  6 Newton it.
steady distributed direct          4      1.18  6 Newton it., |dU|/|U| = 1e-15
steady distributed fieldsplit      4      2.84  6 Newton it., |dU|/|U| = 4e-13, 131 FGMRES it. max
theta 20 steps replicated          4      5.19  9 factorisations
theta distributed direct           4      1.72  9 factorisations, |dU|/|U| = 7e-14
theta distributed fieldsplit       4      6.96  9 factorisations, |dU|/|U| = 9e-10
```

### 3-D: Schaefer-Turek 3D-1Z (Re 20), one process, 4 threads

| mesh | unknowns | solver | Newton it. | wall | peak RSS | C_D | C_L | dp |
|---|---|---|---|---|---|---|---|---|
| `cyl3d_coarse.msh` | 81k | MUMPS | 7 | 268 s | 4.0 GB | 6.2397 | -0.0164 | 0.1758 |
| `cyl3d_coarse.msh` | 81k | fieldsplit / LSC, ILU(1) | 7 | 235 s | 3.1 GB | 6.2397 | -0.0164 | 0.1758 |
| `cylinder3d_tet10.msh` | 132k | MUMPS | 7 | 592 s | 6.6 GB | 6.2273 | 0.0141 | 0.1742 |
| `cyl3d_bl.msh` (cylinder refined, channel coarsened) | 99k | MUMPS, Newton then PTC | 6 | 492 s | | 6.3770 | 0.0540 | 0.1806 |
| `cyl3d_fine2.msh` | 227k | fieldsplit / LSC, ILU(1) | 7 | 992 s | 8.5 GB | 6.2073 | 0.0138 | 0.1718 |

Reference 6.185 / 0.0094 / 0.1710.  In 3-D the factorisation dominates a
Newton iteration (about 40 s of the 38 s per iteration at 81k unknowns,
80 s at 132k) and its memory grows faster than the problem (4.0, 6.6 GB);
the Krylov solver takes about the same time per iteration and less memory,
and is the only option here beyond 150k unknowns on this 15 GB machine.
Newton from the Stokes start needs heavy damping on the first steps and
on the cylinder-refined mesh stalls, where the pseudo-transient
continuation takes over (three steps to |R| = 1e-5, then Newton).
