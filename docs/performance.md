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
at the price of iterations that grow with the Reynolds number.

