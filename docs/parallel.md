# Running faster: threads, direct solvers, MPI and GPUs

This note says where the solver spends its time, what the package does
about it, which option fits which machine, and what the MPI and GPU paths
are and are not.  Measured numbers are in [`performance.md`](performance.md).

## 1. Where the time goes

Every Newton iteration or time step does three things: assemble the
residual and (when needed) the Jacobian element by element, factorise the
Jacobian (a sparse LU of an indefinite saddle-point matrix), and solve with
the factors.  Before any optimisation the 10.8k-node cylinder case spent
65 % of a steady solve in the single-threaded SuperLU factorisation and the
transient runs 60-85 % in assembly, because every Runge-Kutta stage and
every Newton trial point assembled a Jacobian it never used.

## 2. What the code does now

1. **Fixed sparsity pattern** (`backends/pattern.py`).  The union of the
   element dof blocks is built once; every matrix (viscous, mass,
   divergence, the two Jacobian parts, every combination and the
   pseudo-time-step modifications) is a data vector on it, and combining
   them is vector arithmetic.  Element kernels return element-local arrays
   that a `ScatterPlan` sums in a fixed order, so the assembled data are
   bitwise identical for any number of threads.
2. **Residual-only evaluations.**  `momentum_terms(jacobian=False)` for the
   Runge-Kutta stages, the explicit part of the theta scheme and the trial
   points of the damped Newton and pseudo-transient iterations
   (`residual=` callable of `newton_iterate` / `pseudo_transient`).  The
   theta scheme also carries its factorised Jacobian from one time step
   into the next (`reuse_jacobian=True`) and refactorises only when the
   modified Newton iteration stops contracting by a factor 3 per
   iteration, which halves the cost of a factorisation-bound run.
3. **Threaded assembly** (`backends/numba_kernels.py`, `numba_sa.py`).
   The momentum and Spalart-Allmaras element loops are numba `prange`
   kernels; the NumPy einsum kernels remain as the reference and fallback
   (`FlowAssembler(backend="numpy")`, `AA540FEM_ASSEMBLY=numpy`), and the
   test suite checks both agree to 1e-12.
4. **Direct-solver backends** (`linalg/direct.py`).  `factorise(A)` uses
   SuperLU (SciPy, one thread) or PETSc's LU with MUMPS when `petsc4py` is
   installed with it; MUMPS uses the BLAS threads and OpenMP, and under
   `mpirun` it factorises across ranks.  A residual check falls back to
   SuperLU if a backend mishandles the zero pressure block.
5. **Hardware probe and run configuration** (`hardware.py`, `cli.py`).
   `probe()` reports cores (physical, usable under affinity/cgroup limits),
   memory, NVIDIA GPUs, MPI and the installed backends; `RunConfig`
   carries the thread count, memory budget and backend choices and applies
   the thread counts to numba, the BLAS and OpenMP.  Every example takes
   `--machine {50,100}`, `--threads`, `--assembly`, `--direct`, `--linear`,
   `--memory-budget`, `--no-prompt`; on a terminal it asks once whether to
   use 50 % or 100 % of the machine (Enter = 100 %).  `AA540FEM_*`
   environment variables set the same things for batch jobs, and
   `python -m aa540fem.cli` prints the probe.  Before a factorisation,
   `RunConfig.check_memory` compares the estimated factor size with the
   budget and warns above 80 % or refuses above 100 %.

## 3. Which option for which machine

| machine | best option | why |
|---|---|---|
| laptop or desktop, no GPU, 2D cases up to a few 1e5 unknowns | numba threads plus MUMPS (`pip install aa540fem[fast]`), `--machine 100` when the machine is yours, `--machine 50` when you share it | a 2D sparse LU costs O(n^1.5) time and fits in memory; no communication, nothing to partition |
| workstation with an NVIDIA GPU | the same CPU path; the GPU only through the PETSc Krylov path (`--linear petsc-cuda`, section 5) | sparse LU does not map to GPUs, and after step 3 assembly is no longer the bottleneck |
| cluster, 3D, more than about 1e6 unknowns | `mpirun -n R python example.py --threads T` with the PETSc backend: distributed MUMPS first, fieldsplit Krylov beyond that | the factorisation's memory and O(n^2) time in 3D exceed one node |

Rules of thumb: threads x ranks = the cores you chose; MUMPS and numba
never run at the same time, so both get the full thread count; half the
machine is the right choice when a factorisation would otherwise take the
memory an interactive session needs (the prompt shows both budgets).

## 4. MPI

`aa540fem.parallel` (requires `mpi4py`, and `petsc4py` for the solver) has
two stages.

**Replicated assembly, distributed factorisation.**  Every rank reads the
mesh and assembles the whole system (cheap in 2D, and threaded); the PETSc
direct backend inserts each rank's row range into a distributed matrix and
MUMPS factorises it across the ranks; the solution is gathered back to
every rank, so the Newton and pseudo-time loops run identically on all of
them and only rank 0 writes output.  This is the simplest way to put all
cores of one machine, or a few nodes, behind the factorisation, and it
needs no change in the physics code.

**Domain decomposition.**  `parallel/partition.py` partitions the element
graph with METIS (`pymetis`), adds a one-layer halo of ghost nodes and
numbers the owned and ghost dofs per rank; `parallel/distributed.py`
assembles the local element chunk of a rank with the same kernels and
lets PETSc add the off-process contributions (`setValuesCOO` with the
fixed pattern), which is exactly the element-local-then-scatter split of
section 2.  The prototype assembles the viscous matrix this way and checks
it against the serial one; the remaining steps to a fully distributed
Newton are the halo exchange of the state before each assembly, Dirichlet
elimination with `MatZeroRowsColumns`, `allreduce` for the norms and
forces, and `.pvtu` output.

Why not sooner: below roughly 1e6 unknowns the partitioning, halo
exchanges and PETSc setup cost more than the whole serial solve saves, and
a 2D factorisation fits comfortably in one machine's memory.

## 5. GPUs

A GPU helps two things: element assembly (embarrassingly parallel, but
already off the critical path after step 3) and the sparse matrix-vector
products of an iterative solver.  It does not help the sparse LU, which is
what the solver relies on today.  The GPU path is therefore the PETSc
Krylov solver of `linalg/krylov.py` (`solve_flow(method="fieldsplit")`,
`--linear petsc`): flexible GMRES (restart 200) with a block-triangular
(fieldsplit, Schur lower) preconditioner.  The velocity block is
preconditioned by ILU(1) on the CPU or by two V-cycles of algebraic
multigrid (GAMG) where ILU is not available (the device); the pressure
Schur complement `S = C + B F^-1 B^T` by the least-squares commutator
(LSC, `LSCPreconditioner`): `S^-1 ~ L^-1 (B Q^-1 F Q^-1 B^T) L^-1` with
`L = B Q^-1 B^T` and `Q` the lumped velocity mass matrix, `L` factorised
once.  LSC is purely algebraic, so the convection, SUPG and grad-div terms
of `F` are accounted for automatically; the viscosity-scaled pressure mass
matrix (`schur="mass"`) is kept for Stokes-like problems.  Measured on the
converged Jacobians of the validation cases (`docs/performance.md`): 34 to
133 iterations, and 28 on the theta-scheme Jacobian of the shedding
cylinder, where the mass-matrix preconditioner did not converge in 300
on three of the five cases.  In 2D the direct solver is still faster
(a factorisation of 24k unknowns takes 0.4 s); the Krylov path is for
the 3D and GPU cases where it is not possible.  Each solver instance uses
its own PETSc options prefix (PETSc's options database is global) and
accepts extra options through `fieldsplit_factory(options=...)`.

`--linear petsc-cuda` switches PETSc's matrix and vector types to
`aijcusparse`/`cuda`, which PETSc detects at build time; the probe reports
whether the installed PETSc has CUDA and whether a device is present, and
the corresponding tests skip where it is not.  Assembly stays on the CPU;
a `numba.cuda` version of the element kernel would be the next step.

## 6. Installing the optional pieces

```
pip install aa540fem[numba]      # threaded assembly
pip install aa540fem[petsc]      # petsc4py + mpi4py: MUMPS, MPI, Krylov (needs an MPI library)
pip install aa540fem[fast]       # both
pip install aa540fem[hw]         # psutil, threadpoolctl, pynvml: a more complete probe
```

Building PETSc from PyPI compiles it (10-30 minutes).  What worked on
Ubuntu 24.04 with the distribution's MUMPS, ScaLAPACK and METIS:

```
apt-get install openmpi-bin libopenmpi-dev libmumps-dev libscalapack-openmpi-dev libmetis-dev
export MPICC=mpicc MPICXX=mpicxx MPIFORT=mpif90
export PETSC_CONFIGURE_OPTIONS="--with-debugging=0 \
  --with-mumps-include=/usr/include --with-mumps-lib=[-ldmumps,-lmumps_common,-lpord] \
  --with-scalapack-include=/usr/include --with-scalapack-lib=[-lscalapack-openmpi] \
  --with-metis-include=/usr/include --with-metis-lib=[-lmetis]"
pip install petsc && pip install --ignore-installed --no-deps petsc4py mpi4py pymetis
```

(`--download-mumps` needs internet access from the build; the distribution
packages `petsc-dev` / `python3-petsc4py-real` are the shortcut when they
match the Python interpreter in use.)  `python -m aa540fem.cli` confirms
what was found: it lists `mumps` among the installed backends.

Measured on this 4-core machine (see `performance.md`): MUMPS factorises
the 24k-unknown cylinder Jacobian in 0.5 s where SuperLU needs 1.8 s, and
the whole steady solve drops from 9.1 s to 4.8 s; four MPI ranks with the
distributed factorisation take the same 4.7 s, so on one machine the
threaded single process is the right choice and MPI is for problems that
do not fit it.
