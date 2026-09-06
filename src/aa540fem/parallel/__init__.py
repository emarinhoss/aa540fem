"""MPI support (optional: ``mpi4py``, and ``petsc4py`` for the distributed solver).

Two stages, see ``docs/parallel.md``:

- replicated assembly with a distributed factorisation
  (:mod:`aa540fem.parallel.distributed`): every rank assembles the whole
  system, MUMPS factorises across the ranks, the solution is gathered;
- domain decomposition (:mod:`aa540fem.parallel.partition`): METIS
  partition of the element graph with a one-layer halo, the starting point
  of a fully distributed assembly.
"""

from aa540fem.parallel.comm import is_parallel, rank, rank0_print, size, world

__all__ = ["is_parallel", "rank", "rank0_print", "size", "world"]


def __getattr__(name):                        # lazy: flow.py needs petsc4py only when used
    if name in ("DistributedFlowSystem", "solve_flow_distributed",
                "solve_flow_transient_distributed"):
        from aa540fem.parallel import flow

        return getattr(flow, name)
    raise AttributeError(name)
