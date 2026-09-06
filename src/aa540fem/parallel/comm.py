"""Communicator helpers that work with or without ``mpi4py``."""

from __future__ import annotations

import importlib.util


class _SerialComm:
    """Stand-in for ``MPI.COMM_WORLD`` in a serial process."""

    rank = 0
    size = 1

    def Get_rank(self):
        return 0

    def Get_size(self):
        return 1

    def Barrier(self):
        return None

    def bcast(self, obj, root=0):
        return obj

    def allreduce(self, value, op=None):
        return value

    def allgather(self, value):
        return [value]


def mpi4py_available() -> bool:
    return importlib.util.find_spec("mpi4py") is not None


def world():
    """``MPI.COMM_WORLD`` when mpi4py is installed, else a serial stand-in."""
    if mpi4py_available():
        from mpi4py import MPI

        return MPI.COMM_WORLD
    return _SerialComm()


def rank() -> int:
    return world().Get_rank()


def size() -> int:
    return world().Get_size()


def is_parallel() -> bool:
    return size() > 1


def rank0_print(*args, **kwargs):
    if rank() == 0:
        print(*args, **kwargs)
