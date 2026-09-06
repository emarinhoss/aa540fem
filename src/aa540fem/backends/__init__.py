"""Assembly backends: fixed sparsity patterns, NumPy reference kernels and
optional threaded (numba) kernels.

The physics modules ask :func:`assembly_backend` which implementation to use;
the answer comes from the run configuration (:mod:`aa540fem.hardware`) or an
explicit argument, and falls back to NumPy when an optional backend is not
installed.
"""

from __future__ import annotations

import importlib.util
import os
import warnings

from aa540fem.backends.pattern import PatternMatrix, ScatterPlan, SparsityPattern, scatter_vector

ASSEMBLY_BACKENDS = ("numpy", "numba")


def numba_available() -> bool:
    """numba is installed and the numba kernels are present."""
    return (importlib.util.find_spec("numba") is not None
            and importlib.util.find_spec("aa540fem.backends.numba_kernels") is not None)


def assembly_backends() -> list[str]:
    """Installed assembly backends, fastest first."""
    return ["numba", "numpy"] if numba_available() else ["numpy"]


def assembly_backend(requested: str | None = None) -> str:
    """Resolve the assembly backend: argument > ``AA540FEM_ASSEMBLY`` > run
    configuration > fastest installed."""
    name = requested or os.environ.get("AA540FEM_ASSEMBLY") or None
    if name is None:
        try:
            from aa540fem.hardware import get_config

            name = get_config().assembly_backend
        except ImportError:
            name = "auto"
    if name in (None, "auto"):
        return assembly_backends()[0]
    if name not in ASSEMBLY_BACKENDS:
        raise ValueError(f"unknown assembly backend {name!r}; expected one of {ASSEMBLY_BACKENDS}")
    if name == "numba" and not numba_available():
        warnings.warn("numba is not installed; using the NumPy assembly backend", stacklevel=2)
        return "numpy"
    return name


__all__ = ["ASSEMBLY_BACKENDS", "PatternMatrix", "ScatterPlan", "SparsityPattern",
           "assembly_backend", "assembly_backends", "numba_available", "scatter_vector"]
