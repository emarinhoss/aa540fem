"""Conductivity tensor and heat source.  Port of ``conductivity_and_forcing.m``.

This is one of the two files meant to be edited by the user (the other is
``main.py``).  Write the expressions with NumPy so they work on arrays.
"""

from __future__ import annotations

import numpy as np  # noqa: F401  (available for user expressions below)


def conductivity_and_forcing(x, y):
    """Return ``(kxx, kxy, kyx, kyy, f)`` evaluated at points ``(x, y)``.

    ``kij`` are the components of the anisotropic conductivity tensor and
    ``f`` is the volumetric heat source in

        div( kappa . grad T ) + f = 0

    Each return value may be a scalar or an array broadcastable to ``x``.
    Add a parameter named ``T`` (``def conductivity_and_forcing(x, y, T)``)
    for temperature-dependent coefficients, or ``t`` for time dependence.
    """
    kxx = 1.0
    kxy = 0.0
    kyx = 0.0
    kyy = 1.0

    f = 0.0
    # f = (np.sin(np.pi * x) * np.sin(np.pi * y)) ** 40
    return kxx, kxy, kyx, kyy, f
