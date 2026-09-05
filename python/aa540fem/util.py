"""Small helpers shared by the modules."""

from __future__ import annotations

import inspect

import numpy as np


def accepts_time(fn) -> bool:
    """True if ``fn`` is a callable taking a third positional argument for time.

    The third parameter counts if it has no default value or is named ``t``
    or ``time``; ``def kappa(x, y, eps=0.01)`` is therefore *not* treated as
    time dependent, while ``lambda x, y, t: ...`` and ``def f(x, y, t=0)`` are.
    """
    if not callable(fn):
        return False
    try:
        sig = inspect.signature(fn)
    except (TypeError, ValueError):
        return False
    positional = [p for p in sig.parameters.values()
                  if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
    if len(positional) >= 3:
        third = positional[2]
        return third.default is third.empty or third.name in ("t", "time")
    return any(p.kind == p.VAR_POSITIONAL for p in sig.parameters.values())


def call_xyt(fn, x, y, t=0.0):
    """Evaluate ``fn(x, y)`` or ``fn(x, y, t)`` depending on its signature."""
    if accepts_time(fn):
        return fn(x, y, t)
    return fn(x, y)


def values_at(val, x, y, t=0.0):
    """Evaluate a constant or a callable ``val(x, y[, t])`` at points ``(x, y)``."""
    x = np.asarray(x, dtype=float)
    if callable(val):
        return np.broadcast_to(np.asarray(call_xyt(val, x, y, t), dtype=float), x.shape).copy()
    return np.full(x.shape, float(val))
