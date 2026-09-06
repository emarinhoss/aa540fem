"""Small helpers shared by the modules."""

from __future__ import annotations

import inspect

import numpy as np

TEMPERATURE_NAMES = ("T", "temperature")
TIME_NAMES = ("t", "time")


def accepts_time(fn, dim: int = 2) -> bool:
    """True if ``fn`` is a callable taking a positional argument for time after
    the ``dim`` coordinates.

    That parameter counts if it has no default value or is named ``t`` or
    ``time``; ``def kappa(x, y, eps=0.01)`` is therefore *not* treated as
    time dependent, while ``lambda x, y, t: ...`` and ``def f(x, y, t=0)`` are.
    """
    if not callable(fn):
        return False
    try:
        sig = inspect.signature(fn)
    except (TypeError, ValueError):
        return False
    params = list(sig.parameters.values())
    if any(p.name in TIME_NAMES for p in params):
        return True
    positional = [p for p in params if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
    if len(positional) > dim:
        third = positional[dim]
        return third.default is third.empty and third.name not in TEMPERATURE_NAMES
    return any(p.kind == p.VAR_POSITIONAL for p in params)


def accepts_temperature(fn) -> bool:
    """True if ``fn`` has a parameter named ``T`` (or ``temperature``)."""
    if not callable(fn):
        return False
    try:
        sig = inspect.signature(fn)
    except (TypeError, ValueError):
        return False
    return any(name in TEMPERATURE_NAMES for name in sig.parameters)


def call_coeff(fn, x, y, t=0.0, T=None):
    """Evaluate a coefficient callable with the arguments its signature asks for.

    Parameters after ``x, y`` are matched by name: ``t``/``time`` receive the
    time, ``T``/``temperature`` the temperature (an error if ``T`` is None);
    an unnamed required third positional parameter receives the time, as in
    :func:`accepts_time`.
    """
    return call_coeff_nd(fn, (x, y), t, T)


def call_coeff_nd(fn, coords, t=0.0, T=None):
    """:func:`call_coeff` for ``d`` coordinates: ``fn(x, y[, z][, t][, T])``."""
    d = len(coords)
    try:
        params = list(inspect.signature(fn).parameters.values())
    except (TypeError, ValueError):
        return fn(*coords)
    args = []
    kwargs = {}
    for i, p in enumerate(params):
        if p.kind == p.VAR_POSITIONAL:
            args.append(t)
            break
        if p.kind not in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY):
            continue
        if p.name in TEMPERATURE_NAMES:
            if T is None:
                raise ValueError(f"{getattr(fn, '__name__', 'coefficient')} depends on the "
                                 "temperature T; use the nonlinear solver")
            value = T
        elif p.name in TIME_NAMES or (i == d and p.default is p.empty and p.kind != p.KEYWORD_ONLY):
            value = t
        elif i < d:
            value = coords[i]
        else:
            continue
        if p.kind == p.KEYWORD_ONLY:
            kwargs[p.name] = value
        else:
            args.append(value)
    return fn(*args, **kwargs)


def call_xyt(fn, x, y, t=0.0):
    """Evaluate ``fn(x, y)`` or ``fn(x, y, t)`` depending on its signature."""
    return call_coeff(fn, x, y, t)


def values_at(val, x, y, t=0.0, T=None):
    """Evaluate a constant or a callable ``val(x, y[, t][, T])`` at points ``(x, y)``."""
    return values_at_nd(val, (x, y), t, T)


def values_at_nd(val, coords, t=0.0, T=None):
    """:func:`values_at` for a tuple of ``d`` coordinate arrays."""
    coords = tuple(np.asarray(c, dtype=float) for c in coords)
    x = coords[0]
    if callable(val):
        return np.broadcast_to(np.asarray(call_coeff_nd(val, coords, t, T), dtype=float),
                               x.shape).copy()
    return np.full(x.shape, float(val))


def values_rate(val, x, y, t=0.0, T=None, h: float = 1e-6):
    """Time derivative of a prescribed value by a central difference (zero if constant)."""
    return values_rate_nd(val, (x, y), t, T, h)


def values_rate_nd(val, coords, t=0.0, T=None, h: float = 1e-6):
    """:func:`values_rate` for a tuple of ``d`` coordinate arrays."""
    if not accepts_time(val, len(coords)):
        return np.zeros(np.shape(np.asarray(coords[0], dtype=float)))
    step = h * max(1.0, abs(t))
    return (values_at_nd(val, coords, t + step, T)
            - values_at_nd(val, coords, t - step, T)) / (2 * step)
