"""Explicit Runge-Kutta integration: Dormand-Prince 5(4) with step control.

The pair used by MATLAB's ``ode45``: seven stages, first-same-as-last, the
5th-order solution is advanced and the embedded 4th-order one gives the
error estimate.  Works on any ``dy/dt = rhs(t, y)`` with NumPy vectors.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# Dormand & Prince (1980) tableau
DP_C = np.array([0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0])
DP_A = [
    [],
    [1 / 5],
    [3 / 40, 9 / 40],
    [44 / 45, -56 / 15, 32 / 9],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
    [35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
]
DP_B = np.array([35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0.0])
DP_B4 = np.array([5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40])
DP_E = DP_B - DP_B4
ORDER = 5


@dataclass
class RKResult:
    times: np.ndarray
    states: list
    info: dict = field(default_factory=dict)

    @property
    def y(self):
        return self.states[-1]


def rk45(rhs, y0, t0: float, t_end: float, dt: float, *, rtol: float = 1e-4,
         atol: float = 1e-6, dt_max: float | None = None, dt_min: float = 1e-12,
         adaptive: bool = True, output_interval: float | None = None, store_every: int = 1,
         error_mask=None, on_accept=None, store_transform=None,
         verbose: bool = False) -> RKResult:
    """Integrate ``dy/dt = rhs(t, y)`` from ``t0`` to ``t_end``.

    Parameters
    ----------
    dt              : initial step (the fixed step when ``adaptive=False``).
    rtol, atol      : error control; a step is accepted when the RMS of
                      ``err / (atol + rtol * max(|y_n|, |y_n+1|))`` over the
                      ``error_mask`` entries is at most one.
    dt_max, dt_min  : step bounds; a step below ``dt_min`` raises.
    output_interval : store the state exactly at every multiple of this
                      interval (steps are clamped to land on them);
                      otherwise every ``store_every``-th accepted step is
                      stored.  The initial and final states are always stored.
    on_accept       : ``on_accept(step, t, y, h)`` after each accepted step.
    store_transform : ``store_transform(t, y)`` -> object stored instead of ``y``.
    """
    y = np.array(y0, dtype=float)
    t = float(t0)
    if dt <= 0:
        raise ValueError("dt must be positive")
    if dt_max is None:
        dt_max = np.inf
    mask = np.ones(y.size, dtype=bool) if error_mask is None else np.asarray(error_mask)
    if mask.dtype != bool:
        m = np.zeros(y.size, dtype=bool)
        m[mask] = True
        mask = m
    tiny = 1e-10 * max(1.0, abs(t_end - t0))

    def store(t, y):
        times.append(t)
        states.append(y.copy() if store_transform is None else store_transform(t, y))

    times, states = [], []
    store(t, y)
    k = [None] * 7
    k[0] = rhs(t, y)
    nfev = 1
    steps = rejected = 0
    h_min_used, h_max_used = np.inf, 0.0
    next_out = t0 + output_interval if output_interval else None
    h = min(dt, dt_max)

    while t < t_end - tiny:
        target = t_end if next_out is None else min(t_end, next_out)
        h = min(h, dt_max)
        clamped = t + h >= target - tiny
        if clamped:
            h = target - t
        if h < dt_min:
            raise RuntimeError(f"RK45 step size {h:.3e} fell below dt_min = {dt_min:.3e} "
                               f"at t = {t:.6g}")

        for i in range(1, 7):
            yi = y + h * sum(a * k[j] for j, a in enumerate(DP_A[i]) if a != 0.0)
            k[i] = rhs(t + DP_C[i] * h, yi)
        nfev += 6
        y_new = y + h * sum(b * k[i] for i, b in enumerate(DP_B) if b != 0.0)
        err_vec = h * sum(e * k[i] for i, e in enumerate(DP_E) if e != 0.0)
        scale = atol + rtol * np.maximum(np.abs(y[mask]), np.abs(y_new[mask]))
        err = float(np.sqrt(np.mean((err_vec[mask] / scale) ** 2))) if mask.any() else 0.0

        if adaptive and err > 1.0:
            rejected += 1
            h *= max(0.2, 0.9 * err ** (-1.0 / ORDER))
            continue

        # accept
        t = target if clamped else t + h
        y = y_new
        k[0] = k[6]
        steps += 1
        h_min_used = min(h_min_used, h)
        h_max_used = max(h_max_used, h)
        if on_accept is not None:
            on_accept(steps, t, y, h)

        if next_out is not None and abs(t - next_out) <= tiny:
            store(t, y)
            next_out += output_interval
        elif next_out is None and steps % store_every == 0:
            store(t, y)
        if t >= t_end - tiny and (not times or abs(times[-1] - t) > tiny):
            store(t, y)

        if adaptive:
            factor = 5.0 if err == 0.0 else min(5.0, max(0.2, 0.9 * err ** (-1.0 / ORDER)))
            proposed = h * factor
            h = max(dt, proposed) if clamped else proposed
            dt = h
        else:
            h = dt
        if verbose and (steps % 100 == 0):
            print(f"  RK45 step {steps}, t = {t:.6g}, dt = {h:.3e}, rejected {rejected}")

    info = {"scheme": "rk45", "steps": steps, "rejected": rejected, "rhs_evaluations": nfev,
            "dt_min_used": h_min_used, "dt_max_used": h_max_used, "adaptive": adaptive,
            "rtol": rtol, "atol": atol}
    return RKResult(np.asarray(times), states, info)
