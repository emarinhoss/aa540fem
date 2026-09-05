"""Temperature-dependent coefficients: signature detection, Newton, Picard, transient."""

import numpy as np
import pytest

from aa540fem import DIRICHLET, Problem, assemble, error_norms, geometry, solve, solve_transient
from aa540fem.core.mesh import SIDES
from aa540fem.core.util import accepts_temperature, accepts_time, call_coeff

ALL_DIRICHLET = {s: DIRICHLET for s in SIDES}


# --------------------------------------------------------- argument matching
def test_signature_detection_and_dispatch():
    by_T = lambda x, y, T: T
    by_t = lambda x, y, t: t
    both = lambda x, y, T, t: (T, t)
    swapped = lambda x, y, t, T: (t, T)
    optional = lambda x, y, eps=0.1: eps
    legacy = lambda x, y, z: z

    assert accepts_temperature(by_T) and not accepts_time(by_T)
    assert accepts_time(by_t) and not accepts_temperature(by_t)
    assert accepts_time(both) and accepts_temperature(both)
    assert not accepts_time(optional) and not accepts_temperature(optional)
    assert accepts_time(legacy) and not accepts_temperature(legacy)

    assert call_coeff(both, 1, 2, t=3, T=4) == (4, 3)
    assert call_coeff(swapped, 1, 2, t=3, T=4) == (3, 4)
    assert call_coeff(optional, 1, 2, t=5) == 0.1
    assert call_coeff(legacy, 1, 2, t=5) == 5
    assert call_coeff(by_t, 1, 2) == 0.0

    def kw_only(x, y, *, T):
        return T

    assert call_coeff(kw_only, 1, 2, T=7) == 7
    with pytest.raises(ValueError, match="temperature"):
        call_coeff(by_T, 1, 2)


def test_assemble_needs_T_for_nonlinear_material():
    mesh = geometry(1.0, 1.0, 2, "quad")
    with pytest.raises(ValueError, match="temperature"):
        assemble(mesh, material=lambda x, y, T: (1 + T, 0.0, 0.0, 1 + T, 0.0))


# -------------------------------------------------------- Kirchhoff problem
def kirchhoff(elem_type, elems):
    """kappa = 1 + T, T(0) = 0, T(1) = 1  ->  T = -1 + sqrt(1 + 3 x)."""
    return Problem(a=1.0, b=1.0, elems=elems, elem_type=elem_type,
                   bc_type={"left": DIRICHLET, "right": DIRICHLET},
                   bc_val={"left": 0.0, "right": 1.0},
                   material=lambda x, y, T: (1.0 + T, 0.0, 0.0, 1.0 + T, 0.0))


def kirchhoff_exact(x, y):
    return -1.0 + np.sqrt(1.0 + 3.0 * x)


@pytest.mark.parametrize("elem_type, tol", [("quad", 3e-3), ("quad9", 1e-4), ("triangle6", 1e-4)])
def test_newton_solves_kirchhoff_problem(elem_type, tol):
    sol = solve(kirchhoff(elem_type, 8))
    assert sol.info["nonlinear"] == "newton" and sol.info["converged"]
    assert sol.info["iterations"] <= 8
    err = error_norms(sol.mesh, sol.T, kirchhoff_exact)["L2"]
    assert err < tol
    # quadratic convergence: the last reductions are large
    r = sol.info["residuals"]
    assert r[-1] < 1e-9 * r[0]
    assert r[-1] < 1e-3 * r[-2]


def test_newton_converges_faster_than_picard():
    p = kirchhoff("quad9", 6)
    newton = solve(p)
    picard = solve(p, newton=False, max_newton=200)
    assert picard.info["nonlinear"] == "picard" and picard.info["converged"]
    assert newton.info["iterations"] < picard.info["iterations"]
    assert np.allclose(newton.T, picard.T, atol=1e-8)


def test_initial_guess_and_options_are_honoured():
    p = kirchhoff("quad", 6)
    close = solve(p, T0=kirchhoff_exact)
    far = solve(p, T0=5.0)
    assert close.info["iterations"] <= far.info["iterations"]
    assert np.allclose(close.T, far.T, atol=1e-8)
    with pytest.raises(TypeError):
        solve(Problem(elems=2, elem_type="quad"), T0=1.0)     # linear problem, no Newton args


def test_newton_with_gmres_matches_direct_and_cg_is_refused():
    p = kirchhoff("quad9", 4)
    direct = solve(p)
    gm = solve(p, method="gmres", tol=1e-13)
    assert np.allclose(gm.T, direct.T, atol=1e-8)
    with pytest.raises(ValueError, match="gmres"):
        solve(p, method="cg")


# ----------------------------------------------- manufactured nonlinear case
def _sin(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y)


def _grad2(x, y):
    return np.pi ** 2 * (np.cos(np.pi * x) ** 2 * np.sin(np.pi * y) ** 2
                         + np.sin(np.pi * x) ** 2 * np.cos(np.pi * y) ** 2)


def test_manufactured_kappa_and_source_depending_on_T():
    # -div((1 + T^2) grad T) + T^3 = g with T = s(x, y) = sin(pi x) sin(pi y)
    def material(x, y, T):
        s = _sin(x, y)
        lap = -2 * np.pi ** 2 * s
        g = -((1 + s ** 2) * lap + 2 * s * _grad2(x, y)) + s ** 3
        return 1.0 + T ** 2, 0.0, 0.0, 1.0 + T ** 2, g - T ** 3

    p = Problem(a=1.0, b=1.0, elems=12, elem_type="quad9", bc_type=ALL_DIRICHLET,
                bc_val={s: 0.0 for s in SIDES}, material=material)
    sol = solve(p)
    assert sol.info["converged"] and sol.info["iterations"] <= 8
    assert error_norms(sol.mesh, sol.T, _sin)["L2"] < 2e-3
    assert np.abs(sol.T - _sin(sol.mesh.x, sol.mesh.y)).max() < 5e-3


def test_nonlinear_with_convection_and_supg():
    # -div((1 + T) grad T) + u.grad T = g, T = s, u = (1, 0.5)
    u = (1.0, 0.5)

    def material(x, y, T):
        s = _sin(x, y)
        lap = -2 * np.pi ** 2 * s
        conv = u[0] * np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) \
            + u[1] * np.pi * np.sin(np.pi * x) * np.cos(np.pi * y)
        g = -((1 + s) * lap + _grad2(x, y)) + conv
        return 1.0 + T, 0.0, 0.0, 1.0 + T, g

    p = Problem(a=1.0, b=1.0, elems=12, elem_type="quad9", bc_type=ALL_DIRICHLET,
                bc_val={s: 0.0 for s in SIDES}, material=material, velocity=lambda x, y: u)
    sol = solve(p)
    assert sol.info["converged"] and sol.info["iterations"] <= 10
    assert error_norms(sol.mesh, sol.T, _sin)["L2"] < 5e-3


# ------------------------------------------------------------------ transient
def test_transient_nonlinear_manufactured_solution():
    # T = (1 + t) s, kappa = 1 + T, rho_c = 1:  f = s - div((1 + T) grad T)
    exact = lambda x, y, t: (1.0 + t) * _sin(x, y)

    def material(x, y, T, t):
        s = _sin(x, y)
        a = 1.0 + t
        lap = -2 * np.pi ** 2 * a * s
        f = s - ((1.0 + a * s) * lap + a ** 2 * _grad2(x, y))
        return 1.0 + T, 0.0, 0.0, 1.0 + T, f

    p = Problem(a=1.0, b=1.0, elems=12, elem_type="quad9", bc_type=ALL_DIRICHLET,
                bc_val={s: 0.0 for s in SIDES}, material=material)
    sol = solve_transient(p, scheme="theta", dt=0.1, t_end=1.0, theta=0.5, T0=_sin)
    assert sol.info["nonlinear"] == "newton"
    assert max(sol.info["newton_iterations"]) <= 6
    err = np.abs(sol.T - exact(sol.mesh.x, sol.mesh.y, 1.0)).max()
    assert err < 5e-3


def test_transient_nonlinear_reaches_steady_kirchhoff_state():
    p = kirchhoff("quad9", 6)
    sol = solve_transient(p, scheme="theta", dt=0.25, t_end=5.0, theta=1.0, T0=0.0)
    steady = solve(p)
    assert np.allclose(sol.T, steady.T, atol=2e-3)
    assert len(sol.info["newton_iterations"]) == 20
