"""Time integration: accuracy, temporal orders, time-dependent data, output."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from aa540fem import DIRICHLET, NEUMANN, Problem, solve_transient
from aa540fem.geometry import SIDES

ALL_DIRICHLET = {s: DIRICHLET for s in SIDES}
LAMBDA = np.pi ** 2


def decaying_mode(elems=8, elem_type="quad9"):
    """T = exp(-pi^2 t) sin(pi x) on [0, 1] x [0, 0.5]; T = 0 left/right, insulated top/bottom."""
    return Problem(a=1.0, b=0.5, elems=elems, elem_type=elem_type,
                   bc_type={"left": DIRICHLET, "right": DIRICHLET},
                   bc_val={"left": 0.0, "right": 0.0},
                   material=lambda x, y: (1.0, 0.0, 0.0, 1.0, 0.0))


def mode_exact(x, t):
    return np.exp(-LAMBDA * t) * np.sin(np.pi * x)


def test_backward_euler_and_crank_nicolson_accuracy():
    p = decaying_mode()
    T0 = lambda x, y: np.sin(np.pi * x)
    be = solve_transient(p, scheme="theta", dt=0.01, t_end=0.1, theta=1.0, T0=T0)
    cn = solve_transient(p, scheme="theta", dt=0.01, t_end=0.1, theta=0.5, T0=T0)
    exact = mode_exact(be.mesh.x, 0.1)
    err_be = np.abs(be.T - exact).max()
    err_cn = np.abs(cn.T - exact).max()
    assert err_cn < 1e-3
    assert err_be < 3e-2
    assert err_cn < err_be / 10
    assert be.times[0] == 0 and np.isclose(be.times[-1], 0.1) and len(be.snapshots) == 11
    assert np.allclose(be.snapshots[0], T0(be.mesh.x, be.mesh.y))


@pytest.mark.parametrize("theta, expected", [(1.0, 1.0), (0.5, 2.0)])
def test_temporal_convergence_order(theta, expected):
    p = decaying_mode()
    T0 = lambda x, y: np.sin(np.pi * x)
    ref = solve_transient(p, scheme="theta",
        dt=1e-3, t_end=0.1, theta=0.5, T0=T0, store_every=100).T
    errs = [np.abs(solve_transient(p, scheme="theta",
        dt=dt, t_end=0.1, theta=theta, T0=T0, store_every=1000).T
                   - ref).max() for dt in (0.02, 0.01, 0.005)]
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(rates > expected - 0.25), (errs, rates)


def test_manufactured_solution_with_convection_and_time_dependent_dirichlet():
    # T = (1 + t) (sin(pi x) sin(pi y) + x), u = (1, 0.5), kappa = I, rho_c = 1
    u = (1.0, 0.5)
    s = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)
    exact = lambda x, y, t: (1.0 + t) * (s(x, y) + x)

    def material(x, y, t):
        conv = (1.0 + t) * (u[0] * (np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) + 1.0)
                            + u[1] * np.pi * np.sin(np.pi * x) * np.cos(np.pi * y))
        f = (s(x, y) + x) + conv + (1.0 + t) * 2 * np.pi ** 2 * s(x, y)
        return 1.0, 0.0, 0.0, 1.0, f

    p = Problem(a=1.0, b=1.0, elems=12, elem_type="quad9", bc_type=ALL_DIRICHLET,
                bc_val={s_: exact for s_ in SIDES}, material=material,
                velocity=lambda x, y: u, supg=False)
    cn = solve_transient(p, scheme="theta",
        dt=0.1, t_end=1.0, theta=0.5, T0=lambda x, y: exact(x, y, 0.0))
    err = np.abs(cn.T - exact(cn.mesh.x, cn.mesh.y, 1.0)).max()
    assert err < 5e-3
    # Crank-Nicolson is exact in time for solutions linear in t: halving dt changes nothing
    cn2 = solve_transient(p, scheme="theta",
        dt=0.05, t_end=1.0, theta=0.5, T0=lambda x, y: exact(x, y, 0.0))
    assert np.allclose(cn.T, cn2.T, atol=1e-9)
    # right side carries the time-dependent Dirichlet value (1 + t) x = 2 at t = 1
    assert np.allclose(cn.T[cn.mesh.bc_nodes["right"]], 2.0)
    # SUPG on: still accurate
    p.supg = True
    supg = solve_transient(p, scheme="theta",
        dt=0.1, t_end=1.0, theta=0.5, T0=lambda x, y: exact(x, y, 0.0))
    assert np.abs(supg.T - exact(supg.mesh.x, supg.mesh.y, 1.0)).max() < 1e-2


def test_time_dependent_velocity_path_matches_constant_one():
    kw = dict(a=1.0, b=1.0, elems=8, elem_type="quad",
              bc_type={"left": DIRICHLET, "right": DIRICHLET},
              bc_val={"left": 0.0, "right": 1.0},
              material=lambda x, y: (0.05, 0.0, 0.0, 0.05, 0.0))
    const = solve_transient(Problem(velocity=lambda x, y: (1.0, 0.0), **kw),
                            dt=0.05, t_end=0.5, theta=0.5)
    timed = solve_transient(Problem(velocity=lambda x, y, t: (1.0 + 0.0 * t, 0.0), **kw),
                            dt=0.05, t_end=0.5, theta=0.5)
    assert np.allclose(const.T, timed.T, atol=1e-10)
    assert const.T.min() > -1e-3 and const.T.max() <= 1.0 + 1e-9


def test_neumann_flux_heats_insulated_bar():
    # Unit square, kappa = I, rho_c = 1, flux q = 1 on the right, insulated elsewhere:
    # the mean temperature rises at rate q * length / area = 1 per unit time.
    p = Problem(a=1.0, b=1.0, elems=6, elem_type="quad",
                bc_type={"right": NEUMANN, "left": DIRICHLET},
                bc_val={"right": 1.0, "left": 0.0},
                material=lambda x, y: (1.0, 0.0, 0.0, 1.0, 0.0))
    sol = solve_transient(p, scheme="theta", dt=0.1, t_end=5.0, theta=1.0)
    # steady state T = x is reached: n.grad T = 1 on the right, T = 0 on the left
    assert np.allclose(sol.T, sol.mesh.x, atol=1e-3)


def test_series_output(tmp_path):
    pytest.importorskip("meshio")
    p = decaying_mode(elems=3, elem_type="quad")
    sol = solve_transient(p, scheme="theta", dt=0.05, t_end=0.2, theta=0.5,
                          T0=lambda x, y: np.sin(np.pi * x), store_every=2)
    assert len(sol.snapshots) == 3 and np.allclose(sol.times, [0, 0.1, 0.2])
    pvd = sol.save_series(tmp_path / "out" / "mode")
    root = ET.parse(pvd).getroot()
    sets = root.findall("./Collection/DataSet")
    assert [float(d.get("timestep")) for d in sets] == [0.0, 0.1, 0.2]
    for d in sets:
        assert (tmp_path / "out" / d.get("file")).exists()
    assert np.allclose(sol.snapshot(0.11), sol.snapshots[1])
    sol.save(tmp_path / "final.vtu")
    assert (tmp_path / "final.vtu").exists()


def test_invalid_time_step_arguments():
    p = decaying_mode(elems=2, elem_type="quad")
    with pytest.raises(ValueError):
        solve_transient(p, scheme="theta", dt=0.3, t_end=1.0)
    with pytest.raises(ValueError):
        solve_transient(p, scheme="theta", dt=0.1, t_end=1.0, theta=1.5)
