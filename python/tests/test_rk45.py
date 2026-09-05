"""Dormand-Prince RK45: the integrator itself, the scalar solver and the flow solver."""

import numpy as np
import pytest

from aa540fem import DIRICHLET, NEUMANN, Problem, geometry, rk45, solve, solve_transient
from aa540fem.flow import FlowProblem, solve_flow_transient
from aa540fem.geometry import SIDES

ALL_DIRICHLET = {s: DIRICHLET for s in SIDES}


# ------------------------------------------------------------- integrator
def test_fixed_step_is_fifth_order():
    exact = np.exp(1.0)
    errs = []
    for dt in (0.1, 0.05):
        res = rk45(lambda t, y: y, np.array([1.0]), 0.0, 1.0, dt, adaptive=False)
        errs.append(abs(res.y[0] - exact))
        assert np.isclose(res.times[-1], 1.0)
        assert res.info["steps"] == round(1.0 / dt)
    assert errs[0] / errs[1] > 28


def test_adaptive_respects_tolerance_and_counts_rejections():
    # y' = -50 (y - cos t): the controller must shrink the step below the
    # explicit stability limit and reject the far-too-large initial step
    rhs = lambda t, y: -50.0 * (y - np.cos(t))
    res = rk45(rhs, np.array([0.0]), 0.0, 2.0, 1.0, rtol=1e-6, atol=1e-9)
    exact = (50 * (50 * np.cos(2.0) + np.sin(2.0)) - 2500 * np.exp(-100.0)) / 2501
    assert abs(res.y[0] - exact) < 1e-5
    assert res.info["rejected"] >= 1
    assert res.info["dt_max_used"] < 0.1 and res.info["steps"] > 30


def test_output_lands_on_fixed_times_and_callbacks_run():
    seen = []
    res = rk45(lambda t, y: np.cos(t) * np.ones_like(y), np.zeros(3), 0.0, 1.0, 0.07,
               output_interval=0.25, on_accept=lambda n, t, y, h: seen.append((n, t, h)))
    assert np.allclose(res.times, [0, 0.25, 0.5, 0.75, 1.0])
    assert np.allclose([s[0] for s in res.states], np.sin(res.times), atol=1e-8)
    assert len(seen) == res.info["steps"] and seen[-1][1] == pytest.approx(1.0)
    # store_every without an output interval
    res = rk45(lambda t, y: -y, np.ones(1), 0.0, 1.0, 0.1, adaptive=False, store_every=4)
    assert np.allclose(res.times, [0, 0.4, 0.8, 1.0])


def test_error_mask_and_transform_and_dt_min():
    mask = np.array([True, False])
    res = rk45(lambda t, y: np.array([-y[0], 1e6 * np.cos(1e6 * t)]), np.array([1.0, 0.0]),
               0.0, 0.5, 0.1, error_mask=mask, store_transform=lambda t, y: y[0] * 2)
    assert res.states[-1] == pytest.approx(2 * np.exp(-0.5), rel=1e-4)
    with pytest.raises(RuntimeError, match="dt_min"):
        rk45(lambda t, y: -1e8 * y, np.ones(1), 0.0, 1.0, 0.1, dt_min=1e-3)
    with pytest.raises(ValueError):
        solve_transient(Problem(elems=2, elem_type="quad"), dt=0.1, t_end=0.2, scheme="bogus")


# ------------------------------------------------------------- scalar solver
def decaying_mode(elems=4, elem_type="quad"):
    return Problem(a=1.0, b=0.5, elems=elems, elem_type=elem_type,
                   bc_type={"left": DIRICHLET, "right": DIRICHLET},
                   bc_val={"left": 0.0, "right": 0.0},
                   material=lambda x, y: (1.0, 0.0, 0.0, 1.0, 0.0))


def test_scalar_rk45_temporal_order_and_adaptive_accuracy():
    p = decaying_mode()
    T0 = lambda x, y: np.sin(np.pi * x)
    ref = solve_transient(p, dt=1e-4, t_end=0.1, adaptive=False, store_every=10 ** 9, T0=T0).T
    errs = [np.abs(solve_transient(p, dt=dt, t_end=0.1, adaptive=False, store_every=10 ** 9,
                                   T0=T0).T - ref).max() for dt in (0.01, 0.005)]
    assert np.log2(errs[0] / errs[1]) > 4.5
    run = solve_transient(p, dt=0.05, t_end=0.1, T0=T0, rtol=1e-8, atol=1e-10)
    assert run.info["scheme"] == "rk45" and run.info["steps"] > 5
    assert np.abs(run.T - ref).max() < 1e-7
    # the exact continuous solution is met to the spatial accuracy of the mesh
    assert np.abs(run.T - np.exp(-np.pi ** 2 * 0.1) * np.sin(np.pi * run.mesh.x)).max() < 3e-2


def test_scalar_rk45_manufactured_with_convection_and_time_dependent_dirichlet():
    u = (1.0, 0.5)
    s = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)
    exact = lambda x, y, t: (1.0 + t) * (s(x, y) + x)

    def material(x, y, t):
        conv = (1.0 + t) * (u[0] * (np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) + 1.0)
                            + u[1] * np.pi * np.sin(np.pi * x) * np.cos(np.pi * y))
        f = (s(x, y) + x) + conv + (1.0 + t) * 2 * np.pi ** 2 * s(x, y)
        return 1.0, 0.0, 0.0, 1.0, f

    p = Problem(a=1.0, b=1.0, elems=4, elem_type="quad9", bc_type=ALL_DIRICHLET,
                bc_val={s_: exact for s_ in SIDES}, material=material,
                velocity=lambda x, y: u, supg=False)
    run = solve_transient(p, dt=0.01, t_end=1.0, T0=lambda x, y: exact(x, y, 0.0),
                          output_interval=0.25)
    assert np.allclose(run.times, [0, 0.25, 0.5, 0.75, 1.0])
    assert np.abs(run.T - exact(run.mesh.x, run.mesh.y, 1.0)).max() < 2e-2
    assert np.allclose(run.T[run.mesh.bc_nodes["right"]], 2.0)


def test_scalar_rk45_nonlinear_without_newton():
    p = Problem(a=1.0, b=1.0, elems=4, elem_type="quad9",
                bc_type={"left": DIRICHLET, "right": DIRICHLET},
                bc_val={"left": 0.0, "right": 1.0},
                material=lambda x, y, T: (1.0 + T, 0.0, 0.0, 1.0 + T, 0.0))
    run = solve_transient(p, dt=0.01, t_end=1.5, T0=0.0, store_every=10 ** 9)
    steady = solve(p)
    assert "newton_iterations" not in run.info and run.info["scheme"] == "rk45"
    assert np.allclose(run.T, steady.T, atol=2e-3)


def test_scalar_rk45_neumann_heating():
    p = Problem(a=1.0, b=1.0, elems=4, elem_type="quad",
                bc_type={"right": NEUMANN, "left": DIRICHLET},
                bc_val={"right": 1.0, "left": 0.0},
                material=lambda x, y: (1.0, 0.0, 0.0, 1.0, 0.0))
    run = solve_transient(p, dt=0.01, t_end=4.0, store_every=10 ** 9)
    assert np.allclose(run.T, run.mesh.x, atol=2e-3)


# --------------------------------------------------------------- flow solver
def test_flow_rk45_poiseuille_start_up_is_divergence_free_and_exact():
    mesh = geometry(2.0, 1.0, 4, "quad9")
    inflow = lambda x, y: 4 * y * (1 - y)
    prob = FlowProblem(mesh, mu=1.0, rho=1.0,
                       bc={"left": (inflow, 0.0), "bottom": (0.0, 0.0), "top": (0.0, 0.0),
                           "right": "open"})
    run = solve_flow_transient(prob, dt=0.005, t_end=2.0, output_interval=0.5)
    assert run.info["scheme"] == "rk45" and np.allclose(run.times, [0, 0.5, 1.0, 1.5, 2.0])
    B = run.assembler.Bx + run.assembler.By
    for U in run.snapshots:
        assert np.linalg.norm(B @ U) < 1e-9
    sol = run.final
    assert np.allclose(sol.u, inflow(mesh.x, mesh.y), atol=2e-3)
    assert np.allclose(sol.p_nodal, 8 * (2.0 - mesh.x), atol=5e-2)


def test_flow_rk45_taylor_green_matches_theta():
    nu = 0.01
    decay = lambda t: np.exp(-2 * np.pi ** 2 * nu * t)
    ue = lambda x, y, t: -np.cos(np.pi * x) * np.sin(np.pi * y) * decay(t)
    ve = lambda x, y, t: np.sin(np.pi * x) * np.cos(np.pi * y) * decay(t)
    pe = lambda x, y, t: -0.25 * (np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y)) * decay(t) ** 2
    mesh = geometry(1.0, 1.0, 8, "quad9")
    prob = FlowProblem(mesh, mu=nu, rho=1.0, bc={s: (ue, ve) for s in SIDES},
                       pin_value=lambda x, y: pe(x, y, 0.2))
    U0 = lambda x, y: (ue(x, y, 0.0), ve(x, y, 0.0))
    rk = solve_flow_transient(prob, dt=0.01, t_end=0.2, U0=U0, rtol=1e-6, atol=1e-8)
    th = solve_flow_transient(prob, dt=0.005, t_end=0.2, theta=0.5, U0=U0, scheme="theta")
    m = mesh
    err_rk = np.abs(rk.final.u - ue(m.x, m.y, 0.2)).max()
    err_th = np.abs(th.final.u - ue(m.x, m.y, 0.2)).max()
    assert err_rk < 2e-2 and abs(err_rk - err_th) < 2e-3        # both spatially limited
    assert np.abs(rk.final.p_nodal - pe(m.x, m.y, 0.2)).max() < 0.1
    assert np.allclose(rk.final.u, th.final.u, atol=2e-3)


def test_flow_rk45_impulsive_start_is_projected_and_logs_forces():
    mesh = geometry(1.0, 1.0, 6, "quad9")
    walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
    prob = FlowProblem(mesh, mu=0.05, rho=1.0, bc={**walls, "top": (1.0, 0.0)})
    steps = []
    run = solve_flow_transient(prob, dt=0.01, t_end=0.3, U0=lambda x, y: (0.0, 0.0),
                               callback=lambda n, t, sol: steps.append((n, t, sol.forces("top"))))
    B = run.assembler.Bx + run.assembler.By
    assert np.linalg.norm(B @ run.snapshots[0]) < 1e-9
    assert len(steps) == run.info["steps"] and steps[-1][1] == pytest.approx(0.3)
    assert run.final.speed.max() <= 1.0 + 1e-6 and run.final.p[0] == 0.0
    assert steps[-1][2][0] < 0                  # the fluid pulls the moving lid backwards
