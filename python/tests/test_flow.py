"""Incompressible Navier-Stokes with Taylor-Hood elements."""

import numpy as np
import pytest

from aa540fem import error_norms, geometry
from aa540fem.flow import FlowProblem, TaylorHoodSpace, solve_flow, solve_flow_transient

SIDES = ("top", "right", "left", "bottom")


# --------------------------------------------------------------- Poiseuille
def poiseuille(elem_type="quad9", elems=6, length=2.0):
    """Channel of height 1: u = 4 y (1 - y), p = 8 (L - x) for mu = 1."""
    mesh = geometry(length, 1.0, elems, elem_type)
    inflow = lambda x, y: 4 * y * (1 - y)
    bc = {"left": (inflow, 0.0), "bottom": (0.0, 0.0), "top": (0.0, 0.0), "right": "open"}
    return FlowProblem(mesh, mu=1.0, rho=1.0, bc=bc), inflow


@pytest.mark.parametrize("elem_type", ["quad9", "triangle6"])
def test_poiseuille_is_exact(elem_type):
    p, inflow = poiseuille(elem_type)
    sol = solve_flow(p)
    m = sol.mesh
    assert sol.info["converged"] and sol.info["iterations"] <= 5
    assert np.allclose(sol.u, inflow(m.x, m.y), atol=1e-10)
    assert np.allclose(sol.v, 0.0, atol=1e-10)
    assert np.allclose(sol.p_nodal, 8 * (2.0 - m.x), atol=1e-9)
    assert sol.divergence_norm() < 1e-12
    assert not p.pins_pressure


@pytest.mark.parametrize("elem_type", ["quad9", "triangle6"])
def test_forces_on_channel_walls(elem_type):
    # wall shear mu du/dy = 4 over length 2 -> Fx = 8 on each wall; pressure
    # pushes the bottom wall down and the top wall up with int p dx = 16.
    sol = solve_flow(poiseuille(elem_type)[0])
    fx, fy = sol.forces("bottom")
    assert np.isclose(fx, 8.0, atol=1e-9) and np.isclose(fy, -16.0, atol=1e-9)
    fx, fy = sol.forces("top")
    assert np.isclose(fx, 8.0, atol=1e-9) and np.isclose(fy, 16.0, atol=1e-9)
    fx, fy = sol.forces("left")           # inlet: pressure 16 pushes the inlet plane in -x
    assert np.isclose(fx, -16.0, atol=1e-9) and abs(fy) < 1e-9
    with pytest.raises(ValueError):
        sol.forces("nowhere")


def test_stokes_and_transient_reach_poiseuille():
    p, inflow = poiseuille("quad9", 4)
    stokes = solve_flow(p, stokes=True)
    assert stokes.info["stokes"] and stokes.info["iterations"] == 1
    assert np.allclose(stokes.u, inflow(p.mesh.x, p.mesh.y), atol=1e-10)
    run = solve_flow_transient(p, dt=0.1, t_end=3.0, theta=1.0, store_every=10)
    assert np.allclose(run.final.u, inflow(p.mesh.x, p.mesh.y), atol=1e-3)
    assert len(run.snapshots) == 4 and max(run.info["newton_iterations"]) <= 4


# ---------------------------------------------------------------- Kovasznay
def kovasznay(elems, re=40.0):
    lam = re / 2 - np.sqrt(re ** 2 / 4 + 4 * np.pi ** 2)
    ue = lambda x, y: 1 - np.exp(lam * x) * np.cos(2 * np.pi * y)
    ve = lambda x, y: lam / (2 * np.pi) * np.exp(lam * x) * np.sin(2 * np.pi * y)
    pe = lambda x, y: 0.5 * (1 - np.exp(2 * lam * x))
    mesh = geometry(1.5, 1.0, elems, "quad9")
    mesh.points -= [0.5, 0.5]
    prob = FlowProblem(mesh, mu=1 / re, rho=1.0, bc={s: (ue, ve) for s in SIDES}, pin_value=pe)
    return prob, ue, ve, pe


def test_kovasznay_convergence_rates():
    errs = []
    for elems in (8, 16):
        prob, ue, ve, pe = kovasznay(elems)
        assert prob.pins_pressure
        sol = solve_flow(prob)
        assert sol.info["converged"] and sol.info["iterations"] <= 6
        m = sol.mesh
        errs.append((error_norms(m, sol.u, ue)["L2"], error_norms(m, sol.v, ve)["L2"],
                     error_norms(m, sol.p_nodal, pe)["L2"]))
    (u1, v1, p1), (u2, v2, p2) = errs
    assert u2 < 5e-4 and np.log2(u1 / u2) > 2.8       # quadratic velocity: order 3
    assert np.log2(v1 / v2) > 2.8
    assert np.log2(p1 / p2) > 1.5                     # linear pressure: order 2


# ------------------------------------------------------------------- cavity
def test_lid_driven_cavity_matches_ghia():
    mesh = geometry(1.0, 1.0, 32, "quad9")
    walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
    prob = FlowProblem(mesh, mu=0.01, rho=1.0, bc={**walls, "top": (1.0, 0.0)})
    sol = solve_flow(prob)
    assert sol.info["converged"] and sol.info["iterations"] <= 8
    mid = np.isclose(mesh.x, 0.5)
    u, y = sol.u[mid], mesh.y[mid]
    midy = np.isclose(mesh.y, 0.5)
    v = sol.v[midy]
    # Ghia, Ghia & Shin (1982), Re = 100, 129x129 grid
    assert abs(u.min() - (-0.2109)) < 0.013
    assert abs(y[u.argmin()] - 0.4531) < 0.03
    assert abs(v.max() - 0.1753) < 0.01
    assert abs(v.min() - (-0.2453)) < 0.01
    # discrete continuity holds exactly (the pointwise divergence does not: the
    # lid velocity jumps at the corners)
    B = sol.assembler.Bx + sol.assembler.By
    assert np.linalg.norm(B @ sol.U) < 1e-9
    # the pinned pressure node holds its value
    assert sol.p[0] == 0.0


# ---------------------------------------------------------------- transient
def test_taylor_green_vortex_decay():
    nu = 0.01
    decay = lambda t: np.exp(-2 * np.pi ** 2 * nu * t)
    ue = lambda x, y, t: -np.cos(np.pi * x) * np.sin(np.pi * y) * decay(t)
    ve = lambda x, y, t: np.sin(np.pi * x) * np.cos(np.pi * y) * decay(t)
    pe = lambda x, y, t: -0.25 * (np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y)) * decay(t) ** 2
    errs = []
    for elems in (8, 16):
        mesh = geometry(1.0, 1.0, elems, "quad9")
        prob = FlowProblem(mesh, mu=nu, rho=1.0, bc={s: (ue, ve) for s in SIDES},
                           pin_value=lambda x, y: pe(x, y, 0.2))
        assert prob.depends_on_time()
        run = solve_flow_transient(prob, dt=0.02, t_end=0.2, theta=0.5,
                                   U0=lambda x, y: (ue(x, y, 0.0), ve(x, y, 0.0)))
        sol = run.final
        m = sol.mesh
        assert max(run.info["newton_iterations"]) <= 4
        errs.append((np.abs(sol.u - ue(m.x, m.y, 0.2)).max(),
                     np.abs(sol.v - ve(m.x, m.y, 0.2)).max(),
                     np.abs(sol.p_nodal - pe(m.x, m.y, 0.2)).max()))
    (u1, v1, p1), (u2, v2, p2) = errs
    assert u2 < 1e-3 and v2 < 1e-3 and p2 < 2e-2
    assert u1 / u2 > 6 and v1 / v2 > 6 and p1 / p2 > 2.5    # spatial convergence dominates


def test_series_output_and_validation(tmp_path):
    pytest.importorskip("meshio")
    p, _ = poiseuille("quad9", 2)
    run = solve_flow_transient(p, dt=0.5, t_end=1.0, theta=1.0)
    pvd = run.save_series(tmp_path / "flow")
    assert pvd.exists() and (tmp_path / "flow_0002.vtu").exists()
    run.final.save(tmp_path / "final.vtu")
    assert (tmp_path / "final.vtu").exists()

    with pytest.raises(ValueError, match="quadratic"):
        FlowProblem(geometry(1.0, 1.0, 2, "quad"), bc={}).validate()
    with pytest.raises(ValueError, match="Unknown boundary"):
        FlowProblem(geometry(1.0, 1.0, 2, "quad9"), bc={"lid": (1.0, 0.0)}).validate()
    with pytest.raises(ValueError, match="ux, uy"):
        FlowProblem(geometry(1.0, 1.0, 2, "quad9"), bc={"top": 1.0}).validate()
    with pytest.raises(ValueError):
        solve_flow_transient(p, dt=0.3, t_end=1.0)


def test_space_layout():
    mesh = geometry(1.0, 1.0, 3, "quad9")
    space = TaylorHoodSpace(mesh)
    assert space.N == 49 and space.Np == 16 and space.ndof == 2 * 49 + 16
    assert np.array_equal(space.pressure_nodes, np.unique(mesh.conn[:, :4]))
    p = np.arange(space.Np, dtype=float)
    pn = space.pressure_at_nodes(p)
    assert np.allclose(pn[space.pressure_nodes], p)
    # the centre node of the first element gets the mean of its four corners
    c = mesh.conn[0]
    assert np.isclose(pn[c[8]], pn[c[:4]].mean())
    with pytest.raises(ValueError):
        space.dof_p([c[8]])


def test_cylinder_benchmark_drag_and_pressure_drop():
    """Schaefer-Turek 2D-1 (Re = 20): C_D = 5.5795, dp = 0.1175, C_L = 0.0106."""
    pytest.importorskip("meshio")
    import pathlib
    import sys

    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "examples"))
    import cylinder

    sol, res = cylinder.run(cylinder.HERE / "cylinder_tri6.msh", verbose=False)
    assert sol.info["converged"] and sol.info["iterations"] <= 8
    assert abs(res["C_D"] - 5.5795) / 5.5795 < 0.01
    assert abs(res["dp"] - 0.1175) / 0.1175 < 0.01
    assert abs(res["C_L"]) < 0.02          # lift is tiny and mesh sensitive
