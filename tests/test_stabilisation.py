"""SUPG / grad-div / PSPG stabilisation of the flow solver and the pseudo-transient continuation."""

import warnings

import numpy as np
import pytest

from aa540fem import error_norms, geometry, read_mesh
from aa540fem.incompressible import FlowAssembler, FlowProblem, solve_flow, solve_flow_transient

SIDES = ("top", "right", "left", "bottom")
MESHES = __import__("pathlib").Path(__file__).resolve().parent.parent / "examples" / "meshes"


def poiseuille(elem_type, **kw):
    mesh = geometry(2.0, 1.0, 6, elem_type)
    inflow = lambda x, y: 4 * y * (1 - y)
    bc = {"left": (inflow, 0.0), "bottom": (0.0, 0.0), "top": (0.0, 0.0), "right": "open"}
    return FlowProblem(mesh, mu=1.0, rho=1.0, bc=bc, **kw), inflow


@pytest.mark.parametrize("elem_type", ["quad9", "triangle6"])
@pytest.mark.parametrize("pspg", [False, True])
def test_stabilised_poiseuille_stays_exact(elem_type, pspg):
    # residual consistency: the exact solution has zero strong residual
    # (including the viscous term, which needs the shape-function Hessians)
    p, inflow = poiseuille(elem_type, stabilisation=True, pspg=pspg)
    sol = solve_flow(p)
    m = sol.mesh
    assert sol.info["converged"] and sol.info["stabilisation"]
    assert np.allclose(sol.u, inflow(m.x, m.y), atol=1e-9)
    assert np.allclose(sol.p_nodal, 8 * (2.0 - m.x), atol=1e-8)
    fx, fy = sol.forces("bottom")
    assert np.isclose(fx, 8.0, atol=1e-8) and np.isclose(fy, -16.0, atol=1e-8)


def test_jacobian_matches_finite_differences_with_frozen_parameters():
    mesh = geometry(2.0, 1.0, 4, "quad9")
    rng = np.random.default_rng(0)
    for pspg in (False, True):
        prob = FlowProblem(mesh, mu=1e-3, rho=1.0, stabilisation=True, pspg=pspg,
                           bc={"left": (1.0, 0.0), "top": (0.0, 0.0), "bottom": (0.0, 0.0),
                               "right": "open"})
        asm = FlowAssembler(prob)
        B = asm.Bx + asm.By
        U = rng.random(asm.space.ndof)
        d = rng.random(asm.space.ndof)

        def residual(V):
            mt = asm.momentum_terms(V, param_state=U)
            return asm.K @ V + mt.N + mt.S - B.T @ V + B @ V

        mt = asm.momentum_terms(U, param_state=U)      # parameters frozen at U
        J = asm.K + mt.J_N + mt.J_S - B.T + B
        eps = 1e-6
        fd = (residual(U + eps * d) - residual(U - eps * d)) / (2 * eps)
        assert np.linalg.norm(J @ d - fd) < 1e-7 * np.linalg.norm(fd)
        assert mt.J_S.nnz > 0 and np.abs(mt.S).max() > 0


def test_full_jacobian_with_parameter_derivatives_on_smooth_field():
    # tau, gamma and the flow-direction element length are differentiated; on a
    # smooth field (no kinks of |s . grad phi|) the Jacobian is consistent
    mesh = geometry(2.0, 1.0, 6, "quad9")
    inflow = lambda x, y: 4 * y * (1 - y)
    for mu, pspg in ((1e-2, False), (1e-4, False), (1e-4, True)):
        prob = FlowProblem(mesh, mu=mu, rho=1.0, stabilisation=True, pspg=pspg,
                           bc={"left": (inflow, 0.0), "top": (0.0, 0.0), "bottom": (0.0, 0.0),
                               "right": "open"})
        asm = FlowAssembler(prob)
        B = asm.Bx + asm.By
        N = asm.space.N
        U = np.zeros(asm.space.ndof)
        U[:N] = inflow(mesh.x, mesh.y) * (1 + 0.2 * np.sin(np.pi * mesh.x))
        U[N:2 * N] = 0.1 * np.sin(np.pi * mesh.x) * np.sin(np.pi * mesh.y)
        d = np.zeros_like(U)
        d[:N] = np.cos(mesh.x) * mesh.y
        d[N:2 * N] = np.sin(mesh.y) * mesh.x
        d[2 * N:] = np.cos(mesh.x[asm.space.pressure_nodes])

        def residual(V):
            mt = asm.momentum_terms(V)
            return asm.K @ V + mt.N + mt.S - B.T @ V + B @ V

        mt = asm.momentum_terms(U)
        J = asm.K + mt.J_N + mt.J_S - B.T + B
        eps = 1e-6
        fd = (residual(U + eps * d) - residual(U - eps * d)) / (2 * eps)
        assert np.linalg.norm(J @ d - fd) < 1e-7 * np.linalg.norm(fd)


def test_stabilised_kovasznay_keeps_convergence_order():
    re = 40.0
    lam = re / 2 - np.sqrt(re ** 2 / 4 + 4 * np.pi ** 2)
    ue = lambda x, y: 1 - np.exp(lam * x) * np.cos(2 * np.pi * y)
    ve = lambda x, y: lam / (2 * np.pi) * np.exp(lam * x) * np.sin(2 * np.pi * y)
    pe = lambda x, y: 0.5 * (1 - np.exp(2 * lam * x))
    errs = []
    for elems in (8, 16):
        mesh = geometry(1.5, 1.0, elems, "quad9")
        mesh.points -= [0.5, 0.5]
        sol = solve_flow(FlowProblem(mesh, mu=1 / re, rho=1.0, bc={s: (ue, ve) for s in SIDES},
                                     pin_value=pe, stabilisation=True))
        assert sol.info["converged"]
        errs.append(error_norms(mesh, sol.u, ue)["L2"])
    assert np.log2(errs[0] / errs[1]) > 2.5 and errs[1] < 5e-4


def cavity(mu, elems=32, **kw):
    mesh = geometry(1.0, 1.0, elems, "quad9")
    walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
    return FlowProblem(mesh, mu=mu, rho=1.0, bc={**walls, "top": (1.0, 0.0)}, **kw)


def test_cavity_re1000_needs_stabilisation():
    # Galerkin Newton fails on this mesh; stabilised Newton converges and the
    # centreline extrema approach Ghia et al. (129x129 reference; the 32x32 Q2
    # mesh is coarse, so a 15 % band)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        galerkin = solve_flow(cavity(1e-3), continuation="newton", max_newton=15)
    assert not galerkin.info["converged"]
    sol = solve_flow(cavity(1e-3, stabilisation=True))
    assert sol.info["converged"] and sol.info["continuation"] == "newton"
    m = sol.mesh
    u, y = sol.u[np.isclose(m.x, 0.5)], m.y[np.isclose(m.x, 0.5)]
    v = sol.v[np.isclose(m.y, 0.5)]
    assert abs(u.min() / -0.3829 - 1) < 0.15 and abs(y[u.argmin()] - 0.1719) < 0.05
    assert abs(v.max() / 0.3709 - 1) < 0.15
    assert abs(v.min() / -0.5155 - 1) < 0.15


def test_ptc_matches_newton_and_auto_falls_back():
    p = cavity(0.01, elems=12, stabilisation=True)
    newton = solve_flow(p, continuation="newton")
    ptc = solve_flow(p, continuation="ptc")
    assert ptc.info["converged"] and ptc.info["continuation"] == "ptc" and "dtau" in ptc.info
    assert np.allclose(ptc.U, newton.U, atol=1e-7)
    # auto mode: when Newton does not converge (here: not allowed enough
    # iterations) the solve continues with PTC from the initial state
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = solve_flow(p, continuation="auto", max_newton=2)
    assert sol.info["converged"] and sol.info["continuation"] == "newton+ptc"
    assert np.allclose(sol.U, newton.U, atol=1e-7)
    with pytest.raises(ValueError):
        solve_flow(p, continuation="bogus")


def test_mixed_boundary_layer_mesh_cylinder():
    pytest.importorskip("meshio")
    mesh = read_mesh(MESHES / "cylinder_bl.msh")
    assert set(mesh.cells) == {"triangle6", "quad9"}
    assert all(d.min() > 0 for d in mesh.jacobian_at_centroids().values())
    assert set(mesh.tags) == {"inlet", "outlet", "walls", "cylinder"}
    H = 0.41
    inflow = lambda x, y: 4.0 * 0.3 * y * (H - y) / H ** 2
    prob = FlowProblem(mesh, mu=1e-3, rho=1.0, stabilisation=True,
                       bc={"inlet": (inflow, 0.0), "walls": (0.0, 0.0),
                           "cylinder": (0.0, 0.0), "outlet": "open"})
    sol = solve_flow(prob)
    fx, fy = sol.forces("cylinder")
    assert abs(500 * fx - 5.5795) / 5.5795 < 0.005      # boundary layer resolves C_D and C_L
    assert abs(500 * fy - 0.0106) < 0.002
    tr = sol.wall_traction("cylinder")
    assert np.isclose(np.sum(tr["weight"] * tr["tx"]), fx)
    assert np.isclose(np.sum(tr["weight"]), 2 * np.pi * 0.05, rtol=1e-3)


def test_wall_traction_on_poiseuille():
    sol = solve_flow(poiseuille("quad9")[0])
    tr = sol.wall_traction("bottom")
    assert np.allclose(tr["tx"], 4.0, atol=1e-8)        # mu du/dy = 4 along the wall
    assert np.allclose(tr["ny"], -1.0) and np.allclose(tr["nx"], 0.0)
    assert np.all(np.diff(tr["x"]) >= 0)
    assert np.isclose(np.sum(tr["weight"] * tr["ty"]), -16.0)


def test_stabilised_time_stepping_both_schemes():
    p, inflow = poiseuille("quad9", stabilisation=True)
    m = p.mesh
    th = solve_flow_transient(p, dt=0.1, t_end=1.0, theta=1.0, scheme="theta", store_every=10)
    assert np.allclose(th.final.u, inflow(m.x, m.y), atol=2e-2)
    rk = solve_flow_transient(p, dt=0.002, t_end=0.2, output_interval=0.1)
    assert rk.info["steps"] > 5 and np.isfinite(rk.final.u).all()
    B = rk.assembler.Bx + rk.assembler.By
    assert np.linalg.norm(B @ rk.final.U) < 1e-9
