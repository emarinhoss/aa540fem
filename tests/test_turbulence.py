"""Spalart-Allmaras model functions, wall distance, SA assembly and the variable-viscosity flow."""

import numpy as np
import pytest

from aa540fem import geometry
from aa540fem.incompressible import FlowAssembler, FlowProblem, solve_flow
from aa540fem.turbulence import SpalartAllmaras, SpalartAllmarasSolver


def test_model_constants_and_limits():
    m = SpalartAllmaras(nu=1e-5)
    assert np.isclose(m.cw1, 0.1355 / 0.41 ** 2 + (1 + 0.622) / (2 / 3))
    chi = np.array([0.0, 1.0, 7.1, 100.0])
    fv1 = m.fv1(chi)
    assert fv1[0] == 0 and np.isclose(fv1[2], 0.5) and fv1[3] > 0.999
    nt = chi * m.nu
    nu_t = m.eddy_viscosity(nt)
    assert nu_t[0] == 0 and np.allclose(nu_t, nt * fv1)
    assert m.eddy_viscosity(np.array([-1e-6]))[0] == 0.0      # negative branch: no eddy viscosity
    # diffusivity: (nu + nu_tilde)/sigma for positive, f_n < 1 for negative values
    assert np.isclose(m.diffusivity(np.array([3e-5]))[0], 4e-5 / (2 / 3))
    assert 0 < m.diffusivity(np.array([-5e-6]))[0] < m.nu / (2 / 3)
    # far from walls the destruction vanishes and production is c_b1 S nu_tilde
    s = m.source(np.array([3e-5]), np.zeros(1), np.zeros(1), np.array([2.0]), np.array([1e6]))
    assert np.isclose(s[0], 0.1355 * 2.0 * 3e-5, rtol=1e-3)
    # zero working variable: zero source (apart from the gradient term)
    assert m.source(np.zeros(1), np.zeros(1), np.zeros(1), np.ones(1), np.ones(1))[0] == 0.0
    assert np.isclose(m.source(np.zeros(1), np.ones(1), np.zeros(1), np.ones(1), np.ones(1))[0],
                      0.622 / (2 / 3))
    # S_tilde stays positive when the f_v2 correction would make it negative
    st = m.modified_vorticity(np.array([1e-4]), np.array([1e-3]), np.array([1e-3]))
    assert st[0] > 0


def test_wall_distance_rectangle_and_curved_wall():
    mesh = geometry(2.0, 1.0, 4, "quad9")
    assert np.allclose(mesh.wall_distance("bottom"), mesh.y)
    assert np.allclose(mesh.wall_distance(["bottom", "left"]), np.minimum(mesh.x, mesh.y))
    pytest.importorskip("meshio")
    import pathlib

    from aa540fem import read_mesh

    meshes = pathlib.Path(__file__).resolve().parent.parent / "examples" / "meshes"
    ann = read_mesh(meshes / "annulus_tri6.msh")
    d = ann.wall_distance("inner")
    r = np.hypot(ann.x, ann.y)
    assert np.abs(d - (r - 1.0)).max() < 2e-3          # straight-segment approximation of the arc


def test_sa_jacobian_matches_finite_differences():
    mesh = geometry(2.0, 1.0, 4, "quad9")
    nu = 1e-3
    sa = SpalartAllmarasSolver(mesh, SpalartAllmaras(nu), ["bottom"])
    sa.set_velocity(1.0 + 0.3 * mesh.y ** 2, 0.1 * mesh.x * mesh.y)
    rng = np.random.default_rng(1)
    for sign in (1.0, -1.0):
        nt = sign * nu * (0.5 + 5 * rng.random(mesh.n_nodes))
        d = rng.random(mesh.n_nodes) - 0.5
        eps = 1e-7 * nu
        for supg, tol in ((False, 1e-6), (True, 5e-3)):    # SUPG tau is frozen in the Jacobian
            R, J = sa.residual_jacobian(nt, supg=supg)
            fd = (sa.residual_jacobian(nt + eps * d, supg=supg)[0]
                  - sa.residual_jacobian(nt - eps * d, supg=supg)[0]) / (2 * eps)
            assert np.linalg.norm(J @ d - fd) < tol * np.linalg.norm(fd)


def test_sa_steady_solve_uniform_shear_layer():
    # channel with a linear velocity profile: nu_tilde relaxes to a smooth field
    # between the wall value 0 and the inflow level; the solve converges
    mesh = geometry(2.0, 1.0, 6, "quad9")
    nu = 1e-3
    sa = SpalartAllmarasSolver(mesh, SpalartAllmaras(nu), ["bottom"])
    sa.set_velocity(mesh.y, np.zeros_like(mesh.y))
    inflow = mesh.bc_nodes["left"]
    wall = mesh.bc_nodes["bottom"]
    fixed = np.unique(np.concatenate([inflow, wall]))
    vals = np.where(np.isin(fixed, wall), 0.0, 3 * nu)
    res = sa.solve(np.full(mesh.n_nodes, 3 * nu), fixed, vals)
    assert res.converged
    nt = res.T
    assert np.all(nt[wall] == 0) and np.all(nt > -1e-12)
    assert nt.max() < 50 * nu and nt.max() > 3 * nu   # production makes it grow, not diverge


def test_variable_viscosity_flow_jacobian_and_poiseuille():
    mesh = geometry(2.0, 1.0, 4, "quad9")
    rng = np.random.default_rng(2)
    mu_t = 1e-2 * (1 + mesh.y) * (1 + 0.5 * mesh.x)
    for stab in (False, True):
        prob = FlowProblem(mesh, mu=1e-3, rho=1.0, stabilisation=stab, eddy_viscosity=mu_t,
                           bc={"left": (1.0, 0.0), "top": (0.0, 0.0), "bottom": (0.0, 0.0),
                               "right": "open"})
        asm = FlowAssembler(prob)
        B = asm.Bx + asm.By
        U = rng.random(asm.space.ndof)
        d = rng.random(asm.space.ndof)

        def residual(V):
            mt = asm.momentum_terms(V, param_state=U)
            return asm.K @ V + mt.N + mt.S - B.T @ V + B @ V

        mt = asm.momentum_terms(U)
        J = asm.K + mt.J_N + mt.J_S - B.T + B
        eps = 1e-6
        fd = (residual(U + eps * d) - residual(U - eps * d)) / (2 * eps)
        assert np.linalg.norm(J @ d - fd) < 1e-7 * np.linalg.norm(fd)
    # constant eddy viscosity: Poiseuille with the total viscosity is exact
    inflow = lambda x, y: 4 * y * (1 - y)
    prob = FlowProblem(mesh, mu=0.5, rho=1.0, stabilisation=True,
                       eddy_viscosity=np.full(mesh.n_nodes, 0.5),
                       bc={"left": (inflow, 0.0), "top": (0.0, 0.0), "bottom": (0.0, 0.0),
                           "right": "open"})
    sol = solve_flow(prob)
    assert np.allclose(sol.u, inflow(mesh.x, mesh.y), atol=1e-9)
    assert np.allclose(sol.p_nodal, 8 * (2 - mesh.x), atol=1e-8)
    with pytest.raises(ValueError):
        FlowAssembler(FlowProblem(mesh, eddy_viscosity=np.ones(3), bc={}))


def test_rans_coupling_turbulent_skin_friction_exceeds_laminar():
    """Two outer iterations of the SA coupling on the (coarse) laminar-plate mesh:
    the eddy viscosity grows and the skin friction rises above Blasius."""
    pytest.importorskip("meshio")
    import pathlib
    import warnings

    from aa540fem import read_mesh
    from aa540fem.turbulence import solve_rans

    meshes = pathlib.Path(__file__).resolve().parent.parent / "examples" / "meshes"
    mesh = read_mesh(meshes / "flat_plate_bl.msh")
    nu = 1e-5
    prob = FlowProblem(mesh, mu=nu, rho=1.0, stabilisation=True,
                       bc={"inlet": (1.0, 0.0), "top": (1.0, 0.0), "symmetry": (None, 0.0),
                           "plate": (0.0, 0.0), "outlet": "open"})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rans = solve_rans(prob, wall_tags=["plate"], max_outer=2)
    assert len(rans.history) == 2 and all(h[4] for h in rans.history)   # SA solves converged
    assert np.all(rans.nu_tilde[mesh.bc_nodes["plate"]] == 0)
    assert rans.nu_t.max() > 5 * nu and rans.nu_t.min() >= 0
    assert np.allclose(rans.distance[mesh.bc_nodes["plate"]], 0)
    tr = rans.flow.wall_traction("plate")
    rex = tr["x"] / nu
    sel = (rex > 2e4) & (rex < 1e5)
    ratio = 2 * tr["tx"][sel] / (0.664 / np.sqrt(rex[sel]))
    assert ratio.mean() > 1.3
