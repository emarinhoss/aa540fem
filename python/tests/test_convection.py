"""Convection-diffusion: operators, SUPG stabilisation and non-symmetric solvers."""

import numpy as np
import pytest

from aa540fem import (
    DIRICHLET,
    ELEMENTS,
    Problem,
    assemble_operators,
    elem_eqn,
    elem_operators,
    error_norms,
    geometry,
    get_element,
    solve,
    supg_tau,
)
from aa540fem.geometry import SIDES

EPS = 0.01                       # conductivity -> cell Peclet 5 on 10 cells


def boundary_layer(elem_type, elems=10, supg=True, eps=EPS):
    """u = (1, 0), T(0) = 0, T(1) = 1: exact T = (e^{x/eps} - 1) / (e^{1/eps} - 1)."""
    return Problem(a=1.0, b=1.0, elems=elems, elem_type=elem_type,
                   bc_type={"left": DIRICHLET, "right": DIRICHLET},
                   bc_val={"left": 0.0, "right": 1.0},
                   material=lambda x, y: (eps, 0.0, 0.0, eps, 0.0),
                   velocity=lambda x, y: (1.0, 0.0), supg=supg)


def exact_layer(x, y, eps=EPS):
    return np.expm1(x / eps) / np.expm1(1.0 / eps)


# ------------------------------------------------------------------ operators
@pytest.mark.parametrize("name", sorted(ELEMENTS))
def test_mass_matrix_integrates_capacity(name):
    mesh = geometry(2.0, 3.0, 4, name)
    ops = assemble_operators(mesh, Problem(mesh=mesh, rho_c=2.5))
    assert np.isclose(ops.M.sum(), 2.5 * 6.0)
    assert np.allclose((ops.M - ops.M.T).data, 0)
    assert ops.C.nnz == 0


@pytest.mark.parametrize("name", sorted(ELEMENTS))
@pytest.mark.parametrize("supg", [False, True])
def test_convection_rows_sum_to_zero(name, supg):
    # u . grad(sum_j phi_j) = 0, with or without the streamline weighting
    mesh = geometry(2.0, 3.0, 4, name)
    ops = assemble_operators(mesh, Problem(mesh=mesh, velocity=lambda x, y: (1.0, -0.5),
                                           supg=supg))
    assert np.allclose(np.asarray(ops.C.sum(axis=1)).ravel(), 0, atol=1e-12)
    assert ops.C.nnz > 0


@pytest.mark.parametrize("name", sorted(ELEMENTS))
def test_elem_operators_reproduces_elem_eqn(name):
    mesh = geometry(2.0, 3.0, 3, name)
    el = get_element(name)
    xi, eta, w = el.quadrature()
    phi, dxi, deta = el.shape(xi, eta)
    xe, ye = mesh.x[mesh.conn], mesh.y[mesh.conn]
    material = lambda x, y: (2.0, 0.5, 0.5, 1.0, x + y)
    Ke, fe = elem_eqn(xe, ye, phi, dxi, deta, w, material=material)
    K2, C2, M2, f2 = elem_operators(xe, ye, phi, dxi, deta, w, material=material)
    assert np.allclose(Ke, K2)
    assert np.allclose(fe, f2)
    assert np.allclose(C2, 0)
    assert np.all(M2.sum(axis=(1, 2)) > 0)


def test_supg_tau_limits():
    ux = np.array([[0.0, 1.0, 1.0]])
    uy = np.zeros((1, 3))
    k = np.array([[1.0, 1.0, 0.0]])          # zero conductivity -> Pe = inf
    zero = np.zeros_like(k)
    # gradients of a 1-D linear element of length h = 0.5 (two nodes, |dphi/dx| = 2)
    dphi_dx = np.broadcast_to(np.array([-2.0, 2.0]), (1, 3, 2))
    dphi_dy = np.zeros_like(dphi_dx)
    tau = supg_tau(ux, uy, k, zero, zero, k, dphi_dx, dphi_dy)
    assert tau[0, 0] == 0.0                                       # no flow
    pe = 1.0 * 0.5 / 2.0
    assert np.isclose(tau[0, 1], 0.5 / 2 * (1 / np.tanh(pe) - 1 / pe))
    assert np.isclose(tau[0, 2], 0.5 / 2)                         # coth(inf) - 0 = 1
    tau_t = supg_tau(ux, uy, k, zero, zero, k, dphi_dx, dphi_dy, dt=0.1)
    assert tau_t[0, 0] == 0.0 and 0 < tau_t[0, 2] < 0.05


# --------------------------------------------------------------- SUPG steady
def test_supg_is_nodally_exact_on_quads():
    sol = solve(boundary_layer("quad", supg=True))
    assert np.allclose(sol.T, exact_layer(sol.mesh.x, sol.mesh.y), atol=1e-10)


def test_galerkin_oscillates_at_high_cell_peclet():
    sol = solve(boundary_layer("quad", supg=False))
    assert sol.T.min() < -0.1


@pytest.mark.parametrize("name", ["triangle", "triangle6", "quad9"])
def test_supg_damps_oscillations_on_other_elements(name):
    galerkin = solve(boundary_layer(name, supg=False))
    supg = solve(boundary_layer(name, supg=True))
    assert supg.T.min() > -0.2
    assert supg.T.min() > galerkin.T.min()
    assert supg.T.max() < 1.0 + 1e-6
    err_g = np.abs(galerkin.T - exact_layer(galerkin.mesh.x, galerkin.mesh.y)).max()
    err_s = np.abs(supg.T - exact_layer(supg.mesh.x, supg.mesh.y)).max()
    assert err_s < err_g


@pytest.mark.parametrize("name", ["triangle", "quad9"])
def test_supg_converges_with_refinement(name):
    errs = []
    for elems in (10, 20, 40):
        sol = solve(boundary_layer(name, elems=elems, supg=True))
        errs.append(error_norms(sol.mesh, sol.T, exact_layer)["L2"])
    assert errs[2] < errs[1] < errs[0]
    assert errs[2] < errs[0] / 4


@pytest.mark.parametrize("name", ["quad", "triangle"])
def test_skew_advection_stays_bounded(name):
    inlet = lambda x, y: np.where(y > 0.3, 1.0, 0.0)
    kw = dict(a=1.0, b=1.0, elems=20, elem_type=name,
              bc_type={"left": DIRICHLET, "bottom": DIRICHLET},
              bc_val={"left": inlet, "bottom": 0.0},
              material=lambda x, y: (1e-4, 0.0, 0.0, 1e-4, 0.0),
              velocity=lambda x, y: (1.0, 1.0))
    galerkin = solve(Problem(supg=False, **kw))
    supg = solve(Problem(supg=True, **kw))
    assert -0.2 < supg.T.min() and supg.T.max() < 1.1
    assert supg.T.min() > galerkin.T.min()
    assert supg.T.max() < galerkin.T.max()
    # the step is transported along the diagonal: T ~ 1 above y = x + 0.3, ~ 0 below
    m = supg.mesh
    above = (m.y - m.x > 0.45) & (m.x > 0.1)
    below = (m.y - m.x < 0.15) & (m.y > 0.1)
    assert supg.T[above].min() > 0.8 and supg.T[below].max() < 0.2


# ------------------------------------------------------------------- solvers
def test_gmres_matches_direct():
    p = boundary_layer("triangle6", elems=8)
    direct = solve(p)
    gm = solve(p, method="gmres", tol=1e-12)
    assert gm.info["converged"]
    assert np.allclose(gm.T, direct.T, atol=1e-8)


def test_cg_refuses_nonsymmetric_system():
    with pytest.raises(ValueError, match="gmres"):
        solve(boundary_layer("quad"), method="cg")


def test_all_dirichlet_sides_accept_time_independent_problem():
    p = boundary_layer("quad")
    p.bc_type = {s: DIRICHLET for s in SIDES}
    p.bc_val = {s: exact_layer for s in SIDES}
    sol = solve(p)
    assert np.allclose(sol.T, exact_layer(sol.mesh.x, sol.mesh.y), atol=1e-10)
