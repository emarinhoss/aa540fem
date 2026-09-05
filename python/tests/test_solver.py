"""Verification of the Python port against analytical solutions."""

import numpy as np
import pytest

from aa540fem import (
    Problem,
    gauss_legendre_quad,
    gauss_trgl,
    geometry,
    interpfunc,
    quadrature_rule,
    solve,
)
from aa540fem.solver import DIRICHLET, NEUMANN

ELEM_TYPES = [1, 2, 3]


# ----------------------------------------------------------------- quadrature
def test_gauss_legendre_matches_matlab_tables():
    xi, wi = gauss_legendre_quad(1)
    assert np.allclose(xi, [0]) and np.allclose(wi, [2])
    xi, wi = gauss_legendre_quad(2)
    assert np.allclose(np.sort(xi), [-1 / np.sqrt(3), 1 / np.sqrt(3)])
    assert np.allclose(wi, [1, 1])
    xi, wi = gauss_legendre_quad(3)
    assert np.allclose(np.sort(xi), [-np.sqrt(15) / 5, 0, np.sqrt(15) / 5])
    assert np.allclose(wi, np.array([5, 8, 5]) / 9)
    xi, wi = gauss_legendre_quad(4)
    r1 = np.sqrt((3 - 2 * np.sqrt(6 / 5)) / 7)
    r2 = np.sqrt((3 + 2 * np.sqrt(6 / 5)) / 7)
    assert np.allclose(np.sort(xi), [-r2, -r1, r1, r2])
    expected = [(18 - np.sqrt(30)) / 36] * 2 + [(18 + np.sqrt(30)) / 36] * 2
    assert np.allclose(np.sort(wi), np.sort(expected))


@pytest.mark.parametrize("m", [1, 3, 4, 6, 7, 9, 12, 13])
def test_gauss_trgl_integrates_polynomials(m):
    xi, eta, w = gauss_trgl(m)
    assert np.isclose(w.sum(), 1.0)
    # exact: int over unit triangle of x*y dA = 1/24, times 2 for normalised weights
    if m >= 3:
        assert np.isclose(np.sum(w * xi * eta), 2 / 24)


@pytest.mark.parametrize("elem_type, total", [(1, 0.5), (2, 4.0), (3, 4.0)])
def test_quadrature_rule_weights(elem_type, total):
    _, _, w = quadrature_rule(elem_type, 3)
    assert np.isclose(w.sum(), total)


# ------------------------------------------------------------ shape functions
@pytest.mark.parametrize("elem_type", ELEM_TYPES)
def test_shape_functions_partition_of_unity(elem_type):
    xi, eta, _ = quadrature_rule(elem_type, 3)
    phi, dxi, deta = interpfunc(elem_type, xi, eta)
    assert np.allclose(phi.sum(axis=1), 1)
    assert np.allclose(dxi.sum(axis=1), 0)
    assert np.allclose(deta.sum(axis=1), 0)


@pytest.mark.parametrize("elem_type", ELEM_TYPES)
def test_shape_functions_derivatives_by_finite_difference(elem_type):
    rng = np.random.default_rng(0)
    if elem_type == 1:
        pts = rng.random((5, 2)) * 0.4
    else:
        pts = rng.uniform(-0.9, 0.9, (5, 2))
    h = 1e-6
    for xi, eta in pts:
        _, dxi, deta = interpfunc(elem_type, [xi], [eta])
        phi = lambda a, b: interpfunc(elem_type, [a], [b])[0]
        fd_xi = (phi(xi + h, eta) - phi(xi - h, eta)) / (2 * h)
        fd_eta = (phi(xi, eta + h) - phi(xi, eta - h)) / (2 * h)
        assert np.allclose(dxi, fd_xi, atol=1e-6)
        assert np.allclose(deta, fd_eta, atol=1e-6)


def test_shape_functions_are_nodal():
    # Each shape function is 1 at its own node and 0 at the others.
    # Gmsh / VTK local node ordering
    nodes = {
        1: [(0, 0), (1, 0), (0, 1)],
        "triangle6": [(0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)],
        2: [(-1, -1), (1, -1), (1, 1), (-1, 1)],
        3: [(-1, -1), (1, -1), (1, 1), (-1, 1), (0, -1), (1, 0), (0, 1), (-1, 0), (0, 0)],
    }
    for t, pts in nodes.items():
        xi, eta = np.array(pts).T
        phi, _, _ = interpfunc(t, xi, eta)
        assert np.allclose(phi, np.eye(len(pts)))


# ------------------------------------------------------------------- geometry
@pytest.mark.parametrize("elem_type, n_elems, n_nodes", [(1, 2 * 16, 25), (2, 16, 25), (3, 16, 81)])
def test_geometry_sizes(elem_type, n_elems, n_nodes):
    mesh = geometry(2.0, 3.0, 4, elem_type)
    assert mesh.n_elems == n_elems
    assert mesh.n_nodes == n_nodes
    assert mesh.conn.shape[1] == mesh.nodes_per_element
    assert np.allclose(mesh.x[mesh.bc_nodes["left"]], 0)
    assert np.allclose(mesh.x[mesh.bc_nodes["right"]], 2)
    assert np.allclose(mesh.y[mesh.bc_nodes["bottom"]], 0)
    assert np.allclose(mesh.y[mesh.bc_nodes["top"]], 3)


@pytest.mark.parametrize("elem_type", ELEM_TYPES)
def test_geometry_positive_jacobian_and_full_area(elem_type):
    from aa540fem.element import elem_eqn
    mesh = geometry(2.0, 3.0, 5, elem_type)
    xi, eta, w = quadrature_rule(elem_type, 3)
    phi, dxi, deta = interpfunc(elem_type, xi, eta)
    xe, ye = mesh.x[mesh.conn], mesh.y[mesh.conn]
    hs = (xe @ dxi.T) * (ye @ deta.T) - (xe @ deta.T) * (ye @ dxi.T)
    assert np.all(hs > 0)
    # With f = 1 the sum of the load vectors is the domain area.
    _, fe = elem_eqn(xe, ye, phi, dxi, deta, w, material=lambda x, y: (1, 0, 0, 1, 1.0))
    assert np.isclose(fe.sum(), 6.0)


# ---------------------------------------------------------------- full solves
def _problem(elem_type, elems, **kw):
    p = Problem(a=2.0, b=3.0, elems=elems, elem_type=elem_type, order=3)
    for k, v in kw.items():
        setattr(p, k, v)
    return p


@pytest.mark.parametrize("elem_type", ELEM_TYPES)
def test_linear_field_with_neumann_is_exact(elem_type):
    # T = x with kappa = I: q_n = 1 on the right, 0 on top/bottom, T = 0 on the left.
    p = _problem(elem_type, 6,
                 bc_type={"top": NEUMANN, "right": NEUMANN, "left": DIRICHLET, "bottom": NEUMANN},
                 bc_val={"top": 0.0, "right": 1.0, "left": 0.0, "bottom": 0.0})
    sol = solve(p)
    assert np.allclose(sol.T, sol.mesh.x, atol=1e-10)


@pytest.mark.parametrize("elem_type", ELEM_TYPES)
def test_anisotropic_quadratic_field_with_dirichlet(elem_type):
    # T = x^2 + y^2, kappa = [[2, .5], [.5, 1]]  ->  div(kappa grad T) = 6, so f = -6.
    exact = lambda x, y: x ** 2 + y ** 2
    p = _problem(elem_type, 8,
                 bc_type={s: DIRICHLET for s in ("top", "right", "left", "bottom")},
                 bc_val={s: exact for s in ("top", "right", "left", "bottom")},
                 material=lambda x, y: (2.0, 0.5, 0.5, 1.0, -6.0))
    sol = solve(p)
    err = np.abs(sol.T - exact(sol.mesh.x, sol.mesh.y)).max()
    if elem_type == 3:
        assert err < 1e-9          # quadratic elements reproduce T exactly
    else:
        assert err < 0.05


@pytest.mark.parametrize("elem_type", ELEM_TYPES)
def test_anisotropic_quadratic_field_with_neumann(elem_type):
    # Same field as above but with fluxes n.(kappa grad T) on the right and top:
    #   kappa grad T = (4x + y, x + 2y); right (n = +x): 4a + y; top (n = +y): x + 2b.
    a, b = 2.0, 3.0
    exact = lambda x, y: x ** 2 + y ** 2
    p = _problem(elem_type, 8,
                 bc_type={"top": NEUMANN, "right": NEUMANN, "left": DIRICHLET, "bottom": DIRICHLET},
                 bc_val={"top": lambda x, y: x + 2 * b, "right": lambda x, y: 4 * a + y,
                         "left": exact, "bottom": exact},
                 material=lambda x, y: (2.0, 0.5, 0.5, 1.0, -6.0))
    sol = solve(p)
    err = np.abs(sol.T - exact(sol.mesh.x, sol.mesh.y)).max()
    if elem_type == 3:
        assert err < 1e-9
    else:
        assert err < 0.1


@pytest.mark.parametrize("elem_type, rate", [(1, 1.8), (2, 1.8), (3, 3.5)])
def test_convergence_manufactured_solution(elem_type, rate):
    # T = sin(pi x / a) sin(pi y / b), kappa = I, homogeneous Dirichlet.
    a, b = 2.0, 3.0
    exact = lambda x, y: np.sin(np.pi * x / a) * np.sin(np.pi * y / b)
    c = np.pi ** 2 * (1 / a ** 2 + 1 / b ** 2)
    errs = []
    for elems in (4, 8):
        p = _problem(elem_type, elems,
                     bc_type={s: DIRICHLET for s in ("top", "right", "left", "bottom")},
                     bc_val={s: 0.0 for s in ("top", "right", "left", "bottom")},
                     material=lambda x, y: (1.0, 0.0, 0.0, 1.0, c * exact(x, y)))
        sol = solve(p)
        errs.append(np.abs(sol.T - exact(sol.mesh.x, sol.mesh.y)).max())
    observed = np.log2(errs[0] / errs[1])
    assert observed > rate, f"observed rate {observed:.2f} with errors {errs}"


@pytest.mark.parametrize("elem_type", ELEM_TYPES)
def test_default_problem_runs_like_main_m(elem_type):
    # The defaults reproduce main.m: T = 100 on top, 0 on the bottom, insulated sides.
    sol = solve(Problem(elems=10, elem_type=elem_type))
    exact = 100.0 * sol.mesh.y / 6.0
    assert np.allclose(sol.T, exact, atol=1e-8)


def test_under_integration_warns():
    with pytest.warns(UserWarning, match="under-integrates"):
        from aa540fem.solver import assemble
        assemble(geometry(1.0, 1.0, 2, 3), order=1)


@pytest.mark.parametrize("elem_type", ELEM_TYPES + ["triangle6"])
def test_cg_matches_direct(elem_type):
    a, b = 2.0, 3.0
    exact = lambda x, y: np.sin(np.pi * x / a) * np.sin(np.pi * y / b)
    c = np.pi ** 2 * (1 / a ** 2 + 1 / b ** 2)
    p = _problem(elem_type, 12,
                 bc_type={s: DIRICHLET for s in ("top", "right", "left", "bottom")},
                 bc_val={s: 0.0 for s in ("top", "right", "left", "bottom")},
                 material=lambda x, y: (1.0, 0.0, 0.0, 1.0, c * exact(x, y)))
    direct = solve(p)
    cg = solve(p, method="cg", tol=1e-12)
    assert cg.info["converged"]
    assert np.allclose(cg.T, direct.T, atol=1e-8)


def test_unlisted_tags_are_natural():
    # Only top/bottom given: left/right default to zero flux, same as listing them Neumann.
    p = Problem(elems=6, elem_type=2, bc_type={"top": DIRICHLET, "bottom": DIRICHLET},
                bc_val={"top": 100.0, "bottom": 0.0})
    sol = solve(p)
    assert np.allclose(sol.T, 100.0 * sol.mesh.y / 6.0, atol=1e-8)


def test_unknown_method_is_rejected():
    with pytest.raises(ValueError):
        solve(Problem(elems=2, elem_type=2), method="bicg")


def test_pure_neumann_is_rejected():
    with pytest.raises(ValueError):
        all_neumann = {s: NEUMANN for s in ("top", "right", "left", "bottom")}
        solve(Problem(elems=4, elem_type=2, bc_type=all_neumann))
