"""Gradients, fluxes and error norms."""

import numpy as np
import pytest

from aa540fem import DIRICHLET, ELEMENTS, Problem, element_gradient, error_norms, geometry, solve


@pytest.mark.parametrize("name", sorted(ELEMENTS))
def test_gradient_and_flux_of_linear_field(name):
    mesh = geometry(2.0, 3.0, 3, name)
    T = 2 * mesh.x - mesh.y
    kappa = (2.0, 0.5, 0.5, 1.0)
    cf = element_gradient(mesh, T, material=lambda x, y: (*kappa, 0.0))
    assert cf.grad.shape == (mesh.n_elems, 2)
    assert np.allclose(cf.grad, [2.0, -1.0])
    # q = -kappa grad T = -[(2*2 + 0.5*-1), (0.5*2 + 1*-1)] = -[3.5, 0]
    assert np.allclose(cf.flux, [-3.5, 0.0])
    assert np.allclose(cf.centroid, mesh.centroids())


@pytest.mark.parametrize("name, field", [
    ("triangle", lambda x, y: 1 + x - 2 * y),
    ("quad", lambda x, y: 1 + x - 2 * y + 0.5 * x * y),
    ("triangle6", lambda x, y: x ** 2 - y ** 2 + x * y),
    ("quad9", lambda x, y: x ** 2 * y ** 2 - 3 * x),
])
def test_error_norms_vanish_for_representable_fields(name, field):
    mesh = geometry(2.0, 3.0, 3, name)
    T = field(mesh.x, mesh.y)
    h = 1e-6
    grad = lambda x, y: ((field(x + h, y) - field(x - h, y)) / (2 * h),
                         (field(x, y + h) - field(x, y - h)) / (2 * h))
    e = error_norms(mesh, T, field, grad)
    assert e["L2"] < 1e-12
    assert e["H1"] < 1e-8


def test_error_norms_measure_interpolation_error():
    mesh = geometry(1.0, 1.0, 4, "triangle")
    exact = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)
    e = error_norms(mesh, exact(mesh.x, mesh.y), exact)
    assert 1e-3 < e["L2"] < 1e-1
    assert "H1" not in e


def test_solution_flux_helper_uses_problem_material():
    mesh = geometry(1.0, 1.0, 4, "quad")
    p = Problem(mesh=mesh, bc_type={"left": DIRICHLET, "right": DIRICHLET},
                bc_val={"left": 0.0, "right": 1.0},
                material=lambda x, y: (3.0, 0.0, 0.0, 1.0, 0.0))
    sol = solve(p)
    cf = sol.flux()
    assert np.allclose(cf.grad, [1.0, 0.0])
    assert np.allclose(cf.flux, [-3.0, 0.0])
