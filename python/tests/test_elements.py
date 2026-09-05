"""Reference-element registry checks."""

import numpy as np
import pytest

from aa540fem import ELEMENTS, get_element
from aa540fem.elements import QUAD, QUAD9, TRIANGLE, TRIANGLE6

# Natural coordinates of the nodes in Gmsh ordering
NODE_COORDS = {
    "triangle": [(0, 0), (1, 0), (0, 1)],
    "triangle6": [(0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)],
    "quad": [(-1, -1), (1, -1), (1, 1), (-1, 1)],
    "quad9": [(-1, -1), (1, -1), (1, 1), (-1, 1), (0, -1), (1, 0), (0, 1), (-1, 0), (0, 0)],
}


def test_registry_and_aliases():
    assert set(ELEMENTS) == set(NODE_COORDS)
    assert get_element(1) is TRIANGLE
    assert get_element(2) is QUAD
    assert get_element(3) is QUAD9
    assert get_element("triangle6") is TRIANGLE6
    assert get_element(TRIANGLE6) is TRIANGLE6
    with pytest.raises(ValueError):
        get_element(4)
    with pytest.raises(ValueError):
        get_element("hexahedron")


@pytest.mark.parametrize("name", sorted(ELEMENTS))
def test_nodal_and_partition_of_unity(name):
    el = get_element(name)
    xi, eta = np.array(NODE_COORDS[name], dtype=float).T
    phi, dxi, deta = el.shape(xi, eta)
    assert el.n_nodes == len(xi)
    assert np.allclose(phi, np.eye(el.n_nodes))
    q_xi, q_eta, w = el.quadrature()
    phi, dxi, deta = el.shape(q_xi, q_eta)
    assert np.allclose(phi.sum(axis=1), 1)
    assert np.allclose(dxi.sum(axis=1), 0)
    assert np.allclose(deta.sum(axis=1), 0)
    assert np.isclose(w.sum(), 0.5 if el.family == "triangle" else 4.0)


@pytest.mark.parametrize("name", sorted(ELEMENTS))
def test_derivatives_match_finite_differences(name):
    el = get_element(name)
    rng = np.random.default_rng(1)
    pts = rng.random((6, 2)) * 0.4 if el.family == "triangle" else rng.uniform(-0.9, 0.9, (6, 2))
    h = 1e-6
    for xi, eta in pts:
        _, dxi, deta = el.shape([xi], [eta])
        fd_xi = (el.shape([xi + h], [eta])[0] - el.shape([xi - h], [eta])[0]) / (2 * h)
        fd_eta = (el.shape([xi], [eta + h])[0] - el.shape([xi], [eta - h])[0]) / (2 * h)
        assert np.allclose(dxi, fd_xi, atol=1e-6)
        assert np.allclose(deta, fd_eta, atol=1e-6)


@pytest.mark.parametrize("name", sorted(ELEMENTS))
def test_faces_walk_the_boundary_counter_clockwise(name):
    el = get_element(name)
    coords = np.array(NODE_COORDS[name], dtype=float)
    for face in el.faces:
        start, end = coords[face[0]], coords[face[1]]
        if len(face) == 3:
            assert np.allclose(coords[face[2]], 0.5 * (start + end))
        # the centroid must lie to the left of every edge
        edge = end - start
        to_centroid = np.array(el.centroid) - start
        assert edge[0] * to_centroid[1] - edge[1] * to_centroid[0] > 0
    # consecutive faces share their end/start node
    for f, g in zip(el.faces, el.faces[1:] + el.faces[:1]):
        assert f[1] == g[0]


@pytest.mark.parametrize("name", sorted(ELEMENTS))
def test_reverse_flips_orientation(name):
    el = get_element(name)
    coords = np.array(NODE_COORDS[name], dtype=float)
    x, y = coords[:, 0], coords[:, 1]
    _, dxi, deta = el.shape([el.centroid[0]], [el.centroid[1]])

    def det(x, y):
        return (x @ dxi[0]) * (y @ deta[0]) - (x @ deta[0]) * (y @ dxi[0])

    r = list(el.reverse)
    assert sorted(r) == list(range(el.n_nodes))
    assert det(x, y) > 0
    assert det(x[r], y[r]) < 0
    # the reversed element is still nodal: shape function r[k] is 1 at node r[k]
    phi, _, _ = el.shape(x[r], y[r])
    assert np.allclose(phi[:, r], np.eye(el.n_nodes))
