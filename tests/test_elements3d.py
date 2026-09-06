"""3-D reference elements, quadrature, mappings and meshes."""

from math import factorial

import numpy as np
import pytest

from aa540fem.core.elements import ELEMENTS, FACE_CORNERS, PRESSURE_ELEMENT, get_element
from aa540fem.core.mesh import BOX_SIDES, Mesh, box
from aa540fem.core.quadrature import gauss_tetra, quadrature_points
from aa540fem.transport.element import (
    element_metric_nd,
    jacobian_nd,
    map_gradients_nd,
    physical_laplacian_nd,
)

ELEMENTS_3D = sorted(n for n, e in ELEMENTS.items() if e.dim == 3)


def interior_points(el, k=6):
    rng = np.random.default_rng(1)
    if el.family == "tetra":
        return tuple(rng.uniform(0.05, 0.3, (k, 3)).T)
    return tuple(rng.uniform(-0.8, 0.8, (k, 3)).T)


@pytest.mark.parametrize("name", ELEMENTS_3D)
def test_nodal_property_and_partition_of_unity(name):
    el = get_element(name)
    nodes = np.array(el.nodes, dtype=float)
    phi = el.shape(*nodes.T)[0]
    assert np.allclose(phi, np.eye(el.n_nodes))
    phi = el.shape(*interior_points(el))[0]
    assert np.allclose(phi.sum(axis=1), 1.0)
    assert el.dim == 3 and el.n_corners in (4, 8) and len(el.edges) in (6, 12)
    assert sorted(el.reverse) == list(range(el.n_nodes))


@pytest.mark.parametrize("name", ELEMENTS_3D)
def test_derivatives_and_hessians_match_finite_differences(name):
    el = get_element(name)
    pts = np.array(interior_points(el))
    h = 1e-6
    phi, grads = el.shape_at(tuple(pts))
    hess = el.hessian(*pts)
    pairs = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    for k in range(3):
        e = np.zeros((3, 1))
        e[k] = h
        fd = (el.shape(*(pts + e))[0] - el.shape(*(pts - e))[0]) / (2 * h)
        assert np.allclose(fd, grads[k], atol=1e-7)
    for m, (p, q) in enumerate(pairs):
        e = np.zeros((3, 1))
        e[p] = h
        fd = (el.shape(*(pts + e))[q + 1] - el.shape(*(pts - e))[q + 1]) / (2 * h)
        assert np.allclose(fd, hess[m], atol=1e-6)


@pytest.mark.parametrize("degree", [1, 2, 3, 4, 5])
def test_tetrahedral_rule_integrates_monomials(degree):
    xi, eta, zeta, w = gauss_tetra(degree)
    assert (w > 0).all() and abs(w.sum() - 1 / 6) < 1e-15
    for a in range(degree + 1):
        for b in range(degree + 1 - a):
            for c in range(degree + 1 - a - b):
                exact = factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 3)
                assert abs(np.sum(w * xi ** a * eta ** b * zeta ** c) - exact) < 1e-15


def test_hexahedral_rule_and_element_volumes():
    coords, w = quadrature_points("hexahedron", 3)
    assert abs(w.sum() - 8.0) < 1e-14
    for name in ELEMENTS_3D:
        el = get_element(name)
        coords, w = quadrature_points(el.family, el.full_order)
        nodes = np.array(el.nodes, dtype=float)
        _, dnat = el.shape_at(coords)
        hs, _ = jacobian_nd(tuple(nodes[None, :, k] for k in range(3)), dnat)
        volume = np.sum(w * hs[0])
        assert abs(volume - (1 / 6 if el.family == "tetra" else 8.0)) < 1e-14


@pytest.mark.parametrize("name", ELEMENTS_3D)
def test_faces_have_outward_normals_and_cover_the_boundary(name):
    el = get_element(name)
    fel = get_element(el.face_type)
    nodes = np.array(el.nodes, dtype=float)
    centroid = np.array(el.centroid)
    corner_faces = set()
    for face in el.faces:
        assert len(face) == fel.n_nodes
        c = nodes[list(face[:fel.n_corners])]
        normal = np.cross(c[1] - c[0], c[2] - c[0])
        assert np.dot(normal, c[0] - centroid) > 0
        corner_faces.add(frozenset(face[:fel.n_corners]))
        # mid-edge nodes of the face are the midpoints of consecutive face corners
        k = fel.n_corners
        for m in range(k, min(fel.n_nodes, 2 * k)):
            a, b = face[m - k], face[(m - k + 1) % k]
            assert np.allclose(nodes[face[m]], 0.5 * (nodes[a] + nodes[b]))
    assert len(corner_faces) == el.n_faces
    # every edge of the element belongs to exactly two faces
    edge_count = {}
    for face in el.faces:
        k = fel.n_corners
        for i in range(k):
            key = frozenset((face[i], face[(i + 1) % k]))
            edge_count[key] = edge_count.get(key, 0) + 1
    assert set(edge_count) == {frozenset(e) for e in el.edges}
    assert set(edge_count.values()) == {2}


@pytest.mark.parametrize("name", ["tetra10", "hexahedron27"])
def test_laplacian_exact_on_affine_element_and_metric_measures_sizes(name):
    el = get_element(name)
    nodes = np.array(el.nodes, dtype=float)
    scale = np.array([2.0, 1.0, 0.5])                        # stretched, affine
    xe = tuple((scale[k] * nodes[:, k])[None, :] for k in range(3))
    pts = interior_points(el, 4)
    _, dnat = el.shape_at(pts)
    hs, inverse = jacobian_nd(xe, dnat)
    assert np.allclose(hs, np.prod(scale))
    lap = physical_laplacian_nd(inverse, el.hessian(*pts))
    # u = x^2 + 2 y^2 - z^2 interpolates exactly on quadratic elements: lap u = 4
    ue = xe[0] ** 2 + 2 * xe[1] ** 2 - xe[2] ** 2
    assert np.allclose(np.einsum("eqi,ei->eq", lap, ue), 4.0)
    grads = map_gradients_nd(inverse, dnat)
    y = scale[1] * pts[1]                                    # physical y of the points
    assert np.allclose(np.einsum("eqi,ei->eq", grads[1], ue), 4 * y[None, :])
    G = element_metric_nd(inverse, el.family)
    if el.family == "hexahedron":
        # rectangle box of sides 2 scale: G = diag(4 / h^2)
        h = 2 * scale
        assert np.allclose(G[0, 0], np.diag(4.0 / h ** 2))


def test_reverse_flips_the_orientation():
    for name in ELEMENTS_3D:
        el = get_element(name)
        nodes = np.array(el.nodes, dtype=float)
        _, dnat = el.shape_at(tuple([c] for c in el.centroid))
        hs, _ = jacobian_nd(tuple(nodes[None, :, k] for k in range(3)), dnat)
        hs_rev, _ = jacobian_nd(tuple(nodes[list(el.reverse)][None, :, k] for k in range(3)),
                                dnat)
        assert hs[0, 0] > 0 and np.isclose(hs_rev[0, 0], -hs[0, 0])


@pytest.mark.parametrize("name", ["hexahedron27", "tetra10"])
def test_box_mesh_is_consistent(name):
    mesh = box(2.0, 1.0, 0.5, (2, 1, 1), name)
    assert mesh.dim == 3 and mesh.n_nodes == 45
    assert mesh.n_elems == (2 if name == "hexahedron27" else 12)
    assert set(mesh.tags) == set(BOX_SIDES)
    assert (mesh.jacobian_at_centroids()[name] > 0).all() and mesh.check_orientation() == 0
    el = get_element(name)
    nc = FACE_CORNERS[el.face_type]
    volume = 0.0
    for b_faces in mesh.boundary.values():
        assert b_faces.shape[1] == get_element(el.face_type).n_nodes
    for tag, (axis, value) in {"left": (0, 0), "right": (0, 2), "bottom": (1, 0), "top": (1, 1),
                               "front": (2, 0), "back": (2, 0.5)}.items():
        assert np.allclose(mesh.points[mesh.boundary[tag][:, :nc], axis], value)
    # the boundary faces of the whole mesh are exactly the union of the tags
    faces = mesh.boundary_faces()
    assert faces.shape[0] == sum(f.shape[0] for f in mesh.boundary.values())
    # element volumes sum to the box volume (quadrature of the constant 1)
    el = get_element(name)
    coords, w = quadrature_points(el.family, el.full_order)
    _, dnat = el.shape_at(coords)
    conn = mesh.cells[name]
    hs, _ = jacobian_nd(tuple(mesh.points[conn, k] for k in range(3)), dnat)
    volume = np.sum(w[None, :] * hs)
    assert abs(volume - 1.0) < 1e-12
    assert np.allclose(mesh.nodal_size("min"), 0.5)
    assert PRESSURE_ELEMENT[name].n_nodes == el.n_corners


def test_mesh_rejects_mixed_dimensions():
    mesh = box(1.0, 1.0, 1.0, 1)
    with pytest.raises(ValueError, match="mix"):
        Mesh(mesh.points, {"hexahedron27": mesh.cells["hexahedron27"],
                           "quad9": mesh.boundary["left"]}, {})
    with pytest.raises(ValueError):
        Mesh(mesh.points[:, :2], mesh.cells, {})


def test_read_gmsh_tetrahedral_box():
    pytest.importorskip("meshio")
    import pathlib

    from aa540fem import read_mesh

    path = pathlib.Path(__file__).resolve().parent.parent / "examples" / "meshes" / "box_tet10.msh"
    mesh = read_mesh(path)
    assert mesh.dim == 3 and set(mesh.cells) == {"tetra10"}
    assert set(mesh.tags) == set(BOX_SIDES)
    assert (mesh.jacobian_at_centroids()["tetra10"] > 0).all()
    faces = mesh.boundary_faces()
    assert faces.shape[0] == sum(f.shape[0] for f in mesh.boundary.values())
    # tagged faces are outward: corner normals point away from the box centre
    centre = mesh.points.mean(axis=0)
    for tag in mesh.tags:
        c = mesh.points[mesh.boundary[tag][:, :3]]
        n = np.cross(c[:, 1] - c[:, 0], c[:, 2] - c[:, 0])
        assert np.all(np.einsum("ij,ij->i", n, c[:, 0] - centre) > 0)
