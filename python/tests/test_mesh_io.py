"""Unstructured meshes, mesh files and ParaView output."""

import pathlib

import numpy as np
import pytest
from scipy.spatial import Delaunay

from aa540fem import DIRICHLET, NEUMANN, Mesh, Problem, error_norms, geometry, solve

meshio = pytest.importorskip("meshio")
from aa540fem.mesh_io import from_meshio, read_mesh, to_meshio, write_vtk  # noqa: E402

EXAMPLES = pathlib.Path(__file__).resolve().parent.parent / "examples"
ANNULUS = {n: EXAMPLES / f"annulus_{n}.msh" for n in ("tri", "tri6", "quad", "quad9")}


def delaunay_rectangle(a=2.0, b=1.0, n=9, seed=0):
    """Unstructured triangle mesh of [0, a] x [0, b] with tagged sides."""
    rng = np.random.default_rng(seed)
    xs = np.linspace(0, a, n)
    ys = np.linspace(0, b, n // 2 + 1)
    edge = np.vstack([
        np.column_stack([xs, np.zeros_like(xs)]),
        np.column_stack([xs, np.full_like(xs, b)]),
        np.column_stack([np.zeros_like(ys[1:-1]), ys[1:-1]]),
        np.column_stack([np.full_like(ys[1:-1], a), ys[1:-1]]),
    ])
    interior = rng.random((60, 2)) * [a, b] * 0.9 + [0.05 * a, 0.05 * b]
    points = np.vstack([edge, interior])
    mesh = Mesh(points, {"triangle": Delaunay(points).simplices}, {})
    mesh.check_orientation()
    faces = mesh.boundary_faces()
    x, y = mesh.x, mesh.y
    on = lambda arr, v: np.all(np.isclose(arr[faces[:, :2]], v), axis=1)
    mesh.boundary = {
        "left": faces[on(x, 0)], "right": faces[on(x, a)],
        "bottom": faces[on(y, 0)], "top": faces[on(y, b)],
    }
    assert sum(len(e) for e in mesh.boundary.values()) == len(faces)
    return mesh


def test_unstructured_triangles_reproduce_linear_field():
    mesh = delaunay_rectangle()
    exact = lambda x, y: x + 2 * y
    # kappa = I: n.grad T = 1 on the right (n = +x), 2 on the top (n = +y)
    p = Problem(mesh=mesh,
                bc_type={"left": DIRICHLET, "bottom": DIRICHLET, "right": NEUMANN, "top": NEUMANN},
                bc_val={"left": exact, "bottom": exact, "right": 1.0, "top": 2.0})
    sol = solve(p)
    assert np.allclose(sol.T, exact(mesh.x, mesh.y), atol=1e-10)
    grad = lambda x, y: (np.ones_like(x), 2 * np.ones_like(x))
    assert error_norms(mesh, sol.T, exact, grad)["H1"] < 1e-10


def test_check_orientation_flips_clockwise_elements():
    mesh = geometry(1.0, 1.0, 3, "quad9")
    ref = mesh.cells["quad9"].copy()
    mesh.cells["quad9"][::2] = ref[::2][:, [0, 3, 2, 1, 7, 6, 5, 4, 8]]
    assert (mesh.jacobian_at_centroids()["quad9"] < 0).sum() == len(ref[::2])
    assert mesh.check_orientation() == len(ref[::2])
    assert np.array_equal(mesh.cells["quad9"], ref)
    assert mesh.check_orientation() == 0


def test_boundary_faces_of_structured_mesh():
    mesh = geometry(2.0, 1.0, 4, "triangle6")
    faces = mesh.boundary_faces()
    tagged = np.vstack(list(mesh.boundary.values()))
    assert faces.shape == tagged.shape
    assert set(map(tuple, np.sort(faces, axis=1))) == set(map(tuple, np.sort(tagged, axis=1)))


def test_unknown_tag_lists_available_tags():
    with pytest.raises(ValueError, match="inner"):
        solve(Problem(mesh=read_mesh(ANNULUS["tri"]), bc_type={"top": 0}, bc_val={"top": 0.0}))


@pytest.mark.parametrize("name", sorted(ANNULUS))
def test_read_gmsh_annulus(name):
    mesh = read_mesh(ANNULUS[name])
    assert set(mesh.tags) == {"inner", "outer"}
    assert all(d.min() > 0 for d in mesh.jacobian_at_centroids().values())
    r = np.hypot(mesh.x, mesh.y)
    assert np.allclose(r[mesh.bc_nodes["inner"]], 1.0)
    assert np.allclose(r[mesh.bc_nodes["outer"]], 2.0)
    # every node is used and the corner nodes of the boundary edges are on the circles
    assert np.unique(np.concatenate([c.ravel() for c in mesh.cells.values()])).size == mesh.n_nodes


def _annulus_error(name):
    mesh = read_mesh(ANNULUS[name])
    p = Problem(mesh=mesh, bc_type={"inner": DIRICHLET, "outer": DIRICHLET},
                bc_val={"inner": 0.0, "outer": 1.0})
    sol = solve(p)
    exact = lambda x, y: np.log(np.hypot(x, y)) / np.log(2)
    return error_norms(mesh, sol.T, exact)["L2"], sol


def test_annulus_curved_quadratic_elements_are_more_accurate():
    e_tri, _ = _annulus_error("tri")
    e_tri6, _ = _annulus_error("tri6")
    e_quad, _ = _annulus_error("quad")
    e_quad9, _ = _annulus_error("quad9")
    assert e_tri < 1e-2 and e_quad < 1e-2
    assert e_tri6 < e_tri / 10
    assert e_quad9 < e_quad / 10


def test_cg_matches_direct_on_unstructured_mesh():
    mesh = read_mesh(ANNULUS["tri6"])
    p = Problem(mesh=mesh, bc_type={"inner": DIRICHLET, "outer": DIRICHLET},
                bc_val={"inner": 0.0, "outer": 1.0})
    direct = solve(p)
    cg = solve(p, method="cg", tol=1e-12)
    assert cg.info["converged"]
    assert np.allclose(cg.T, direct.T, atol=1e-8)


def test_vtk_round_trip(tmp_path):
    _, sol = _annulus_error("quad9")
    out = tmp_path / "annulus.vtu"
    sol.save(out)
    back = meshio.read(out)
    assert np.allclose(back.point_data["T"], sol.T)
    assert [c.type for c in back.cells] == ["quad9"]
    assert np.array_equal(back.cells[0].data, sol.mesh.cells["quad9"])
    flux = np.vstack(back.cell_data["flux"])
    assert flux.shape == (sol.mesh.n_elems, 2)
    # re-import through the reader path as well
    again = from_meshio(back)
    assert again.n_nodes == sol.mesh.n_nodes and again.tags == ["boundary"]


def test_write_vtk_validates_cell_data(tmp_path):
    mesh = geometry(1.0, 1.0, 2, 2)
    with pytest.raises(ValueError):
        write_vtk(tmp_path / "bad.vtu", mesh, cell_data={"q": np.zeros((3, 2))})
    m = to_meshio(mesh, point_data={"T": mesh.x})
    assert m.points.shape == (mesh.n_nodes, 3)
