"""Turbulence in 3-D: wall distance to faces, the extrusion utility, the SA solver
and the RANS coupling against their 2-D results on extruded meshes."""

import numpy as np
import pytest

from aa540fem import geometry
from aa540fem.core.mesh import box
from aa540fem.incompressible import FlowProblem
from aa540fem.turbulence import SpalartAllmaras, SpalartAllmarasSolver

NU = 1e-3


def test_wall_distance_to_faces():
    mesh = box(2.0, 1.0, 0.5, (4, 2, 1), "hexahedron27")
    assert np.allclose(mesh.wall_distance("bottom"), mesh.y)
    assert np.allclose(mesh.wall_distance(["bottom", "left"]), np.minimum(mesh.x, mesh.y))
    tet = box(1.0, 1.0, 1.0, 2, "tetra10")
    assert np.allclose(tet.wall_distance("front"), tet.z)
    assert np.allclose(tet.wall_distance(["top", "right"]), np.minimum(1 - tet.y, 1 - tet.x))


def test_point_triangle_distance_regions():
    from aa540fem.core.wall_distance import point_triangle_distance

    A = np.array([[0.0, 0.0, 0.0]])
    B = np.array([[1.0, 0.0, 0.0]])
    C = np.array([[0.0, 1.0, 0.0]])
    P = np.array([[0.2, 0.2, 0.5],      # above the face
                  [-1.0, -1.0, 0.0],    # vertex A region
                  [2.0, 0.0, 0.0],      # vertex B region
                  [0.5, -1.0, 0.0],     # edge AB region
                  [1.0, 1.0, 0.0],      # edge BC region
                  [-1.0, 0.5, 1.0]])    # edge AC region, lifted
    d = point_triangle_distance(P, A, B, C)[:, 0]
    assert np.allclose(d, [0.5, np.sqrt(2), 1.0, 1.0, np.sqrt(0.5), np.sqrt(2)])


def test_extrusion_of_a_quad9_mesh():
    mesh2 = geometry(2.0, 1.0, 3, "quad9")
    mesh3 = mesh2.extrude(0.5, layers=2)
    assert mesh3.dim == 3 and set(mesh3.cells) == {"hexahedron27"}
    assert mesh3.n_elems == 2 * mesh2.n_elems and mesh3.n_nodes == 5 * mesh2.n_nodes
    assert (mesh3.jacobian_at_centroids()["hexahedron27"] > 0).all()
    assert set(mesh3.tags) == {"top", "right", "left", "bottom", "front", "back"}
    for tag in ("top", "right", "left", "bottom"):
        assert mesh3.boundary[tag].shape == (2 * mesh2.boundary[tag].shape[0], 9)
    assert np.allclose(mesh3.z[mesh3.bc_nodes["front"]], 0.0)
    assert np.allclose(mesh3.z[mesh3.bc_nodes["back"]], 0.5)
    # the tagged faces are exactly the boundary faces of the extruded mesh
    faces = mesh3.boundary_faces()
    assert faces.shape[0] == sum(f.shape[0] for f in mesh3.boundary.values())
    # volume from the quadrature of the constant 1
    from aa540fem.core.elements import get_element
    from aa540fem.core.quadrature import quadrature_points
    from aa540fem.transport.element import jacobian_nd

    el = get_element("hexahedron27")
    coords, w = quadrature_points(el.family, el.full_order)
    _, dnat = el.shape_at(coords)
    conn = mesh3.cells["hexahedron27"]
    hs, _ = jacobian_nd(tuple(mesh3.points[conn, k] for k in range(3)), dnat)
    assert abs(np.sum(w[None, :] * hs) - 1.0) < 1e-12


def extruded_pair(elems=4):
    """A 2-D quad9 channel and its one-layer extrusion, with the same shear flow."""
    mesh2 = geometry(2.0, 1.0, elems, "quad9")
    mesh3 = mesh2.extrude(0.25, layers=1)
    u2 = (1.0 + 0.3 * mesh2.y ** 2, 0.1 * mesh2.x * mesh2.y)
    u3 = (1.0 + 0.3 * mesh3.y ** 2, 0.1 * mesh3.x * mesh3.y, np.zeros(mesh3.n_nodes))
    return mesh2, mesh3, u2, u3


@pytest.mark.parametrize("element_length", ["streamline", "metric"])
def test_sa_jacobian_matches_finite_differences_in_3d(element_length):
    _, mesh3, _, u3 = extruded_pair(2)
    sa = SpalartAllmarasSolver(mesh3, SpalartAllmaras(NU), ["bottom"],
                               element_length=element_length, backend="numpy")
    sa.set_velocity(*u3)
    rng = np.random.default_rng(1)
    for sign in (1.0, -1.0):
        nt = sign * NU * (0.5 + 5 * rng.random(mesh3.n_nodes))
        d = rng.random(mesh3.n_nodes) - 0.5
        eps = 1e-7 * NU
        for supg, tol in ((False, 1e-6), (True, 5e-3)):
            R, J = sa.residual_jacobian(nt, supg=supg)
            fd = (sa.residual_jacobian(nt + eps * d, supg=supg)[0]
                  - sa.residual_jacobian(nt - eps * d, supg=supg)[0]) / (2 * eps)
            assert np.linalg.norm(J @ d - fd) < tol * np.linalg.norm(fd)


@pytest.mark.parametrize("element_length", ["streamline", "metric"])
def test_sa_numba_matches_numpy_in_3d(element_length):
    pytest.importorskip("numba")
    _, mesh3, _, u3 = extruded_pair(2)
    solvers = [SpalartAllmarasSolver(mesh3, SpalartAllmaras(NU), ["bottom"],
                                     element_length=element_length, backend=b)
               for b in ("numpy", "numba")]
    rng = np.random.default_rng(3)
    nt = NU * (rng.random(mesh3.n_nodes) * 6 - 1.0)              # both signs
    for sa in solvers:
        sa.set_velocity(*u3)
    for supg in (False, True):
        (R0, J0), (R1, J1) = (s.residual_jacobian(nt, supg=supg) for s in solvers)
        assert np.abs(R0 - R1).max() < 1e-12 * max(1.0, np.abs(R0).max())
        assert abs(J0 - J1).max() < 1e-12 * max(1.0, abs(J0).max())


@pytest.mark.parametrize("element_length", ["streamline", "metric"])
def test_sa_residual_of_an_extruded_field_matches_2d(element_length):
    """On the extruded mesh the 3-D SA residual of a z-independent field equals the
    2-D residual scaled by the depth on every plane (mid-plane nodes carry the
    quadratic weight 4/6, the end planes 1/6 each)."""
    mesh2, mesh3, u2, u3 = extruded_pair(3)
    sa2 = SpalartAllmarasSolver(mesh2, SpalartAllmaras(NU), ["bottom"],
                                element_length=element_length)
    sa3 = SpalartAllmarasSolver(mesh3, SpalartAllmaras(NU), ["bottom"],
                                element_length=element_length)
    sa2.set_velocity(*u2)
    sa3.set_velocity(*u3)
    assert np.allclose(sa3.distance, np.tile(sa2.distance, 3))
    nt2 = NU * (1.0 + 2.0 * mesh2.y + 0.5 * np.sin(mesh2.x))
    nt3 = np.tile(nt2, 3)
    N = mesh2.n_nodes
    depth = 0.25
    # Galerkin terms: exact proportionality
    R2 = sa2.residual_jacobian(nt2, supg=False)[0]
    R3 = sa3.residual_jacobian(nt3, supg=False)[0]
    for plane, weight in ((0, 1 / 6), (1, 4 / 6), (2, 1 / 6)):
        assert np.allclose(R3[plane * N:(plane + 1) * N], weight * depth * R2, atol=1e-14)
    # with SUPG the stabilisation parameter sees the third dimension (the
    # hexahedral metric and the streamline length differ from the quad ones), so
    # the 3-D residual is close to, not identical with, the scaled 2-D one
    R2 = sa2.residual_jacobian(nt2)[0]
    R3 = sa3.residual_jacobian(nt3)[0]
    assert np.abs(R3[N:2 * N] - 4 / 6 * depth * R2).max() < 0.1 * np.abs(R2).max() * depth


def test_sa_steady_solve_extruded_equals_2d():
    mesh2 = geometry(2.0, 1.0, 5, "quad9")
    mesh3 = mesh2.extrude(0.3, layers=1)
    N = mesh2.n_nodes
    for supg, tol in ((False, 1e-8), (True, 2e-2)):      # SUPG parameters are mesh dependent
        results = []
        for mesh, vel in ((mesh2, (mesh2.y, np.zeros(N))),
                          (mesh3, (mesh3.y, np.zeros(3 * N), np.zeros(3 * N)))):
            sa = SpalartAllmarasSolver(mesh, SpalartAllmaras(NU), ["bottom"])
            sa.supg = supg
            sa.set_velocity(*vel)
            inflow, wall = mesh.bc_nodes["left"], mesh.bc_nodes["bottom"]
            fixed = np.unique(np.concatenate([inflow, wall]))
            vals = np.where(np.isin(fixed, wall), 0.0, 3 * NU)
            res = sa.solve(np.full(mesh.n_nodes, 3 * NU), fixed, vals, rtol=1e-8)
            assert res.converged
            results.append(res.T)
        nt2, nt3 = results
        for plane in range(3):
            assert np.abs(nt3[plane * N:(plane + 1) * N] - nt2).max() < tol * nt2.max()


def test_rans_coupling_extruded_equals_2d():
    """Two outer iterations of the SA coupling on a coarse channel: the 3-D
    extrusion with symmetry planes reproduces the 2-D eddy viscosity and wall force."""
    import warnings

    from aa540fem.turbulence import solve_rans

    mesh2 = geometry(2.0, 1.0, 4, "quad9")
    mesh3 = mesh2.extrude(0.5, layers=1)
    N = mesh2.n_nodes
    inflow2 = lambda x, y: 4 * y * (1 - y)
    inflow3 = lambda x, y, z: 4 * y * (1 - y)
    prob2 = FlowProblem(mesh2, mu=NU, rho=1.0, stabilisation=True,
                        bc={"left": (inflow2, 0.0), "top": (0.0, 0.0), "bottom": (0.0, 0.0),
                            "right": "open"})
    prob3 = FlowProblem(mesh3, mu=NU, rho=1.0, stabilisation=True,
                        bc={"left": (inflow3, 0.0, 0.0), "top": (0.0, 0.0, 0.0),
                            "bottom": (0.0, 0.0, 0.0), "front": (None, None, 0.0),
                            "back": (None, None, 0.0), "right": "open"})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r2 = solve_rans(prob2, wall_tags=["bottom"], max_outer=2)
        r3 = solve_rans(prob3, wall_tags=["bottom"], max_outer=2)
    # the SUPG parameters of both equations see the third dimension: agreement
    # to a fraction of a per cent, not to round-off
    scale = max(r2.nu_t.max(), NU)
    assert np.abs(r3.nu_t[:N] - r2.nu_t).max() < 1e-2 * scale
    assert np.abs(r3.nu_t[N:2 * N] - r2.nu_t).max() < 1e-2 * scale
    assert np.abs(r3.nu_t[N:2 * N] - r3.nu_t[:N]).max() < 1e-9 * scale   # z-independent
    fx2, _ = r2.flow.forces("bottom")
    fx3, _, fz3 = r3.flow.forces("bottom")
    assert fx3 == pytest.approx(0.5 * fx2, rel=1e-2) and abs(fz3) < 1e-10
    assert np.abs(r3.flow.w).max() < 1e-9
