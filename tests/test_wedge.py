"""Quadratic prisms: extrusion of mixed triangle6/quad9 meshes into wedge18 and
hexahedron27 cells, graded layers, boundary tags with two face types, the flow
and turbulence kernels on prisms, and Gmsh prism files."""

import itertools

import numpy as np
import pytest

from aa540fem import geometry
from aa540fem.core.elements import get_element
from aa540fem.core.mesh import Mesh, graded_layers
from aa540fem.incompressible import FlowProblem, solve_flow
from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.turbulence import SpalartAllmaras, SpalartAllmarasSolver


def mixed_channel(n=2, length=2.0):
    """Quads on the left half and triangles on the right half of ``[0, length] x [0, 1]``."""
    q = geometry(length / 2, 1.0, n, "quad9")
    t = geometry(length / 2, 1.0, n, "triangle6")
    pts = np.vstack([q.points, t.points + [length / 2, 0.0]])
    cells = {"quad9": q.cells["quad9"], "triangle6": t.cells["triangle6"] + q.n_nodes}
    bnd = {"left": q.boundary["left"], "right": t.boundary["right"] + q.n_nodes,
           "bottom": np.vstack([q.boundary["bottom"], t.boundary["bottom"] + q.n_nodes]),
           "top": np.vstack([q.boundary["top"], t.boundary["top"] + q.n_nodes])}
    _, first, inv = np.unique(np.round(pts, 12), axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    return Mesh(pts[first], {k: inv[v] for k, v in cells.items()},
                {k: inv[v] for k, v in bnd.items()})


def test_graded_layers():
    assert np.allclose(graded_layers(1.0, 4), [0, 0.25, 0.5, 0.75, 1.0])
    z = graded_layers(0.41, 5, 1.5)
    t = np.diff(z)
    assert z[0] == 0 and np.isclose(z[-1], 0.41) and (t > 0).all()
    assert np.allclose(t, t[::-1]) and np.allclose(t[1:3] / t[0:2], 1.5)
    with pytest.raises(ValueError):
        graded_layers(1.0, 0)


def test_mixed_extrusion_geometry():
    m2 = mixed_channel(2)
    assert set(m2.cells) == {"quad9", "triangle6"} and m2.check_orientation() == 0
    m3 = m2.extrude(0.5, layers=3, grading=1.5)
    assert set(m3.cells) == {"hexahedron27", "wedge18"}
    assert m3.n_nodes == 7 * m2.n_nodes and m3.n_elems == 3 * m2.n_elems
    assert all((d > 0).all() for d in m3.jacobian_at_centroids().values())
    # the z planes carry both face types, the extruded edges only quadrilaterals
    for tag in ("front", "back"):
        assert isinstance(m3.boundary[tag], dict)
        assert set(m3.boundary[tag]) == {"quad9", "triangle6"}
    assert m3.boundary["bottom"].shape == (3 * m2.boundary["bottom"].shape[0], 9)
    assert np.allclose(m3.z[m3.bc_nodes["front"]], 0.0)
    assert np.allclose(m3.z[m3.bc_nodes["back"]], 0.5)
    assert m3.bc_nodes["front"].size == m2.n_nodes
    faces = m3.boundary_faces()
    tagged = sum(b.shape[0] for tag in m3.tags for b in m3.face_blocks(tag))
    assert sum(f.shape[0] for f in faces.values()) == tagged
    # wall distance to the mixed planes, graded interfaces, explicit interfaces
    assert np.allclose(m3.wall_distance(["front", "back"]), np.minimum(m3.z, 0.5 - m3.z))
    assert np.allclose(np.unique(m3.z)[::2], graded_layers(0.5, 3, 1.5))
    m3b = m2.extrude(None, layers=[0.0, 0.1, 0.4])
    assert np.allclose(np.unique(m3b.z), [0, 0.05, 0.1, 0.25, 0.4])
    with pytest.raises(ValueError):
        m2.extrude(1.0, layers=[0.0, 0.5, 0.4])
    # volume from the quadrature of the constant 1 over both blocks
    from aa540fem.core.quadrature import quadrature_points
    from aa540fem.transport.element import jacobian_nd

    volume = 0.0
    for name, conn in m3.cells.items():
        el = get_element(name)
        coords, w = quadrature_points(el.family, el.full_order)
        _, dnat = el.shape_at(coords)
        hs, _ = jacobian_nd(tuple(m3.points[conn, k] for k in range(3)), dnat)
        volume += np.sum(w[None, :] * hs)
    assert abs(volume - 1.0) < 1e-12


@pytest.mark.parametrize("stabilisation", [False, True])
def test_poiseuille_is_exact_on_prisms_and_hexahedra(stabilisation):
    m3 = mixed_channel(2).extrude(0.5, layers=2, grading=1.3)
    mu = 0.1
    inflow = lambda x, y, z: 4 * y * (1 - y)
    prob = FlowProblem(m3, mu=mu, rho=1.0, stabilisation=stabilisation,
                       bc={"left": (inflow, 0.0, 0.0), "top": (0.0, 0.0, 0.0),
                           "bottom": (0.0, 0.0, 0.0), "front": (None, None, 0.0),
                           "back": (None, None, 0.0), "right": "open"})
    sol = solve_flow(prob, rtol=1e-12)
    assert np.abs(sol.u - 4 * m3.y * (1 - m3.y)).max() < 1e-12
    assert np.abs(sol.v).max() < 1e-12 and np.abs(sol.w).max() < 1e-12
    assert np.abs(sol.p_nodal - 8 * mu * (2.0 - m3.x)).max() < 1e-12
    # wall shear mu du/dy = 4 mu over the area 2 x 0.5 on the plate of both face types
    fx, fy, fz = sol.forces("bottom")
    assert np.isclose(fx, 4 * mu * 1.0) and abs(fz) < 1e-12
    # the z planes (mixed faces): no shear, pressure force in -z on the front plane
    fx, fy, fz = sol.forces("front")
    assert abs(fx) < 1e-12 and np.isclose(fz, -8 * mu * 2.0)       # -(mean p) x area
    tr = sol.wall_traction("bottom")
    assert np.allclose(tr["tx"], 4 * mu) and np.allclose(tr["ny"], -1.0)
    assert np.isclose(tr["weight"].sum(), 1.0)


def test_numba_matches_numpy_on_prisms():
    pytest.importorskip("numba")
    m3 = mixed_channel(2).extrude(0.5, layers=1)
    rng = np.random.default_rng(1)
    N = m3.n_nodes
    for stab, pspg, length in itertools.product([False, True], [False, True],
                                               ["metric", "streamline"]):
        prob = FlowProblem(m3, mu=1e-3, rho=1.2, stabilisation=stab, pspg=pspg,
                           element_length=length, eddy_viscosity=1e-3 * (1 + m3.y),
                           body_force=lambda x, y, z: (np.sin(x), np.cos(y), np.sin(z)),
                           bc={"left": (1.0, 0.0, 0.0), "bottom": (0.0, 0.0, 0.0)})
        ref, fast = FlowAssembler(prob, backend="numpy"), FlowAssembler(prob, backend="numba")
        U = np.zeros(ref.space.ndof)
        U[:N] = 1 - np.exp(-m3.y / 0.2) + 0.1 * np.sin(m3.x)
        U[N:2 * N] = 0.05 * np.sin(m3.y) * np.cos(m3.x)
        U[2 * N:3 * N] = 0.03 * np.sin(m3.z) * np.cos(m3.y)
        U[3 * N:] = np.cos(m3.x[ref.space.pressure_nodes])
        U_old = U + 0.05 * rng.standard_normal(U.size)
        for kw in ({}, {"dt": 0.05, "U_old": U_old, "param_state": U_old}, {"jacobian": False}):
            a, b = ref.momentum_terms(U, **kw), fast.momentum_terms(U, **kw)
            for x, y in ((a.N, b.N), (a.S, b.S), (a.J_N, b.J_N), (a.J_S, b.J_S)):
                if x is None:
                    assert y is None
                    continue
                assert abs(x - y).max() <= 1e-12 * max(1.0, abs(y).max())


@pytest.mark.parametrize("element_length", ["metric", "streamline"])
def test_spalart_allmaras_on_prisms(element_length):
    m3 = mixed_channel(2).extrude(0.5, layers=1)
    rng = np.random.default_rng(2)
    vel = (1 + 0.3 * m3.y ** 2, 0.1 * m3.x * m3.y, 0.05 * m3.z)
    sa = SpalartAllmarasSolver(m3, SpalartAllmaras(1e-3), ["bottom"],
                               element_length=element_length, backend="numpy")
    sa.set_velocity(*vel)
    assert np.allclose(sa.distance, m3.y)
    nt = 1e-3 * (rng.random(m3.n_nodes) * 6 - 1)
    d = rng.random(m3.n_nodes) - 0.5
    eps = 1e-10
    for supg, tol in ((False, 1e-6), (True, 5e-3)):
        R, J = sa.residual_jacobian(nt, supg=supg)
        fd = (sa.residual_jacobian(nt + eps * d, supg=supg)[0]
              - sa.residual_jacobian(nt - eps * d, supg=supg)[0]) / (2 * eps)
        assert np.linalg.norm(J @ d - fd) < tol * np.linalg.norm(fd)
    pytest.importorskip("numba")
    fast = SpalartAllmarasSolver(m3, SpalartAllmaras(1e-3), ["bottom"],
                                 element_length=element_length, backend="numba")
    fast.set_velocity(*vel)
    for supg in (False, True):
        (R0, J0), (R1, J1) = (s.residual_jacobian(nt, supg=supg) for s in (sa, fast))
        assert np.abs(R0 - R1).max() < 1e-11 * max(1.0, np.abs(R0).max())
        assert abs(J0 - J1).max() < 1e-11 * max(1.0, abs(J0).max())


def test_gmsh_prism_is_renumbered_to_vtk_order(tmp_path):
    gmsh = pytest.importorskip("gmsh")
    pytest.importorskip("meshio")
    from aa540fem import read_mesh

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("prism")
        geo = gmsh.model.geo
        p = [geo.addPoint(*c, 1.0) for c in [(0, 0, 0), (1, 0, 0), (0, 1, 0)]]
        lines = [geo.addLine(p[0], p[1]), geo.addLine(p[1], p[2]), geo.addLine(p[2], p[0])]
        surface = geo.addPlaneSurface([geo.addCurveLoop(lines)])
        geo.extrude([(2, surface)], 0, 0, 2.0, [1], recombine=True)
        geo.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", 5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 5)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(2)
        path = tmp_path / "prism.msh"
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    mesh = read_mesh(path)
    assert set(mesh.cells) == {"wedge18"} and mesh.n_elems == 1 and mesh.n_nodes == 18
    el = get_element("wedge18")
    conn, P = mesh.cells["wedge18"], mesh.points
    # a straight prism: mid-edge nodes and face centres sit where the element says
    for i, (a, b) in enumerate(el.edges):
        assert np.allclose(P[conn[:, 6 + i]], 0.5 * (P[conn[:, a]] + P[conn[:, b]]))
    for face, fel in el.face_elements():
        if fel.n_nodes == 9:
            assert np.allclose(P[conn[:, face[8]]], P[conn[:, list(face[:4])]].mean(axis=1))
    assert (mesh.jacobian_at_centroids()["wedge18"] > 0).all()
    # the extruded surfaces come along untagged: one boundary with both face types
    assert mesh.tags == ["boundary"] and set(mesh.boundary["boundary"]) == {"triangle6", "quad9"}
    assert mesh.face_blocks("boundary")[0].shape[0] + mesh.face_blocks("boundary")[1].shape[0] == 5


def test_prism_in_the_registry():
    el = get_element("wedge18")
    assert el.dim == 3 and el.n_corners == 6 and el.n_faces == 5
    assert el.face_types.count("quad9") == 3 and el.face_types.count("triangle6") == 2
    from aa540fem.core.elements import PRESSURE_ELEMENT

    assert PRESSURE_ELEMENT["wedge18"].name == "wedge"
