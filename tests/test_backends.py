"""The numba assembly kernels reproduce the NumPy reference kernels."""

import itertools

import numpy as np
import pytest

from aa540fem import geometry, read_mesh
from aa540fem.backends import assembly_backend, assembly_backends
from aa540fem.incompressible import FlowAssembler, FlowProblem

MESHES = __import__("pathlib").Path(__file__).resolve().parent.parent / "examples" / "meshes"
numba = pytest.importorskip("numba")


def meshes():
    yield "quad9", geometry(2.0, 1.0, 3, "quad9")
    yield "triangle6", geometry(2.0, 1.0, 3, "triangle6")
    try:
        yield "cylinder_bl", read_mesh(MESHES / "cylinder_bl.msh")
    except ImportError:                                   # meshio missing
        pass


def states(asm, rng):
    N, ndof = asm.space.N, asm.space.ndof
    mesh = asm.mesh
    U = np.zeros(ndof)
    U[:N] = 1 - np.exp(-mesh.y / 0.2) + 0.1 * np.sin(mesh.x)
    U[N:2 * N] = 0.05 * np.sin(mesh.y) * np.cos(mesh.x)
    U[2 * N:] = np.cos(mesh.x[asm.space.pressure_nodes])
    U_old = U + 0.05 * rng.standard_normal(ndof)
    par = U + 0.1 * rng.standard_normal(ndof)
    return U, U_old, par


def compare(a, b, tol):
    if a is None or b is None:
        assert a is None and b is None
        return
    if hasattr(a, "toarray"):                              # sparse: never densify
        scale = max(1.0, abs(b).max())
        diff = a - b
        assert (abs(diff).max() if diff.nnz else 0.0) <= tol * scale
        return
    a, b = np.asarray(a), np.asarray(b)
    scale = max(1.0, np.abs(b).max())
    assert np.abs(a - b).max() <= tol * scale


OPTIONS = list(itertools.product([False, True], [False, True], ["metric", "streamline"]))


@pytest.mark.parametrize("name,mesh", list(meshes()))
@pytest.mark.parametrize("stabilisation,pspg,element_length", OPTIONS)
def test_numba_matches_numpy(name, mesh, stabilisation, pspg, element_length):
    rng = np.random.default_rng(1)
    tags = list(mesh.boundary)
    for eddy, body in ((None, None),
                       (1e-3 * (1 + mesh.y), lambda x, y: (np.sin(x), np.cos(y)))):
        prob = FlowProblem(mesh, mu=1e-3, rho=1.2, stabilisation=stabilisation, pspg=pspg,
                           element_length=element_length, eddy_viscosity=eddy, body_force=body,
                           bc={tags[0]: (1.0, 0.0), tags[1]: (0.0, 0.0)})
        ref = FlowAssembler(prob, backend="numpy")
        fast = FlowAssembler(prob, backend="numba")
        assert fast.backend == "numba"
        for A, B in ((ref.K, fast.K), (ref.M, fast.M), (ref.B, fast.B), (ref.BT, fast.BT)):
            compare(A, B, 1e-13)
        U, U_old, par = states(ref, rng)
        for kw in ({}, {"dt": 0.05, "U_old": U_old}, {"param_state": par},
                   {"jacobian": False}, {"dt": 0.05, "U_old": U_old, "param_state": par}):
            a = ref.momentum_terms(U, **kw)
            b = fast.momentum_terms(U, **kw)
            compare(a.N, b.N, 1e-13)
            compare(a.S, b.S, 1e-12)
            compare(a.J_N, b.J_N, 1e-13)
            compare(a.J_S, b.J_S, 1e-12)
        compare(ref.body_load(), fast.body_load(), 1e-14)


def test_thread_count_does_not_change_the_result():
    mesh = read_mesh(MESHES / "cylinder_bl.msh") if (MESHES / "cylinder_bl.msh").exists() \
        else geometry(2.0, 1.0, 6, "quad9")
    tags = list(mesh.boundary)
    prob = FlowProblem(mesh, mu=1e-3, rho=1.0, stabilisation=True,
                       bc={tags[0]: (1.0, 0.0), tags[1]: (0.0, 0.0)})
    asm = FlowAssembler(prob, backend="numba")
    U = np.linspace(-1, 1, asm.space.ndof)
    old = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        one = asm.momentum_terms(U)
        numba.set_num_threads(max(2, old))
        many = asm.momentum_terms(U)
    finally:
        numba.set_num_threads(old)
    assert np.array_equal(one.N, many.N) and np.array_equal(one.S, many.S)
    assert np.array_equal(one.JN_data, many.JN_data) and np.array_equal(one.JS_data, many.JS_data)


def test_backend_selection(monkeypatch):
    assert assembly_backends()[0] == "numba"
    monkeypatch.setenv("AA540FEM_ASSEMBLY", "numpy")
    assert assembly_backend() == "numpy"
    assert assembly_backend("numba") == "numba"
    with pytest.raises(ValueError):
        assembly_backend("cuda")


@pytest.mark.parametrize("element_length", ["metric", "streamline"])
def test_spalart_allmaras_kernel_matches_numpy(element_length):
    from aa540fem.turbulence.spalart_allmaras import SpalartAllmaras, SpalartAllmarasSolver

    mesh = geometry(2.0, 1.0, 4, "quad9")
    nu = 1e-4
    solvers = [SpalartAllmarasSolver(mesh, SpalartAllmaras(nu), ["bottom"],
                                     element_length=element_length, backend=b)
               for b in ("numpy", "numba")]
    for s in solvers:
        s.set_velocity(1 - np.exp(-mesh.y / 0.2), 0.1 * np.sin(mesh.x) * mesh.y)
    rng = np.random.default_rng(2)
    nt = 3 * nu * (1 + mesh.y) + 0.5 * nu * rng.standard_normal(mesh.n_nodes)   # some negative
    nt[mesh.bc_nodes["bottom"]] = 0.0
    for supg in (False, True):
        R0, J0 = solvers[0].residual_jacobian(nt, supg=supg)
        R1, J1 = solvers[1].residual_jacobian(nt, supg=supg)
        compare(R0, R1, 1e-11)
        compare(J0, J1, 1e-10)
    compare(solvers[0].M, solvers[1].M, 1e-14)
