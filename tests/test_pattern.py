"""Fixed sparsity pattern, scatter plans and the pattern-based Dirichlet elimination."""

import numpy as np
import pytest
import scipy.sparse as sp

from aa540fem import geometry, read_mesh
from aa540fem.backends import ScatterPlan, SparsityPattern, scatter_vector
from aa540fem.backends.pattern import add_matrices
from aa540fem.incompressible import FlowAssembler, FlowProblem
from aa540fem.linalg.dirichlet import DirichletEliminator, PatternDirichlet, eliminate

MESHES = __import__("pathlib").Path(__file__).resolve().parent.parent / "examples" / "meshes"


def problem(mesh):
    return FlowProblem(mesh, mu=1e-2, rho=1.0, stabilisation=True,
                       bc={tag: (0.0, 0.0) for tag in list(mesh.boundary)[:2]})


@pytest.mark.parametrize("mesh_name", ["structured", "cylinder_bl"])
def test_scatter_reproduces_coo_assembly(mesh_name):
    if mesh_name == "structured":
        mesh = geometry(2.0, 1.0, 4, "quad9")
    else:
        pytest.importorskip("meshio")
        mesh = read_mesh(MESHES / "cylinder_bl.msh")          # mixed quad9 / triangle6
    asm = FlowAssembler(problem(mesh))
    rng = np.random.default_rng(0)
    vals = [rng.random((b.conn.shape[0], b.L, b.L)) for b in asm.blocks]
    data = asm.plan.assemble(vals)
    rows = np.concatenate([np.repeat(b.ldof, b.L, axis=1).ravel() for b in asm.blocks])
    cols = np.concatenate([np.tile(b.ldof, (1, b.L)).ravel() for b in asm.blocks])
    ref = sp.coo_matrix((np.concatenate([v.ravel() for v in vals]), (rows, cols)),
                        shape=(asm.space.ndof,) * 2).tocsr()
    assert abs(asm.pattern.matrix(data) - ref).max() < 1e-14 * abs(ref).max()
    # the assembled system matrices live on the pattern too
    for A in (asm.K, asm.M, asm.Bx, asm.By, asm.BT):
        assert np.allclose(asm.pattern.data_of(A), A.data)
    assert abs(asm.BT - (asm.Bx + asm.By).T).max() == 0.0


def test_scatter_map_rejects_entries_outside_the_pattern():
    pattern = SparsityPattern.from_element_dofs(4, [np.array([[0, 1], [2, 3]])])
    assert pattern.nnz == 8
    with pytest.raises(ValueError):
        pattern.scatter_map(np.array([[0]]), np.array([[3]]))
    with pytest.raises(ValueError):
        pattern.matrix(np.zeros(3))


def test_scatter_vector_matches_add_at():
    dofs = np.array([[0, 1, 2], [2, 3, 0]])
    vals = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    ref = np.zeros(5)
    np.add.at(ref, dofs.ravel(), vals.ravel())
    assert np.array_equal(scatter_vector(5, dofs, vals), ref)


def test_pattern_dirichlet_matches_generic_elimination():
    mesh = geometry(1.0, 1.0, 3, "quad9")
    asm = FlowAssembler(problem(mesh))
    fixed, vals = asm.dirichlet()
    vals = vals + np.linspace(0.5, 1.5, fixed.size)
    J = asm.pattern.matrix(asm.K_data - asm.BT_data + asm.B_data + 0.3 * asm.M_data)
    generic = DirichletEliminator(J, fixed)
    fast = eliminate(J, fixed)
    assert isinstance(fast, PatternDirichlet)
    assert abs(fast.K_bc - generic.K_bc).max() == 0.0
    F = np.linspace(-1, 1, asm.space.ndof)
    assert np.allclose(fast.apply_rhs(F, vals), generic.apply_rhs(F, vals), atol=1e-14)
    # the cache is per (pattern, nodes) and does not leak between node sets
    other = eliminate(J, fixed[::2])
    assert abs(other.K_bc - DirichletEliminator(J, fixed[::2]).K_bc).max() == 0.0


def test_add_matrices_uses_the_shared_pattern():
    mesh = geometry(1.0, 1.0, 2, "triangle6")
    asm = FlowAssembler(problem(mesh))
    A = add_matrices(asm.K, asm.M, 0.25)
    assert getattr(A, "pattern", None) is asm.pattern
    assert abs(A - (asm.K + 0.25 * asm.M)).max() == 0.0
    plain = add_matrices(asm.K, sp.eye(asm.space.ndof, format="csr"), 2.0)
    assert getattr(plain, "pattern", None) is None


def test_plan_partial_blocks():
    mesh = geometry(1.0, 1.0, 2, "quad9")
    asm = FlowAssembler(problem(mesh))
    b = asm.blocks[0]
    ne, n = b.conn.shape[0], b.n
    vals = np.ones((ne, 2 * n, 2 * n))
    data = asm.plan_uu.assemble([vals])
    A = asm.pattern.matrix(data)
    assert A[:2 * asm.space.N, :2 * asm.space.N].sum() == pytest.approx(vals.sum())
    assert abs(A[2 * asm.space.N:, :]).sum() == 0.0      # pressure rows untouched
    plan = ScatterPlan(asm.pattern, [asm.plan.maps[0]])
    with pytest.raises(ValueError):
        plan.assemble([np.ones((ne, 3, 3))])
