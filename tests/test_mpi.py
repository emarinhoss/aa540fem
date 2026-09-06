"""Partitioning (serial) and the MPI prototypes.

The MPI tests need ``mpi4py`` and ``petsc4py`` and run on several ranks:

    mpirun -n 4 python -m pytest tests/test_mpi.py -q

In a serial process they run on one rank (the distributed code paths with
``size == 1``), so the ordinary test suite covers them too.
"""

import numpy as np
import pytest

from aa540fem import geometry, read_mesh
from aa540fem.incompressible import FlowAssembler, FlowProblem
from aa540fem.parallel.partition import element_parts, partition_mesh

MESHES = __import__("pathlib").Path(__file__).resolve().parent.parent / "examples" / "meshes"


def cavity(elems=6):
    mesh = geometry(1.0, 1.0, elems, "quad9")
    walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
    return FlowProblem(mesh, mu=0.05, rho=1.0, bc={**walls, "top": (1.0, 0.0)},
                       stabilisation=True)


@pytest.mark.parametrize("n_parts", [1, 3])
def test_partition_covers_the_mesh(n_parts):
    pytest.importorskip("meshio")
    mesh = read_mesh(MESHES / "cylinder_bl.msh")               # mixed quad9 / triangle6
    parts = partition_mesh(mesh, n_parts)
    assert len(parts) == n_parts
    per_block = {name: np.concatenate([p.elements[name] for p in parts])
                 for name in mesh.cells}
    for name, conn in mesh.cells.items():
        assert np.array_equal(np.sort(per_block[name]), np.arange(conn.shape[0]))
    owned = np.concatenate([p.owned for p in parts])
    assert np.array_equal(np.sort(owned), np.arange(mesh.n_nodes))      # disjoint, complete
    for p in parts:
        assert np.all(p.owner[p.owned] == p.part) and np.all(p.owner[p.ghost] != p.part)
        assert np.array_equal(np.sort(p.local[p.local >= 0]), np.arange(p.n_local))
        for name, local_conn in p.local_cells(mesh).items():
            assert local_conn.min() >= 0 and local_conn.max() < p.n_local
    if n_parts > 1:
        sizes = [sum(idx.size for idx in p.elements.values()) for p in parts]
        assert max(sizes) < 2 * min(sizes)                                 # balanced


def test_fallback_partition_without_metis(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def no_pymetis(name, *args, **kwargs):
        if name == "pymetis":
            raise ImportError
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_pymetis)
    mesh = geometry(2.0, 1.0, 4, "quad9")
    parts = element_parts(mesh, 2)
    assert set(parts) == {0, 1} and abs((parts == 0).sum() - (parts == 1).sum()) <= 1


def test_distributed_factorisation_matches_serial():
    pytest.importorskip("mpi4py")
    from aa540fem.linalg.direct import petsc_available

    if not petsc_available():
        pytest.skip("petsc4py not available")
    from aa540fem.linalg.dirichlet import eliminate
    from aa540fem.parallel.comm import world
    from aa540fem.parallel.distributed import DistributedFactorisation

    asm = FlowAssembler(cavity())
    fixed, vals = asm.dirichlet()
    U = np.zeros(asm.space.ndof)
    U[fixed] = vals
    R, J = asm.steady_residual_jacobian(U, asm.body_load())
    elim = eliminate(J, fixed)
    A, b = elim.K_bc, elim.apply_rhs(-R, np.zeros(fixed.size))
    f = DistributedFactorisation(A, comm=world())
    x = f.solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-10
    # every rank holds the same full solution
    comm = world()
    ref = comm.bcast(x, root=0)
    assert np.allclose(x, ref)


def test_distributed_assembly_matches_serial():
    pytest.importorskip("mpi4py")
    from aa540fem.linalg.direct import petsc_available

    if not petsc_available():
        pytest.skip("petsc4py not available")
    from aa540fem.backends.numpy_kernels import linear_local
    from aa540fem.core.mesh import Mesh
    from aa540fem.parallel.comm import world
    from aa540fem.parallel.distributed import distributed_matrix, gather_csr

    comm = world()
    prob = cavity(8)
    mesh = prob.mesh
    serial = FlowAssembler(prob)                                  # the reference, on every rank
    parts = partition_mesh(mesh, comm.Get_size())
    part = parts[comm.Get_rank()]
    # a sub-mesh of the part's elements in local numbering, assembled with the serial kernels
    nodes = np.concatenate([part.owned, part.ghost])
    local_cells = part.local_cells(mesh)
    sub = Mesh(mesh.points[nodes], local_cells, {})
    sub_prob = FlowProblem(sub, mu=prob.mu, rho=prob.rho, bc={})
    sub_asm = FlowAssembler(sub_prob)
    values, dofs = [], []
    for b in sub_asm.blocks:
        mu_q, dmu = sub_asm.viscosity(b)
        K_e, M_e, B_e, BT_e = linear_local(b, mu_q, dmu, prob.rho, False)
        values.append(K_e + M_e + sum(B_e) - BT_e)
        gnodes = nodes[b.conn]                                    # local -> global nodes
        dofs.append(serial.space.local_dofs(gnodes, gnodes[:, :b.nc]))
    values = np.concatenate(values) if values else np.zeros((0, 1, 1))
    dofs = np.concatenate(dofs) if dofs else np.zeros((0, 1), dtype=int)
    A = distributed_matrix(serial.space.ndof, dofs, values, comm=comm)
    gathered = gather_csr(A, comm)
    ref = serial.K + serial.M + serial.Bx + serial.By - serial.BT
    assert abs(gathered - ref).max() < 1e-12 * abs(ref).max()


def _petsc_or_skip():
    pytest.importorskip("mpi4py")
    from aa540fem.linalg.direct import petsc_available

    if not petsc_available():
        pytest.skip("petsc4py not available")


def cylinder():
    pytest.importorskip("meshio")
    mesh = read_mesh(MESHES / "cylinder_bl.msh")
    inflow = lambda x, y: 4.0 * 0.3 * y * (0.41 - y) / 0.41 ** 2
    return FlowProblem(mesh, mu=1e-3, rho=1.0, stabilisation=True,
                       bc={"inlet": (inflow, 0.0), "walls": (0.0, 0.0), "cylinder": (0.0, 0.0),
                           "outlet": "open"})


@pytest.mark.parametrize("case,method", [("cavity", "direct"), ("cavity", "fieldsplit"),
                                         ("cylinder", "direct"), ("cylinder", "fieldsplit")])
def test_distributed_newton_matches_serial(case, method):
    """Every rank assembles its own elements; the distributed Newton reproduces the
    serial solution (same iterations, states equal to the solver tolerance)."""
    _petsc_or_skip()
    from aa540fem.incompressible import solve_flow
    from aa540fem.parallel.flow import DistributedFlowSystem

    prob = cavity(10) if case == "cavity" else cylinder()
    ref = solve_flow(prob)
    system = DistributedFlowSystem(prob)
    if system.size > 1:
        assert system.part.owned.size < prob.mesh.n_nodes            # really distributed
    sol = system.solve_steady(method=method)
    assert sol.info["converged"]
    assert abs(sol.info["iterations"] - ref.info["iterations"]) <= 1
    assert np.abs(sol.U - ref.U).max() < 1e-7 * np.abs(ref.U).max()
    if method == "fieldsplit":
        assert max(sol.info["linear_iterations"]) < 200


def test_distributed_theta_scheme_matches_serial():
    _petsc_or_skip()
    from aa540fem.incompressible import solve_flow_transient
    from aa540fem.parallel.flow import DistributedFlowSystem

    prob = cavity(8)
    ref = solve_flow_transient(prob, dt=0.05, t_end=0.2, scheme="theta", startup_steps=1,
                               rtol=1e-8)
    run = DistributedFlowSystem(prob).solve_transient(0.05, 0.2, startup_steps=1, rtol=1e-8)
    assert run.info["factorisations"] <= ref.info["factorisations"]
    assert np.abs(run.snapshots[-1] - ref.snapshots[-1]).max() < 1e-7


def test_distributed_state_round_trip():
    """set_state / gather are inverse on every rank and the Dirichlet values are imposed."""
    _petsc_or_skip()
    from aa540fem.parallel.flow import DistributedFlowSystem

    prob = cavity(6)
    system = DistributedFlowSystem(prob)
    rng = np.random.default_rng(0)
    U = rng.standard_normal(system.space.ndof)
    system.set_state(U)
    back = system.gather()
    fixed, vals = FlowAssembler(prob).dirichlet()
    assert np.allclose(back[fixed], vals)
    free = np.ones(U.size, dtype=bool)
    free[fixed] = False
    assert np.allclose(back[free], U[free])


def test_distributed_newton_in_3d_matches_serial():
    _petsc_or_skip()
    from aa540fem.core.mesh import box
    from aa540fem.incompressible import solve_flow
    from aa540fem.parallel.flow import DistributedFlowSystem

    mesh = box(2.0, 1.0, 1.0, (3, 2, 2), "tetra10")
    inflow = lambda x, y, z: 4.0 * y * (1.0 - y)
    prob = FlowProblem(mesh, mu=0.05, rho=1.0, stabilisation=True,
                       bc={"left": (inflow, 0.0, 0.0), "bottom": (0.0, 0.0, 0.0),
                           "top": (0.0, 0.0, 0.0), "front": (None, None, 0.0),
                           "back": (None, None, 0.0), "right": "open"})
    ref = solve_flow(prob)
    sol = DistributedFlowSystem(prob).solve_steady()
    assert sol.info["converged"]
    assert np.abs(sol.U - ref.U).max() < 1e-8 * np.abs(ref.U).max()
    assert np.abs(sol.u - 4.0 * mesh.y * (1.0 - mesh.y)).max() < 1e-8
