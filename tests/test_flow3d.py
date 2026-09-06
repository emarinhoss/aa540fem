"""Three-dimensional Navier-Stokes: exact extruded flows, Jacobians, forces, time stepping."""

import numpy as np
import pytest

from aa540fem.core.mesh import box
from aa540fem.incompressible import FlowAssembler, FlowProblem, solve_flow, solve_flow_transient
from aa540fem.linalg.dirichlet import eliminate

MU = 0.1


def channel(elem_type, elems=(2, 2, 1), stabilisation=False, umax=1.0, size=(2.0, 1.0, 0.5)):
    """Plane Poiseuille flow between the walls y = 0 and y = 1, extruded in z:
    parabolic inflow, no-slip walls, symmetry planes (only u_z fixed) in z and
    an open outlet.  The exact solution is quadratic, so Taylor-Hood
    reproduces it to round-off."""
    mesh = box(*size, elems, elem_type)
    inflow = lambda x, y, z: 4.0 * umax * y * (1.0 - y)
    return FlowProblem(mesh, mu=MU, rho=1.0, stabilisation=stabilisation,
                       bc={"left": (inflow, 0.0, 0.0), "bottom": (0.0, 0.0, 0.0),
                           "top": (0.0, 0.0, 0.0), "front": (None, None, 0.0),
                           "back": (None, None, 0.0), "right": "open"})


def exact(mesh, umax=1.0, length=2.0):
    u = 4.0 * umax * mesh.y * (1.0 - mesh.y)
    p = 8.0 * MU * umax * (length - mesh.x)               # dp/dx = -8 mu umax, p(outlet) = 0
    return u, p


@pytest.mark.parametrize("elem_type", ["hexahedron27", "tetra10"])
@pytest.mark.parametrize("stabilisation", [False, True])
def test_extruded_poiseuille_is_exact(elem_type, stabilisation):
    prob = channel(elem_type, stabilisation=stabilisation)
    sol = solve_flow(prob)
    u, p = exact(prob.mesh)
    assert sol.info["converged"]
    tol = 1e-9 if stabilisation else 1e-12
    assert np.abs(sol.u - u).max() < tol
    assert np.abs(sol.v).max() < tol and np.abs(sol.w).max() < tol
    assert np.abs(sol.p_nodal - p).max() < tol
    assert sol.divergence_norm() < 1e-10
    assert sol.speed.max() == pytest.approx(1.0)


def test_exact_solution_holds_on_a_refined_tetrahedral_channel():
    prob = channel("tetra10", elems=(4, 3, 2), stabilisation=True, size=(2.0, 1.0, 1.0))
    sol = solve_flow(prob)
    u, p = exact(prob.mesh)
    assert sol.info["converged"] and sol.info["iterations"] <= 6
    assert np.abs(sol.u - u).max() < 1e-8 and np.abs(sol.p_nodal - p).max() < 1e-7


@pytest.mark.parametrize("elem_type", ["hexahedron27", "tetra10"])
def test_wall_traction_and_forces_of_the_channel(elem_type):
    prob = channel(elem_type, elems=(2, 2, 2))
    sol = solve_flow(prob)
    # wall shear on the bottom wall: mu du/dy = 4 mu, in +x (fluid drags the wall along)
    tr = sol.wall_traction("bottom")
    assert np.allclose(tr["tx"], 4.0 * MU, atol=1e-9)
    assert np.allclose(tr["ty"], -tr["p"], atol=1e-9)      # pressure pushes the wall down
    assert np.allclose(tr["nz"], 0.0) and np.allclose(tr["ny"], -1.0)
    fx, fy, fz = sol.forces("bottom")
    area = 2.0 * 0.5                                        # wall length x depth
    assert fx == pytest.approx(4.0 * MU * area, rel=1e-9)
    assert fz == pytest.approx(0.0, abs=1e-9)
    assert tr["weight"].sum() == pytest.approx(area)
    # the two walls carry equal and opposite shear
    assert sol.forces("top")[0] == pytest.approx(fx, rel=1e-9)


@pytest.mark.parametrize("elem_type", ["hexahedron27", "tetra10"])
def test_jacobian_matches_finite_differences_in_3d(elem_type):
    """Stabilised (SUPG + PSPG + grad-div, metric parameters) Jacobian against
    central differences of the residual at a non-trivial state."""
    prob = channel(elem_type, elems=(1, 1, 1))
    prob.stabilisation = True
    prob.pspg = True
    asm = FlowAssembler(prob, backend="numpy")
    rng = np.random.default_rng(0)
    U = rng.standard_normal(asm.space.ndof)
    F = asm.body_load()
    R, J = asm.steady_residual_jacobian(U, F)
    h = 1e-6
    for k in rng.choice(asm.space.ndof, 12, replace=False):
        e = np.zeros(asm.space.ndof)
        e[k] = h
        fd = (asm.steady_residual(U + e, F) - asm.steady_residual(U - e, F)) / (2 * h)
        col = J[:, k].toarray().ravel()
        assert np.abs(fd - col).max() < 1e-5 * max(1.0, np.abs(col).max())


def test_dirichlet_and_pin_in_3d():
    mesh = box(1.0, 1.0, 1.0, 1)
    walls = {t: (0.0, 0.0, 0.0) for t in ("left", "right", "bottom", "front", "back")}
    prob = FlowProblem(mesh, mu=0.01, rho=1.0, bc={**walls, "top": (1.0, 0.0, 0.0)},
                       stabilisation=True)
    assert prob.pins_pressure
    asm = FlowAssembler(prob)
    fixed, vals = asm.dirichlet()
    assert fixed.max() >= asm.space.n_vel                   # the pressure pin
    assert np.sum(vals == 1.0) == mesh.bc_nodes["top"].size
    with pytest.raises(ValueError, match="3 velocity components"):
        FlowProblem(mesh, bc={"top": (1.0, 0.0)}).validate()
    sol = solve_flow(prob)
    assert sol.info["converged"] and sol.divergence_norm() < 1e-8


def test_theta_and_rk45_reach_the_steady_channel():
    prob = channel("hexahedron27", elems=(2, 2, 1), stabilisation=True)
    u, _ = exact(prob.mesh)
    U0 = lambda x, y, z: (0.0 * x, 0.0 * y, 0.0 * z)
    run = solve_flow_transient(prob, dt=0.5, t_end=6.0, scheme="theta", theta=1.0, U0=U0,
                               rtol=1e-8)
    assert np.abs(run.final.u - u).max() < 1e-6
    rk = solve_flow_transient(prob, dt=0.02, t_end=0.2, scheme="rk45", U0=U0,
                              output_interval=0.1, adaptive=False)
    assert rk.times[-1] == pytest.approx(0.2) and np.isfinite(rk.final.U).all()
    # the projection keeps the discrete continuity equation B u = 0 exactly
    asm = FlowAssembler(prob)
    assert np.abs((asm.B @ rk.final.U)[asm.space.n_vel:]).max() < 1e-10


def test_fieldsplit_solver_in_3d():
    pytest.importorskip("petsc4py")
    from aa540fem.linalg.direct import petsc_available

    if not petsc_available():
        pytest.skip("petsc4py not available")
    prob = channel("tetra10", elems=(2, 2, 1), stabilisation=True)
    direct = solve_flow(prob)
    krylov = solve_flow(prob, method="fieldsplit")
    assert krylov.info["converged"]
    assert np.abs(krylov.U - direct.U).max() < 1e-6 * np.abs(direct.U).max()


def test_eliminated_system_size_and_save(tmp_path):
    pytest.importorskip("meshio")
    prob = channel("hexahedron27", elems=(1, 1, 1))
    asm = FlowAssembler(prob)
    fixed, _ = asm.dirichlet()
    _, J = asm.steady_residual_jacobian(np.zeros(asm.space.ndof), asm.body_load())
    A = eliminate(J, fixed).K_bc
    assert A.shape == (asm.space.ndof, asm.space.ndof) == (3 * 27 + 8, 3 * 27 + 8)
    sol = solve_flow(prob)
    out = sol.save(tmp_path / "channel.vtu")
    import meshio

    m = meshio.read(out)
    assert m.points.shape == (27, 3) and m.point_data["velocity"].shape == (27, 3)
