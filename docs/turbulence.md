# Turbulence modelling: the Spalart-Allmaras model

This note documents the Reynolds-averaged turbulence model implemented in
`src/aa540fem/turbulence/`: the equations, the constants, how the model is
discretised and coupled to the incompressible solver, how it is validated,
and its limitations.  Section references point to the sources listed at the
end.

## 1. Reynolds averaging and the eddy-viscosity hypothesis

The Reynolds-averaged Navier-Stokes (RANS) equations for a steady
incompressible mean flow `u` with pressure `p` read

    rho (u . grad) u = - grad p + div[ (mu + mu_t) (grad u + grad u^T) ],   div u = 0,

where the Reynolds stresses have been modelled with the Boussinesq
hypothesis `-rho <u'_i u'_j> = mu_t (d_i u_j + d_j u_i) - (2/3) rho k delta_ij`
(the isotropic part is absorbed in the pressure).  Everything about the
turbulence is therefore condensed into one field, the eddy viscosity
`mu_t = rho nu_t`, and a turbulence model is a recipe for `nu_t` [1, 2].

## 2. The Spalart-Allmaras model

Spalart and Allmaras [3] proposed a single transport equation for a
working variable `nu_tilde`, built term by term for aerodynamic boundary
layers and wakes and calibrated on the flat plate and mixing layer.  The
form implemented here is the "SA-neg" variant of Allmaras, Johnson and
Spalart [4], as written on the NASA Turbulence Modeling Resource [5], with
the trip term `f_t2` switched off ("SA-neg-noft2", the standard choice for
fully turbulent computations).

### 2.1 Transport equation (`nu_tilde >= 0`)

    d(nu_tilde)/dt + u . grad(nu_tilde)
        = c_b1 (1 - f_t2) S_tilde nu_tilde                          (production)
          - [ c_w1 f_w - (c_b1 / kappa^2) f_t2 ] (nu_tilde / d)^2      (destruction)
          + (1/sigma) [ div((nu + nu_tilde) grad nu_tilde)
                        + c_b2 |grad nu_tilde|^2 ]                    (diffusion)

with `d` the distance to the nearest wall, `nu` the laminar kinematic
viscosity, and

    nu_t = nu_tilde f_v1,        f_v1 = chi^3 / (chi^3 + c_v1^3),   chi = nu_tilde / nu,
    S_tilde = Omega + (nu_tilde / (kappa^2 d^2)) f_v2,               f_v2 = 1 - chi / (1 + chi f_v1),
    f_w = g [ (1 + c_w3^6) / (g^6 + c_w3^6) ]^(1/6),   g = r + c_w2 (r^6 - r),
    r = min( nu_tilde / (S_tilde kappa^2 d^2), 10 ),
    f_t2 = c_t3 exp(-c_t4 chi^2)   (off by default).

`Omega = |d_x u_y - d_y u_x|` is the magnitude of the vorticity.  The
`f_v1`, `f_v2` functions make `nu_t` vanish like `y^4` at the wall so that
the model is integrated to the wall without wall functions (which requires
the first cell at `y+` of order one).

The modified vorticity is guarded so that `S_tilde` never becomes negative
or small [4, 5]: with `Sbar = nu_tilde f_v2 / (kappa^2 d^2)`,

    S_tilde = Omega + Sbar                                           if Sbar >= -c_v2 Omega,
    S_tilde = Omega + Omega (c_v2^2 Omega + c_v3 Sbar) / ((c_v3 - 2 c_v2) Omega - Sbar)   otherwise.

### 2.2 Negative branch (`nu_tilde < 0`)

The original model can produce negative `nu_tilde` transients that are not
well posed.  The negative modification [4] replaces, for `nu_tilde < 0`,

    production   = c_b1 (1 - c_t3) Omega nu_tilde,
    destruction  = - c_w1 (nu_tilde / d)^2,
    diffusivity  = (nu + f_n nu_tilde) / sigma,   f_n = (c_n1 + chi^3) / (c_n1 - chi^3),
    nu_t         = 0,

which keeps the equation dissipative and drives the variable back to zero.

### 2.3 Constants

| constant | value | | constant | value |
|---|---|---|---|---|
| `c_b1` | 0.1355 | | `c_v1` | 7.1 |
| `sigma` | 2/3 | | `c_t3`, `c_t4` | 1.2, 0.5 |
| `c_b2` | 0.622 | | `c_v2`, `c_v3` | 0.7, 0.9 |
| `kappa` | 0.41 | | `c_n1` | 16 |
| `c_w1` | `c_b1/kappa^2 + (1 + c_b2)/sigma` = 3.2391 | | `r` limit | 10 |
| `c_w2`, `c_w3` | 0.3, 2 | | | |

All are in `SpalartAllmaras` (`spalart_allmaras.py`) as dataclass fields.

### 2.4 Boundary conditions

- Walls: `nu_tilde = 0`.
- Inflow and far field: `nu_tilde = 3 nu`, the value recommended by the
  TMR [5] (it gives a freestream `nu_t / nu` of about 0.21); imposed on
  every boundary that has a full velocity Dirichlet condition and is not a
  wall.
- Outflow and symmetry lines: natural (zero diffusive flux).

## 3. Implementation

### 3.1 Where the pieces live

| file | contents |
|---|---|
| `turbulence/spalart_allmaras.py` | `SpalartAllmaras` (model functions of section 2), `SpalartAllmarasSolver` (assembly of the steady equation and its solution) |
| `turbulence/rans.py` | `solve_rans`: the segregated flow / turbulence iteration, `RANSSolution` |
| `core/wall_distance.py` | nodal distance to tagged walls |
| `incompressible/assembler.py` | variable viscosity `mu + mu_t(x)` in the momentum equations |
| `linalg/continuation.py` | pseudo-transient continuation shared by the flow and the SA solve |
| `examples/turbulent_flat_plate.py`, `examples/meshes/flat_plate_turb.msh` | validation case |
| `tests/test_turbulence.py` | unit tests |

### 3.2 Wall distance

`Mesh.wall_distance(tags)` computes for every node the minimum distance to
the boundary faces of the wall tags: in 2-D each half of a quadratic edge
is a straight segment (exact on straight walls; the error on a curved wall
is the sagitta of a half edge, checked on the annulus in the tests), in
3-D each quadratic triangle is split into four and each quadrilateral into
eight flat sub-triangles and the point-triangle distance (Ericson's
region test, `core/wall_distance.py`) is taken.  The computation is
vectorised in chunks of nodes; for the meshes here it takes seconds.  At the quadrature points
the distance is interpolated from the corner nodes with the linear shape
functions, a convex combination that stays between the nodal values: the
quadratic interpolant undershoots to negative values in the distorted cells
around the ends of a wall, and since the distance enters as
`(nu_tilde / d)^2` that produced destruction terms of 1e11 at the start of
the plate case.  (The floor of 1e-12 remains as a guard; `nu_tilde = 0` on
the wall nodes.)

### 3.3 Discretisation of the SA equation

The working variable lives on the same nodes and quadratic elements as
the velocity (`quad9`, `triangle6`), so it reuses the precomputed
quadrature blocks of the flow solver (`incompressible/space._Block`: shape
functions, physical gradients, weights).  For a frozen velocity field
`(u, v)` the steady weak form assembled by `SpalartAllmarasSolver.residual_jacobian` is

    R_i = int phi_i (u . grad nu_tilde)                      convection (Galerkin)
        + int D(nu_tilde) grad phi_i . grad nu_tilde          diffusion, D = (nu + f_n nu_tilde) / sigma
        - int phi_i s(nu_tilde, grad nu_tilde; Omega, d)      production - destruction + (c_b2/sigma)|grad nu_tilde|^2
        + int tau (u . grad phi_i) [ u . grad nu_tilde - s ]  SUPG stabilisation.

The model source `s` is evaluated at the quadrature points from the
interpolated `nu_tilde`, its gradient, the vorticity of the velocity field
(from the velocity gradients) and the wall distance.  The equation is
convection dominated everywhere except in the near-wall region, so the
convective term is stabilised with the same streamline-upwind
Petrov-Galerkin weighting and cell measure used for the momentum equations:
`tau = [u . G u + D^2 G:G / 2]^(-1/2)` with the element metric tensor `G`
(section 3.4; `u . G u = (2|u|/h_s)^2` for the cell size `h_s` in the flow
direction, `D^2 G:G / 2 = (4 D / h^2)^2` on a square cell and the smallest
cell dimension on a stretched one [6, 9, 10]); the diffusion part of the
strong residual is omitted in the SUPG term, the usual simplification.

### 3.4 Newton linearisation

The Jacobian is analytic for the convection and diffusion operators, and
uses central finite differences at the quadrature points for the model
functions: `ds/d(nu_tilde)`, `ds/d(grad nu_tilde)` and `dD/d(nu_tilde)`.
This avoids differentiating `f_w`, `S_tilde` and the negative branch by
hand (error prone in every published implementation) at the cost of a few
extra vectorised evaluations per assembly.  With SUPG off the Jacobian
matches a finite-difference directional derivative to 1e-8 (test); with
SUPG on the stabilisation parameter of the turbulence equation is frozen,
which leaves a small inconsistency that slows Newton but does not change
the converged solution.

The momentum equations do differentiate their stabilisation parameters,
and the way the cell size enters them decided whether the wall-resolved
plate could be solved at all.  Tezduyar's flow-direction element length
`h = 2 / sum_i |s . grad phi_i|`, `s = u/|u|` [6], is the natural choice on
isotropic meshes, but on a boundary-layer cell of aspect ratio 1000 it
jumps from the cell length to the cell height when the velocity rotates by
a milliradian, and its derivative changes sign at every such rotation; the
Newton direction then stops being a descent direction long before the
residual is small, and the pseudo-transient continuation crawls at one per
cent per step.  The element metric tensor `G = J^-T T J^-1` (Shakib [9],
Bazilevs et al. [10]; `J` the Jacobian of the isoparametric map, `T` a
constant that refers triangles to an equilateral reference element) gives
the smooth alternative

    tau   = [ (2/dt)^2 + u . G u + nu^2 G:G / 2 ]^(-1/2)
    gamma = (h_s |u| / 2) min(1, Re_h / 3),   h_s |u| / 2 = |u|^2 / sqrt(u . G u)

with the same limits as before (`u . G u = (2|u|/h_s)^2`; the viscous term
equals `(4 nu/h^2)^2` on a square cell and uses the smallest cell dimension
on a stretched one, as Tezduyar's `h_RGN` does in boundary layers).
`FlowProblem(element_length="metric")` is the default; `"streamline"`
keeps the previous definition.  On the plate mesh at Re 1e4 the first flow
solve of the RANS start-up converges in 9 pseudo-time steps with the
metric form (the step grows from 0.25 to 3e4 cell CFL numbers) where the
streamline form was still at a residual of 2e-3 after 100 steps.

### 3.5 Steady solve: pseudo-transient continuation

Starting from the freestream level, the SA equation is solved by
pseudo-transient continuation (`linalg/continuation.pseudo_transient`,
shared with the flow solver): backward-Euler pseudo-time steps
`R(x) + M (x - x_k) / dtau = 0` with the consistent mass matrix `M`, up to a
few Newton iterations per step, and the step grown by switched evolution
relaxation `dtau_{k+1} = dtau_k |R_k| / |R_{k+1}|` [7].  Near the wall the
destruction term `(nu_tilde/d)^2` makes the equation stiff; the implicit
pseudo-time treatment handles it, and the iteration turns into plain Newton
once the residual is small.

### 3.6 Coupling to the flow: `solve_rans`

The coupling is segregated:

1. solve the flow with the current eddy viscosity (initially from the
   freestream `nu_tilde`, i.e. essentially laminar) using
   `solve_flow(continuation="auto")`, which falls back to pseudo-transient
   continuation when Newton from rest diverges at high Reynolds number;
2. solve the SA equation for the frozen velocity field;
3. under-relax `nu_tilde` (factor 0.7), update `mu_t = rho nu_tilde f_v1`
   and re-solve the flow from the previous state (Newton converges in a few
   iterations from there);
4. stop when the relative change of `nu_t` between outer iterations drops
   below the tolerance (1e-3 by default).

The flow assembler takes the nodal `mu_t` through
`FlowProblem.eddy_viscosity`.  The viscous operator becomes
`int (mu + mu_t) grad u : grad v`, and because the viscosity now varies the
divergence of the full stress `div[mu_eff (grad u + grad u^T)]` differs from
the Laplacian form by `grad(u)^T . grad(mu_eff)` (using `div u = 0`); this
term is kept as a volume term (assembled into `K` with the gradient of the
nodal `mu_t`), which preserves the do-nothing outflow condition
`mu du/dn - p n = 0` of the Laplacian form.  The strong momentum residual
used by the stabilisation carries the same terms, and the finite-difference
Jacobian check of the momentum equations with a spatially varying `mu_t`
passes to 1e-10 (test).

### 3.7 Start-up at high Reynolds number

Starting the wall-resolved plate at `Re = 1e6` from uniform flow does not
work: the first cells are 2e-5 thick, the impulsive start puts a wall shear
of order `1 / 2e-5` into the residual and neither Newton nor the
pseudo-transient continuation recovers from it.  `solve_rans` therefore
offers three devices, all used by the validation case:

- **a smooth initial profile** (`U0` may be a callable `(x, y) -> (u, v)`):
  the example starts from `u = 1 - exp(-y / delta(x))` with a boundary-layer
  thickness `delta` of a few per cent of the plate length, which already has
  the right wall shear to within an order of magnitude;
- **a viscosity ramp** (`viscosity_ramp=(100, 10, 1)`): the coupled problem
  is first converged loosely at 100 and 10 times the laminar viscosity, each
  stage starting from the previous flow and `nu_tilde`, before the target
  viscosity is solved; the eddy viscosity of the previous stage is a good
  guess for the next because `nu_t / nu` in the log layer depends only
  weakly on the Reynolds number;
- **local pseudo-time stepping** in both continuations
  (`solve_flow(local_timestep=True)`, `SpalartAllmarasSolver.solve(local_timestep=True)`):
  the pseudo-time step of every node is scaled by its convective time scale
  `h / U` with `h` the longest edge of the surrounding cells
  (`Mesh.nodal_size("max")`), so `dtau0` is a CFL number and the wall cells
  and the coarse far field advance at their own pace.  The choice of `h`
  matters on stretched cells: with a global step the far field never
  moves, and with the thin dimension of the wall cells (their viscous
  scale) the corrections at the leading edge need thousands of steps to
  convect along the plate, so the residual crawled at one per cent per
  step; with the streamwise dimension the second ramp stage converges in
  11 pseudo-time steps.  The implicit steps do not need the viscous limit.

The pseudo-transient continuation of the flow also projects its starting
velocity onto the discretely divergence-free space (the saddle-point
projection of the RK45 integrator, `steady.project_divergence_free`).  From
a velocity that violates continuity, such as the impulsive start or a
damped Newton iterate, the first pseudo-time step needs a pressure jump of
order `1 / dtau` to enforce it, and the stabilisation terms, quadratic in
velocity and pressure, turn that jump into a residual that does not shrink
with the step: every step is rejected.  With the projection the laminar
plate at Re 1e5 converges from rest in 10 pseudo-time steps.  In the
`"auto"` mode the continuation starts from Newton's last iterate only if
Newton halved the residual; a couple of heavily damped Newton steps that
barely reduced it leave a worse state than the initial one (on the second
ramp stage the first pseudo-time step from them is rejected, while from the
initial state the stage converges in 11 steps).

The sub-solves inside the outer iteration are converged only to a relative
residual of 1e-5 (flow) and 1e-4 (turbulence): the outer iteration changes
both fields again anyway, and the final flow solve tightens automatically
once `nu_t` has settled.

## 4. Validation: turbulent flat plate

`examples/turbulent_flat_plate.py`: zero-pressure-gradient plate of length
2 in a `[-0.5, 2.5] x [0, 1]` domain, `U = 1`, `nu = 1e-6` (Re per unit
length 1e6, `Re_x` up to 2e6), inflow on the left and top, symmetry lines
upstream and downstream of the plate, do-nothing outflow, `nu_tilde_inf = 3 nu`.
The mesh (`make_meshes.make_turbulent_flat_plate()`) has 30 quadrilateral
layers from the wall, the first one 2e-5 thick (`y+` about 1 with the computed
wall shear), inside a triangular mesh; the streamwise spacing shrinks from
0.01 in the middle of the plate to 1e-4 at its two ends; 23k nodes.  The
refinement at the ends is essential: the leading and trailing edges are
singular points of the flow (the pressure and the wall shear behave like
`x^(-1/2)`), and with cells a thousand times longer than high there the
Newton and pseudo-transient solvers stalled on the final ramp stage, with
the whole residual concentrated in the half-dozen cells around the two
points.  A first cell of aspect ratio 50 at the ends removes the problem.

Comparisons:

- skin friction against the correlation of White [8, eq. 6-78]
  `Cf = 0.455 / ln^2(0.06 Re_x)` (accurate to a few per cent for
  1e5 < Re_x < 1e9) and the 1/7-power law `Cf = 0.0576 Re_x^(-1/5)`;
- the law of the wall at `x = 1.5`: `u+ = y+` in the viscous sublayer and
  `u+ = ln(y+)/kappa + B` with `kappa = 0.41`, `B = 5.0` [1, 8] in the log
  layer; the model was calibrated to reproduce it, so this is the sharpest
  check of the near-wall implementation.

### Results

Start-up with the viscosity ramp 100, 10, 1: the three stages take 5, 6 and
11 outer iterations; every flow solve after the first converges with plain
Newton in 2 to 5 iterations and the whole case runs in 4 minutes with the
numba kernels and the MUMPS factorisation (9 minutes with NumPy and
SuperLU; 23k nodes, one 4-core machine).  The final eddy viscosity reaches `nu_t / nu` of
about 200 in the wake at the outlet.

Skin friction along the plate (window 3e5 < Re_x < 2e6, i.e. the last 85 %
of the plate; the two nodes at the trailing edge, where the wall shear is
singular, are excluded):

| Re_x   | Cf, FEM + SA | White [8] | 1/7-power law |
|--------|--------------|-----------|---------------|
| 3.0e5  | 0.00413      | 0.00473   | 0.00462       |
| 5.0e5  | 0.00384      | 0.00428   | 0.00417       |
| 1.0e6  | 0.00347      | 0.00376   | 0.00363       |
| 1.5e6  | 0.00327      | 0.00350   | 0.00335       |
| 1.9e6  | 0.00318      | 0.00336   | 0.00320       |

Over the window the computed skin friction lies 4.5 to 12.8 % below White's
correlation (mean 7.7 %) and 0 to 10.6 % below the 1/7-power law (mean
4.3 %), the deficit shrinking downstream.  Both correlations assume a
boundary layer that is turbulent from the leading edge, while the model
without trip term starts laminar: the computed `Cf` follows Blasius up to
`Re_x` of about 1e4 and rises to the turbulent level between 1e4 and 4e4 as
`nu_tilde` grows from its freestream value, so the momentum thickness at a
given `x` is not that of the correlations.  The comparison that removes
this dependence on the origin is the Coles-Fernholz law at the momentum
thickness Reynolds number of the computed profile: at `x = 1.5` the profile
gives `Re_theta = 2884` and `Cf = 0.00327` against the Coles-Fernholz
value `2 [ln(Re_theta)/0.384 + 4.127]^(-2) = 0.00323` (Nagib, Chauhan and
Monkewitz [11]), a difference of 1.1 %.

Law of the wall at `x = 1.5` (`u_tau = 0.0404`, first node at `y+ = 0.40`,
62 nodes across the layer): the viscous sublayer sits on `u+ = y+`; in the
log layer `30 < y+ < 300` the computed `u+` exceeds `ln(y+)/0.41 + 5.0` by
0.15 to 0.87 (20 nodes, i.e. within 4 % of `u+`), the excess growing towards
the outer edge where the wake component sets in; in the buffer layer
`5 < y+ < 30` the profile runs up to 1.3 below the two-piece law, which is
the crudeness of the law there rather than of the solution.  The freestream
is reached at `u+ = 24.9` (`delta_99 = 0.026`).

The figure `turbulent_plate.png` written by the example shows both
comparisons.

### 3.8 Three dimensions

The solver and its numba kernel take the dimension from the mesh: the
gradients are tuples of `d` arrays, the vorticity magnitude is `|curl u|`
in 3-D (`|dv/dx - du/dy|` in 2-D), the SUPG parameter uses the `d x d`
metric tensor or the streamline length in `d` dimensions, and
`set_velocity` takes the `d` components (`solve_rans` passes
`flow.velocity`).  The consistency checks (`tests/test_turbulence3d.py`)
extrude a 2-D quadratic mesh into one layer of 27-node hexahedra and
18-node prisms (`Mesh.extrude`) with symmetry planes: the Galerkin part of the 3-D
residual of a z-independent field is exactly the 2-D residual times the
depth (weights 1/6, 4/6, 1/6 on the three node planes), the steady SA
solve without SUPG reproduces the 2-D field to 1e-8, and with SUPG the
3-D solution and the coupled RANS eddy viscosity and wall force agree with
2-D to a fraction of a per cent (the stabilisation parameters see the third
mesh dimension, so exact equality is not expected).  The numba kernel
matches the NumPy one to 1e-12 in 3-D.

The wall-resolved plate itself runs in 3-D with
`examples/turbulent_flat_plate.py --extrude 0.1`: the mesh is extruded one
layer (prisms from the triangles, hexahedra from the boundary-layer
quadrilaterals, symmetry planes `u_z = 0`) and the 3-D coupling starts
from the converged 2-D velocity, pressure and `nu_tilde` copied onto the
node planes (`solve_rans(..., U0=, nu_tilde0=)`, no viscosity ramp).  The
full-length mesh gives 280k unknowns, which the direct solver cannot
factorise in 15 GB, and the fieldsplit Krylov solver stalls on the `y+ = 1`
cells (in 2-D as well), while coarsening the streamwise spacing makes the
Re 1e6 case diverge; the half-length twin `flat_plate_turb_short.msh`
(plate `0 <= x <= 1`, the same cells, 13k nodes, 125k unknowns in 3-D,
4.3 GB) is the case that fits.  Its 3-D solution converges in one outer
iteration (274 s) and gives, at `x = 0.75` (`Re_theta = 1639`),
`Cf = 0.00364` against the Coles-Fernholz value 0.00365 (-0.3 %), the
first node at `y+ = 0.43`, and `u+` within 0.2 to 1.8 of the log law for
`30 < y+ < 300`; the 2-D solution on the same mesh gives `Cf = 0.00364`
at `Re_theta = 1638` and `u+` within 0.25 to 1.80 of the log law, i.e. the
3-D discretisation reproduces the 2-D one to the printed digits (the
prism and hexahedral SUPG parameters differ from the 2-D ones by the
extrusion depth, see the consistency tests above).  The SA solves of the coupling stop at 1e-4 of the
residual of the freestream start with the current velocity (an absolute
target next to the relative one): a start near the solution was otherwise
asked for 1e-4 of an already small residual and ran out of pseudo-time
steps.

## 4b. The lifting airfoil: why alpha = 10 does not converge

`examples/turbulent_airfoil.py` reproduces the NASA TMR NACA 0012 case at
Re 6e6 at alpha = 0 (Cd 0.00741 against CFL3D 0.00819, Cp within about
0.01, see `README.md`) and fails to converge at alpha = 10.  The diagnosis
is recorded here because every part of it is a property of the
discretisation rather than of one solver setting.

1. **The residual norm measures the far field, not the boundary layer.**
   The mesh spans 100 chords and its first cell is 9.1e-7 chords, so the
   element volumes differ by more than ten orders of magnitude.  The
   validated alpha = 0 state has |R| = 96.6, of which the cells inside
   d < 1e-2 contribute about 1e-4.  At incidence the far-field circulation
   adjustment dominates |R| completely.
2. **Therefore relative tolerances short-circuit the lifting case.**  The
   far-field part of the residual is nearly linear, so one near-Newton step
   removes it and a relative test declares convergence after a single
   iteration.  With `local_timestep=True` the pseudo-time step of a
   boundary-layer node is scaled by `h/u_ref` and advances about 3e-4 per
   step, so nothing develops there by marching either: a cold start at
   alpha = 10 converges to CL 0.949, Cd 0.0524 with Cd_f = 5e-5, i.e. the
   first cells still hold their initial profile.  The signature to watch
   for is Cd_f near zero with flow sub-solves accepting in one iteration.
3. **Absolute tolerances are equally useless**, for the same reason from
   the other side: no absolute level distinguishes a converged boundary
   layer from a frozen one, and asking for one (atol 3e-5 on this mesh)
   simply burns the whole step budget every stage.
4. **What limits the solve is the stabilisation linearisation.**  With a
   fresh Jacobian per iteration, a backward-Euler step from the alpha = 0
   state goes |R| 99 -> 1.01 -> 0.93 and then stops with no descent
   direction.  The residual left over lives entirely in the boundary layer
   (0.82 within d < 1e-4, 2e-5 outside d > 1e-2, pressure 1e-14).  This is
   the frozen SUPG tau of section 3.4 in the flow equations: tau ~ h/|u| is
   least consistent exactly where h is smallest and |u| varies steepest, so
   the Newton direction stops being a descent direction while the boundary
   layer is still unconverged.  Freezing the Jacobian across iterations
   makes it worse (the step stalls at |R| ~ 60 instead of 0.93).
5. **Consequences for continuation.**  Viscosity continuation converges the
   upper rungs honestly if every stage is given a step budget rather than a
   tolerance (factor 100, 30 and 10 reach |R| ~ 1e-5, CL 1.011 at Re 6e5),
   but no state below a factor of about 3 is marchable: the post-SA state
   diverges to |R| 3.1e4 at every pseudo-time step size tried.  Nor is the
   ramp robust: walking down in 1.25x steps instead of 3.3x never reaches
   the fine rungs, because on a regenerated mesh differing by 40 nodes
   (gmsh is not bitwise reproducible) the factor-30 stage itself diverges,
   at the same outer iteration and to the same digits on repeated runs,
   where it had converged on the original mesh.  The signal common to all
   of these is the change in eddy viscosity between outer iterations
   staying of order one: the segregated coupling has no stable fixed point
   in the transitional regime on this mesh, which no tolerance or guard
   threshold can repair.  Angle
   continuation from the converged alpha = 0 solution fails differently -
   a rotation-blended restart introduces kinks at the wall-distance medial
   axis, and a smooth time-ramped far-field direction (a time-dependent
   Dirichlet condition, which the theta scheme supports) clears its first
   step and then stalls on the same boundary-layer floor.  Global
   pseudo-time stepping is not an escape either: `dtau` is then limited by
   the 9.1e-7 cells and collapses to 1e-9.

The two changes that would unblock it are a consistent linearisation of
the stabilisation terms (differentiating tau with respect to the velocity)
and a residual measure normalised by cell volume, so that convergence is
judged where the physics is.  Both are solver changes rather than tuning,
which is why this release documents the case instead of reporting a number
for it.

## 5. Limitations and next steps

- Steady RANS; the model is used without the trip term, i.e. fully
  turbulent from the leading edge (no transition prediction).  The 3-D
  validation is the extruded plate (section 3.8), a 2-D flow solved with
  the 3-D code; a genuinely three-dimensional turbulent reference case
  (a wing-body junction, a swept plate) is still to be run.
- Wall-resolved only: the first cell must be at `y+` of order one.  Wall
  functions would allow coarser meshes and are a boundary-condition
  addition in `SpalartAllmarasSolver`.
- The turbulence equation freezes its SUPG parameter in the Jacobian, and
  the coupling is segregated with under-relaxation; a monolithic Newton on
  flow plus turbulence would converge faster but needs the
  cross-derivatives.
- The direct solver limits the mesh to roughly a hundred thousand
  unknowns in 3-D; wall-resolved meshes at flight Reynolds numbers need
  the fieldsplit Krylov solver (`docs/parallel.md`) and the MPI domain
  decomposition.
- The rotation/curvature correction ("SA-RC") and the quadratic
  constitutive relation ("SA-QCR"), useful for wing-body junctions and
  vortices, are not implemented.

## References

1. S. B. Pope, *Turbulent Flows*, Cambridge University Press, 2000.
2. D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries, 2006.
3. P. R. Spalart and S. R. Allmaras, "A one-equation turbulence model for
   aerodynamic flows", *La Recherche Aerospatiale* 1 (1994) 5-21; also
   AIAA Paper 92-0439, 1992.
4. S. R. Allmaras, F. T. Johnson and P. R. Spalart, "Modifications and
   clarifications for the implementation of the Spalart-Allmaras
   turbulence model", *7th International Conference on Computational Fluid
   Dynamics*, ICCFD7-1902, Big Island, Hawaii, 2012.
5. C. L. Rumsey (ed.), *NASA Langley Turbulence Modeling Resource*,
   https://turbmodels.larc.nasa.gov (Spalart-Allmaras page and the
   "2D Zero Pressure Gradient Flat Plate" verification case).
6. T. E. Tezduyar, "Stabilized finite element formulations for
   incompressible flow computations", *Advances in Applied Mechanics* 28
   (1992) 1-44.
7. C. T. Kelley and D. E. Keyes, "Convergence analysis of pseudo-transient
   continuation", *SIAM Journal on Numerical Analysis* 35 (1998) 508-523.
8. F. M. White, *Viscous Fluid Flow*, 3rd ed., McGraw-Hill, 2006 (chapter
   6: the turbulent flat-plate skin-friction correlations and the law of
   the wall).
9. F. Shakib, T. J. R. Hughes and Z. Johan, "A new finite element
   formulation for computational fluid dynamics: X. The compressible Euler
   and Navier-Stokes equations", *Computer Methods in Applied Mechanics and
   Engineering* 89 (1991) 141-219 (the element metric tensor in the
   stabilisation parameter).
10. Y. Bazilevs, V. M. Calo, J. A. Cottrell, T. J. R. Hughes, A. Reali and
    G. Scovazzi, "Variational multiscale residual-based turbulence modeling
    for large eddy simulation of incompressible flows", *Computer Methods in
    Applied Mechanics and Engineering* 197 (2007) 173-201 (the metric form
    `tau = [4/dt^2 + u.G u + C_I nu^2 G:G]^(-1/2)`).
11. H. M. Nagib, K. A. Chauhan and P. A. Monkewitz, "Approach to an
    asymptotic state for zero pressure gradient turbulent boundary layers",
    *Philosophical Transactions of the Royal Society A* 365 (2007)
    755-770 (the Coles-Fernholz skin-friction law with kappa = 0.384,
    C = 4.127).
