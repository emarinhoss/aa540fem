# Python port of the AA 540 anisotropic heat conduction FEM code

This directory is a NumPy/SciPy port of the MATLAB code in `../project`.
It solves, on the rectangle `[0, a] x [0, b]`,

    div( kappa(x, y) . grad T ) + f(x, y) = 0

where `kappa` is a full 2x2 conductivity tensor, with Dirichlet
(`T = T0`) or Neumann (`n . kappa grad T = q_n`) conditions on each side,
using 3-node linear triangles, 4-node bilinear quadrilaterals or 9-node
biquadratic quadrilaterals.

## Running

```bash
pip install -r requirements.txt
python main.py                  # writes temperature.png
python main.py --show           # interactive window
python -m pytest tests          # verification against analytical solutions
```

As in the MATLAB version, the two files meant to be edited are `main.py`
(domain size, mesh, element type, boundary conditions, quadrature order) and
`aa540fem/conductivity_and_forcing.py` (conductivity tensor and heat source).

The solver can also be driven from Python:

```python
from aa540fem import Problem, solve

problem = Problem(a=4, b=6, elems=50, elem_type=3, order=2,
                  bc_type={"top": 0, "right": 1, "left": 1, "bottom": 0},
                  bc_val={"top": 100, "right": 0, "left": 0, "bottom": 0})
solution = solve(problem)
solution.T            # nodal temperatures
solution.grid         # temperatures on the (ny, nx) grid of solution.mesh.X / .Y
```

Boundary values and the `material` callable may be functions of `(x, y)`.

## Layout

| MATLAB file                  | Python module                          |
|------------------------------|----------------------------------------|
| `main.m`                     | `main.py` (inputs) + `aa540fem/solver.py` |
| `geometry.m`                 | `aa540fem/geometry.py`                 |
| `interpfunc_{3,4,9}.m`       | `aa540fem/shape_functions.py`          |
| `gauss_legendre_quad.m`, `gauss_trgl.m` | `aa540fem/quadrature.py`    |
| `elem_eqn.m`                 | `aa540fem/element.py`                  |
| `dirichlet.m`                | `aa540fem/boundary.py`                 |
| `conductivity_and_forcing.m` | `aa540fem/conductivity_and_forcing.py` |

## Differences from the MATLAB code

The port keeps the structure and the user interface of the original but
fixes several problems found while translating it:

* **Shape-function derivatives.** `interpfunc_4.m` had the wrong sign on
  half of the derivatives and `interpfunc_9.m` had a wrong `dPhi/dxi` for
  node 1.  The Python versions are built as tensor products of 1-D Lagrange
  polynomials and are checked against finite differences in the tests.
* **Quadrature.** The MATLAB code only sampled the diagonal `xi = eta` of
  the quadrilateral and used the 1-D rule for triangles too.  The port uses
  the full tensor-product Gauss-Legendre rule on quadrilaterals and the
  FSELIB `gauss_trgl` rule on triangles (`order` selects the rule).  The
  MATLAB default `order = 1` under-integrates 4- and 9-node elements and
  makes the system singular, so `order = None` now selects the lowest
  order that fully integrates the element (1, 2 or 3) and a warning is
  issued for anything lower.
* **9-node mesh.** The connectivity was hard-coded for a 100-element mesh
  and elements overlapped; it is now derived from the mesh size.  Nodes are
  numbered with `x` varying fastest for every element type.
* **`elems` for triangles** is the number of cells per side (each split in
  two triangles), the same meaning as for quadrilaterals.
* **Neumann conditions** are implemented (the MATLAB branch was empty), as
  the boundary integral of the prescribed flux.  A zero flux, the default,
  reproduces the MATLAB behaviour.
* The global system is stored as a SciPy sparse matrix, so large 9-node
  meshes fit in memory.
