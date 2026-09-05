"""Post-processing: gradients, fluxes and error norms."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aa540fem.core.elements import get_element
from aa540fem.core.mesh import Mesh
from aa540fem.core.util import call_xyt
from aa540fem.transport.element import physical_gradients
from aa540fem.transport.material import conductivity_and_forcing as _default_material


@dataclass
class CellField:
    """Per-element quantities evaluated at the element centroids."""

    centroid: np.ndarray    # (n_elems, 2)
    grad: np.ndarray        # (n_elems, 2)  grad T
    flux: np.ndarray        # (n_elems, 2)  q = -kappa . grad T


def element_gradient(mesh: Mesh, T, material=None, t: float = 0.0) -> CellField:
    """Gradient and heat flux of the nodal field ``T`` at element centroids.

    Elements are concatenated in ``mesh.cells`` block order, matching the
    cell ordering written by :func:`aa540fem.mesh_io.write_vtk`.
    """
    if material is None:
        material = _default_material
    T = np.asarray(T, dtype=float)
    cents, grads, fluxes = [], [], []
    for name, conn in mesh.cells.items():
        el = get_element(name)
        phi, dxi, deta = el.shape([el.centroid[0]], [el.centroid[1]])
        xe, ye = mesh.x[conn], mesh.y[conn]
        _, dphi_dx, dphi_dy = physical_gradients(xe, ye, dxi, deta)
        Te = T[conn]
        gx = np.einsum("ei,ei->e", dphi_dx[:, 0, :], Te)
        gy = np.einsum("ei,ei->e", dphi_dy[:, 0, :], Te)
        X = xe @ phi[0]
        Y = ye @ phi[0]
        k11, k12, k21, k22, _ = (np.broadcast_to(np.asarray(v, dtype=float), X.shape)
                                 for v in call_xyt(material, X, Y, t))
        cents.append(np.column_stack([X, Y]))
        grads.append(np.column_stack([gx, gy]))
        fluxes.append(-np.column_stack([k11 * gx + k12 * gy, k21 * gx + k22 * gy]))
    return CellField(np.vstack(cents), np.vstack(grads), np.vstack(fluxes))


def _accurate_order(el):
    return 7 if el.family == "triangle" else el.min_order + 2


def error_norms(mesh: Mesh, T, exact, exact_grad=None, order=None) -> dict:
    """L2 and H1-seminorm errors of the nodal field ``T`` against ``exact(x, y)``.

    ``exact_grad(x, y)`` should return ``(dT/dx, dT/dy)``; without it the H1
    error is not computed.  Returns ``{"L2": ..., "H1": ...}`` with absolute
    errors integrated by a quadrature of higher order than the assembly.
    """
    T = np.asarray(T, dtype=float)
    l2 = 0.0
    h1 = 0.0
    for name, conn in mesh.cells.items():
        el = get_element(name)
        xi, eta, w = el.quadrature(_accurate_order(el) if order is None else order)
        phi, dxi, deta = el.shape(xi, eta)
        xe, ye = mesh.x[conn], mesh.y[conn]
        hs, dphi_dx, dphi_dy = physical_gradients(xe, ye, dxi, deta)
        X = xe @ phi.T
        Y = ye @ phi.T
        Te = T[conn]
        Th = Te @ phi.T                                   # (n_elems, nq)
        wh = w[None, :] * hs
        l2 += np.sum(wh * (Th - exact(X, Y)) ** 2)
        if exact_grad is not None:
            gx = np.einsum("eqi,ei->eq", dphi_dx, Te)
            gy = np.einsum("eqi,ei->eq", dphi_dy, Te)
            ex, ey = exact_grad(X, Y)
            h1 += np.sum(wh * ((gx - ex) ** 2 + (gy - ey) ** 2))
    out = {"L2": float(np.sqrt(l2))}
    if exact_grad is not None:
        out["H1"] = float(np.sqrt(h1))
    return out
