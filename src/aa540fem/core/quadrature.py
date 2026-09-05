"""Gauss quadrature rules for line, quadrilateral and triangle elements."""

from __future__ import annotations

import numpy as np

_TRGL_ORDERS = (1, 3, 4, 6, 7, 9, 12, 13)

# Lowest quadrature order that integrates the stiffness matrix of each
# element type without rank deficiency (for constant conductivity):
# 1 point for linear triangles, 2x2 for 4-node and 3x3 for 9-node quads.
MIN_ORDER = {1: 1, 2: 2, 3: 3}


def default_order(elem_type: int) -> int:
    """Default quadrature order for an element type (see ``MIN_ORDER``)."""
    try:
        return MIN_ORDER[elem_type]
    except KeyError:
        raise ValueError(f"Unknown element type {elem_type!r}; expected 1, 2 or 3") from None


def gauss_legendre_quad(order: int) -> tuple[np.ndarray, np.ndarray]:
    """1-D Gauss-Legendre abscissas and weights on [-1, 1].

    Port of ``gauss_legendre_quad.m``.  The MATLAB version tabulated orders
    2, 3 and 4 and fell back to the 1-point rule for anything else; here any
    order >= 1 is accepted, and anything smaller falls back to 1 point.
    """
    order = int(order)
    if order < 1:
        order = 1
    xi, wi = np.polynomial.legendre.leggauss(order)
    return xi, wi


def gauss_trgl(m: int = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Abscissas ``(xi, eta)`` and weights ``w`` for Gauss integration over a
    flat triangle in the xi-eta plane, in barycentric coordinates.

    Port of ``gauss_trgl.m`` (FSELIB, Pozrikidis, chapter 4).  The weights sum
    to one, i.e. they are normalised by the area of the reference triangle.

    Parameters
    ----------
    m : order of the quadrature, one of 1, 3, 4, 6, 7, 9, 12, 13.
        Any other value falls back to 3, as in the original.
    """
    if m not in _TRGL_ORDERS:
        print(" Gauss_trgl: Chosen number of points is not available; Will take m=3")
        m = 3

    third = 1.0 / 3.0

    if m == 1:
        xi = [third]
        eta = [third]
        w = [1.0]

    elif m == 3:
        xi = [1 / 6, 2 / 3, 1 / 6]
        eta = [1 / 6, 1 / 6, 2 / 3]
        w = [third, third, third]

    elif m == 4:
        xi = [third, 0.2, 0.6, 0.2]
        eta = [third, 0.2, 0.2, 0.6]
        w = [-27 / 48, 25 / 48, 25 / 48, 25 / 48]

    elif m == 6:
        al = 0.816847572980459
        be = 0.445948490915965
        ga = 0.108103018168070
        de = 0.091576213509771
        o1 = 0.109951743655322
        o2 = 0.223381589678011
        xi = [de, al, de, be, ga, be]
        eta = [de, de, al, be, be, ga]
        w = [o1, o1, o1, o2, o2, o2]

    elif m == 7:
        al = 0.797426958353087
        be = 0.470142064105115
        ga = 0.059715871789770
        de = 0.101286507323456
        o1 = 0.125939180544827
        o2 = 0.132394152788506
        xi = [de, al, de, be, ga, be, third]
        eta = [de, de, al, be, be, ga, third]
        w = [o1, o1, o1, o2, o2, o2, 0.225]

    elif m == 9:
        al = 0.124949503233232
        qa = 0.165409927389841
        rh = 0.797112651860071
        de = 0.437525248383384
        ru = 0.037477420750088
        o1 = 0.205950504760887
        o2 = 0.063691414286223
        xi = [de, al, de, qa, ru, rh, qa, ru, rh]
        eta = [de, de, al, ru, qa, qa, rh, rh, ru]
        w = [o1, o1, o1, o2, o2, o2, o2, o2, o2]

    elif m == 12:
        al = 0.873821971016996
        be = 0.249286745170910
        ga = 0.501426509658179
        de = 0.063089014491502
        rh = 0.636502499121399
        qa = 0.310352451033785
        ru = 0.053145049844816
        o1 = 0.050844906370207
        o2 = 0.116786275726379
        o3 = 0.082851075618374
        xi = [de, al, de, be, ga, be, qa, ru, rh, qa, ru, rh]
        eta = [de, de, al, be, be, ga, ru, qa, qa, rh, rh, ru]
        w = [o1, o1, o1, o2, o2, o2, o3, o3, o3, o3, o3, o3]

    else:  # m == 13
        al = 0.479308067841923
        be = 0.065130102902216
        ga = 0.869739794195568
        de = 0.260345966079038
        rh = 0.638444188569809
        qa = 0.312865496004875
        ru = 0.048690315425316
        o1 = 0.175615257433204
        o2 = 0.053347235608839
        o3 = 0.077113760890257
        o4 = -0.149570044467670
        xi = [de, al, de, be, ga, be, qa, ru, rh, qa, ru, rh, third]
        eta = [de, de, al, be, be, ga, ru, qa, qa, rh, rh, ru, third]
        w = [o1, o1, o1, o2, o2, o2, o3, o3, o3, o3, o3, o3, o4]

    return np.asarray(xi), np.asarray(eta), np.asarray(w)


_FAMILY = {1: "triangle", 2: "quad", 3: "quad"}


def quadrature_rule(family, order: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Quadrature points ``(xi, eta)`` and weights ``w`` for a reference element.

    ``family`` is ``"triangle"`` or ``"quad"`` (the legacy element type
    numbers 1, 2 and 3 are accepted too).

    * triangle: ``gauss_trgl(order)`` with the weights scaled by the
      reference-triangle area (1/2) so that ``sum(w) == 1/2``.
    * quad: tensor product of the 1-D Gauss-Legendre rule of the given order
      on [-1, 1]^2, ``sum(w) == 4``.
    """
    family = _FAMILY.get(family, family)
    if family == "triangle":
        xi, eta, w = gauss_trgl(order)
        return xi, eta, 0.5 * w
    if family == "quad":
        x1, w1 = gauss_legendre_quad(order)
        xi, eta = np.meshgrid(x1, x1, indexing="ij")
        w = np.outer(w1, w1)
        return xi.ravel(), eta.ravel(), w.ravel()
    raise ValueError(f"Unknown element family {family!r}; expected 'triangle' or 'quad'")
