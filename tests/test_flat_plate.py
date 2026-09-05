"""Laminar flat plate at Re_L = 1e5 against the Blasius solution."""

import pathlib
import sys
import time

import numpy as np
import pytest

pytest.importorskip("meshio")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "examples"))
import flat_plate  # noqa: E402


def test_blasius_solution():
    fpp0, profile = flat_plate.blasius()
    assert abs(fpp0 - 0.33206) < 1e-4
    assert abs(profile(5.0) - 0.99155) < 1e-3          # f'(5) of the Blasius profile


def test_skin_friction_and_profiles_match_blasius():
    t0 = time.time()
    sol, res = flat_plate.run(nu=1e-5, verbose=False)
    assert sol.info["converged"]
    assert time.time() - t0 < 120
    rex, cf, cf_blasius = res["Re_x"], res["Cf"], res["Cf_blasius"]
    sel = (rex > 1e4) & (rex < 1e5)
    assert np.abs(cf[sel] / cf_blasius[sel] - 1).max() < 0.06
    assert np.abs(cf[sel] / cf_blasius[sel] - 1).mean() < 0.04
    for x, err in res["profile_errors"].items():
        assert err < 0.02, (x, err)
    # near the leading edge the Navier-Stokes skin friction exceeds Blasius
    lead = (rex > 1e3) & (rex < 3e3)
    assert (cf[lead] > cf_blasius[lead]).all()
