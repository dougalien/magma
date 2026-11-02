"""The Holland & Powell equilibrium solver (magma.hp.minimize) against MAGEMin equilibria."""
import json
from pathlib import Path

import numpy as np
import pytest

from magma.hp.minimize import COMPONENTS, equilibrate

POINTS = json.loads((Path(__file__).parent / "data" / "hp_magemin_points.json").read_text())["points"]


def _find(p, t_c):
    return next(pt for pt in POINTS if pt["P"] == p and pt["T_C"] == t_c and pt["bulk"][10] > 0)


@pytest.mark.parametrize("p,t_c", [(3.0, 950), (5.0, 800), (4.0, 800)],
                         ids=["wet-basalt-amphibole", "andesite-8-phases", "wet-andesite-amp-bi"])
def test_assemblage_and_potentials_match_magemin(p, t_c):
    pt = _find(p, t_c)
    r = equilibrate(pt["bulk"], t_c + 273.15, p)
    assert r.converged
    assert sorted(ph.name for ph in r.phases) == sorted(ph["name"] for ph in pt["phases"])
    ref = np.array([pt["gamma"][c] for c in COMPONENTS])
    present = np.array(pt["bulk"]) > 0
    assert np.max(np.abs(r.gamma - ref)[present]) < 0.01


def test_mass_balance():
    pt = _find(5.0, 800)
    r = equilibrate(pt["bulk"], 1073.15, 5.0)
    total = sum(ph.oxide_moles for ph in r.phases)
    assert np.allclose(total, pt["bulk"], atol=1e-8)
