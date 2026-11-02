"""Holland & Powell solid solutions (magma.hp.solutions) against MAGEMin and their own definitions."""
import json
from pathlib import Path

import numpy as np
import pytest

from magma.hp.endmembers import OXIDES
from magma.hp.solutions import SOLUTIONS

IG = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"]
POINTS = json.loads((Path(__file__).parent / "data" / "hp_magemin_points.json").read_text())["points"]
ALIASES = {"qndm": "usp"}     # MAGEMin's name for the spinel end member made from ulvospinel
CASES = [(pt, ph) for pt in POINTS for ph in pt["phases"] if ph["name"] in SOLUTIONS and not ph.get("pure")]


@pytest.mark.parametrize("pt,ph", CASES, ids=[f'{ph["name"]}-{pt["P"]}kbar-{pt["T_C"]}C' for pt, ph in CASES])
def test_mu_at_magemin_equilibrium(pt, ph):
    """At a MAGEMin equilibrium, each end member's mu equals sum(oxide moles x oxide chemical potential)."""
    s = SOLUTIONS[ph["name"]]
    assert tuple(ALIASES.get(e, e) for e in ph["em"]) == s.names
    t, p = pt["T_C"] + 273.15, pt["P"]
    prop = np.array(ph["frac"])
    gamma = np.array([pt["gamma"][o] for o in IG])
    comp = s.compositions()[:, [OXIDES.index(o) for o in IG]]
    mu = s.mu(prop, t, p)
    ok = np.isfinite(mu)          # an end member with an empty site (no Cr) has mu = -inf
    # end members held at MAGEMin's lower bound (p <= 1e-5, as in the fluid) are not at equilibrium;
    # they are checked against MAGEMin's own mu below instead
    eq = ok & (np.abs(prop) > 1e-5)
    assert np.max(np.abs(mu[eq] - (comp @ gamma)[eq])) < 0.01
    em_mu = np.array(ph["mu"])
    factor = (em_mu[eq] @ mu[eq]) / (mu[eq] @ mu[eq])
    assert np.max(np.abs(mu[ok] - em_mu[ok] / factor)) < 0.01


@pytest.mark.parametrize("name", sorted(SOLUTIONS))
def test_mu_is_derivative_of_g(name):
    """mu_i = d(N G)/dn_i, so the activities and excess terms are mutually consistent.

    Feldspar and the fluid are the exceptions for the ideal part: feldspar's Al-avoidance activities
    (Holland and Powell 2003) and the fluid's p_i (1 - p_H2O) activities are not the exact derivative of
    one entropy, so only their excess terms are held to this.
    """
    s = SOLUTIONS[name]
    pt, ph = next((pt, ph) for pt, ph in CASES if ph["name"] == name) if any(ph["name"] == name for _, ph in CASES) \
        else (None, None)
    if ph is None:
        pytest.skip("no MAGEMin point with this phase")
    t, p = pt["T_C"] + 273.15, pt["P"]
    n0 = np.array(ph["frac"])
    mu = s.mu(n0, t, p)

    if name in ("fsp", "fl"):
        mu = s.mu_excess(n0, t, p)

        def total(n):
            return n.sum() * s.g_excess(n / n.sum(), t, p)
    else:
        def total(n):
            return n.sum() * s.gibbs(n / n.sum(), t, p)

    h = 1e-6
    for i in range(s.n):
        if not np.isfinite(mu[i]):
            continue
        up, dn = n0.copy(), n0.copy()
        up[i] += h
        dn[i] -= h
        assert (total(up) - total(dn)) / (2 * h) == pytest.approx(mu[i], abs=1e-3)
