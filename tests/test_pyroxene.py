"""Clinopyroxene (Sack and Ghiorso 1994c as in MELTS)."""
import numpy as np
import pytest

from magma.composition import PRESETS
from magma.equilibrium import Equilibrator
from magma.thermo import pyroxene as px
from magma.thermo.phases import CLINOPYROXENE

from test_equilibrium import MELTS_LIQ, _mu

T, P = 1493.15, 1000.0


def test_reference_states_match_melts_definitions():
    # values from evaluating ThermoEngine CpxBerman.m's EN_G, ES_G, JD_G at 1220 C, 1000 bar
    g = px.end_member_g(T, P)
    assert g[1] == pytest.approx(1612.8756645, abs=1e-4)      # clinoenstatite
    assert g[5] == pytest.approx(-4902.7959586, abs=1e-4)     # essenite, with its own ordering
    assert g[[0, 2, 3, 4, 6]] == pytest.approx(0.0, abs=1e-6)


def test_solution_gibbs_matches_melts_definition():
    # CpxBerman.m's G macro at these r, s (1220 C, 1000 bar)
    r, s = np.array([0.12, 0.05, -0.02, 0.08, 0.03, 0.15]), np.array([0.01, 0.05])
    assert px.g_rs(r, s, T, P) == pytest.approx(-9987.2881567, abs=1e-4)


def test_chemical_potentials_are_partial_derivatives():
    t, p = 1450.0, 1.0
    x = np.array([0.48, 0.23, 0.17, 0.07, -0.045, 0.087, 0.012])
    x /= x.sum()
    mu = CLINOPYROXENE.mu(x, t, p)
    g = lambda n: n.sum() * CLINOPYROXENE.gibbs(n / n.sum(), t, p)
    h = 1e-6
    num = [(g(x + h * e) - g(x - h * e)) / (2 * h) for e in np.eye(7)]
    assert mu == pytest.approx(num, abs=0.01)
    assert x @ mu == pytest.approx(CLINOPYROXENE.gibbs(x, t, p), abs=1e-6)


def test_augite_affinity_against_melts_example():
    a, x = CLINOPYROXENE.affinity(_mu(MELTS_LIQ, T, P), T, P)
    assert a == pytest.approx(-1333.92, abs=150)               # MELTS: augite 1333.92 J undersaturated
    assert x[0] > 0.3 and x[4] < 0                             # diopside-rich, negative buffonite (CaTs)


def test_morb_crystallizes_augite():
    r = Equilibrator(PRESETS["MORB"], fo2="QFM").equilibrate(1165.0, 1.0)
    names = {s.name for s in r.solids}
    assert names == {"feldspar", "olivine", "clinopyroxene"}
    cpx = next(s for s in r.solids if s.name == "clinopyroxene")
    assert "Wo" in CLINOPYROXENE.formula(np.array(list(cpx.x.values())))
    assert 3 < cpx.grams < 20
