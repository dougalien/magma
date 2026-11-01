"""Spinel against the MELTS model definitions (SpinelBerman.m, evaluated independently) and the
MELTS example liquid, which MELTS reports as spinel-saturated at 1220 C and 1 kbar."""
import numpy as np
import pytest

from magma.composition import PRESETS
from magma.path import crystallize
from magma.thermo import spinel as sp

from test_equilibrium import MELTS_LIQ, _mu

T, P = 1473.15, 1000.0            # where the model definitions were evaluated
T_EX = 1493.15                     # the MELTS example, 1220 C


def test_reference_states_match_melts_definitions():
    ref = [-125520.0, -40264.088742032895, -60072.32961969408, -27411.372644696956, -104006.82590548985]
    assert sp.end_member_g(T, P) == pytest.approx(ref, abs=1e-3)


def test_solution_gibbs_matches_melts_definition():
    z = [0.3, 0.2, 0.1, 0.15, 0.1, 0.05, 0.03]
    assert sp.g_and_grad(z, T, P)[0] == pytest.approx(-74313.2742353875, abs=1e-3)


def test_chemical_potentials_are_partial_derivatives():
    x = np.array([0.15, 0.25, 0.1, 0.4, 0.1])
    mu = sp.SPINEL.mu(x, T, P)
    for j in range(5):
        h = 1e-6
        n = x.copy()
        n[j] += h
        g1 = n.sum() * sp.SPINEL.gibbs(n / n.sum(), T, P)
        n[j] -= 2 * h
        g0 = n.sum() * sp.SPINEL.gibbs(n / n.sum(), T, P)
        assert (g1 - g0) / (2 * h) == pytest.approx(mu[j], abs=2.0)


def test_spinel_saturation_against_melts_example():
    a, x = sp.SPINEL.affinity(_mu(MELTS_LIQ, T_EX, P), T_EX, P)
    assert abs(a) < 500                                   # MELTS: saturated (ours +160 J)
    cr, hc, mt, spl, usp = x
    assert 2 * cr == pytest.approx(0.45, abs=0.05)         # MELTS: Cr 0.45, Mg 0.79 per 4 O
    assert spl == pytest.approx(0.79, abs=0.05)


def test_morb_reaches_titanomagnetite_without_chromium():
    steps = crystallize(PRESETS["MORB"], 1130, 1100, 10)
    last = steps[-1]
    names = {s.name for s in last.solids}
    assert "spinel" in names and "feldspar" in names and "clinopyroxene" in names
    spl = next(s for s in last.solids if s.name == "spinel")
    assert spl.x["ulvospinel"] > 0.5 and spl.x["chromite"] == 0.0
