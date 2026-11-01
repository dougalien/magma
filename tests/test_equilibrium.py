import numpy as np
import pytest

from magma.composition import PRESETS, Composition
from magma.equilibrium import Equilibrator
from magma.path import crystallize
from magma.thermo.liquid import COMPONENTS, mu_mixing, to_components
from magma.thermo.minerals import MINERALS
from magma.thermo.phases import FELDSPAR, OLIVINE, mu_in_liquid
from magma.thermo.standard_state import liquid_gibbs

# ThermoEngine's MELTS v1.0.2 fractionation example (Documentation/source/MELTS-v1.0.2-fractionation.rst)
MELTS_BULK = Composition({"SiO2": 48.68, "TiO2": 1.01, "Al2O3": 17.64, "Fe2O3": 0.89, "Cr2O3": 0.0425,
                          "FeO": 7.59, "MgO": 9.10, "CaO": 12.45, "Na2O": 2.65, "K2O": 0.03,
                          "P2O5": 0.08, "H2O": 0.2})
MELTS_LIQ = Composition({"SiO2": 48.53, "TiO2": 1.02, "Al2O3": 17.33, "Fe2O3": 0.90, "Cr2O3": 0.04,
                         "FeO": 7.67, "MgO": 9.20, "CaO": 12.34, "Na2O": 2.65, "K2O": 0.03,
                         "P2O5": 0.08, "H2O": 0.2})


def _mu(liq, t, p):
    return np.array([liquid_gibbs(c, t, p) for c in COMPONENTS]) + mu_mixing(to_components(liq).x, t)


def test_affinities_against_melts_example():
    t, p = 1493.15, 1000.0
    mu = _mu(MELTS_LIQ, t, p)
    fo = MINERALS["forsterite"].gibbs(t, p) - mu_in_liquid("forsterite", mu)
    assert fo == pytest.approx(3494.95, abs=25)          # MELTS: forsterite 3494.95 J undersaturated
    a_pl, x_pl = FELDSPAR.affinity(mu, t, p)             # MELTS: plagioclase present, An81
    assert abs(a_pl) < 300 and x_pl[1] == pytest.approx(0.81, abs=0.02)
    a_ol, _ = OLIVINE.affinity(mu, t, p)                  # MELTS: olivine 118.94 J undersaturated
    assert -500 < a_ol < 0


def test_end_member_potentials_consistent_with_gibbs():
    t, p = 1450.0, 1.0
    x = np.array([0.3, 0.68, 0.02])
    mu = FELDSPAR.mu(x, t, p)
    assert float(x @ mu) == pytest.approx(FELDSPAR.gibbs(x, t, p), rel=1e-9)


def test_melts_example_equilibrium():
    r = Equilibrator(MELTS_BULK, fo2=None).equilibrate(1220.0, 1000.0)
    by = {s.name: s for s in r.solids}
    assert set(by) == {"feldspar", "spinel"}              # MELTS: plagioclase + 0.02 g spinel
    assert by["feldspar"].x["anorthite"] == pytest.approx(0.81, abs=0.02)
    assert 0.0 < by["feldspar"].grams < 3.0               # MELTS: 1.47 g
    assert by["spinel"].grams < 0.2


def test_morb_cooling_mass_balance_and_order():
    path = crystallize(PRESETS["MORB"], 1200, 1150, 10)
    assert not path[0].solids                             # above the liquidus
    for r in path:
        assert r.liquid_grams + sum(s.grams for s in r.solids) == pytest.approx(100.0, rel=1e-9)
    assert path[-1].melt_fraction < 0.9
