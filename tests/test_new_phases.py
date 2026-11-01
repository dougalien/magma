"""Pigeonite, garnet and the water fluid."""
import math

import pytest

from magma.composition import PRESETS, Composition
from magma.equilibrium import Equilibrator
from magma.melting import iter_melt
from magma.path import crystallize
from magma.thermo.phases import g_ex_garnet
from magma.thermo.water_fluid import WATER_FLUID


def test_morb_grows_augite_and_pigeonite_together():
    r = crystallize(PRESETS["MORB"], 1160, 1140, 10)[-1]
    by = {s.name: s for s in r.solids}
    assert {"clinopyroxene", "pigeonite", "feldspar"} <= set(by)
    ca = lambda s: 1 - s.x["jadeite"] - s.x["clinoenstatite"]      # Ca per formula
    assert ca(by["clinopyroxene"]) > 0.6 and ca(by["pigeonite"]) < 0.35


def test_garnet_matches_melts_definition():
    # GarnetBerman.m's H - T S + P V macros (with ideal mixing on 3 sites) at this point
    t, p, x = 1600.0, 20000.0, (0.15, 0.12, 0.73)
    g = g_ex_garnet(x, t, p) + 3 * 8.3143 * t * sum(v * math.log(v) for v in x)
    assert g == pytest.approx(-27637.50201615857, abs=1e-6)


def test_garnet_lherzolite_at_3_gpa():
    m = next(iter_melt(PRESETS["KLB-1 peridotite"], 1525, 1525, 25, p_bar=30000, mode="batch"))
    names = {s.name for s in m.result.solids}
    assert {"garnet", "olivine", "orthopyroxene", "clinopyroxene"} <= names and "pigeonite" not in names
    gt = next(s for s in m.result.solids if s.name == "garnet")
    assert 0.6 < gt.x["pyrope"] < 0.85


def _water_affinity(h2o, p, t_c=1200.0):
    wt = dict(PRESETS["MORB"].wt, H2O=h2o)
    e = Equilibrator(Composition(wt))
    liq = e._liquid({}, t_c, p)
    return WATER_FLUID.affinity(e._mu_liquid(liq, t_c + 273.15, p), t_c + 273.15, p)[0]


def test_water_solubility_in_morb():
    # about 3 wt% at 1 kbar and 1200 C (MELTS-like; Dixon et al. 1995 measure 2-3 wt%)
    assert _water_affinity(2.0, 1000.0) < 0 < _water_affinity(3.5, 1000.0)
    assert _water_affinity(3.5, 2000.0) < 0                         # more soluble at higher P


def test_hydrous_morb_saturates_on_cooling():
    wt = dict(PRESETS["MORB"].wt, H2O=2.5)
    path = crystallize(Composition(wt), 1100, 1020, 20, p_bar=1000)
    assert not any(s.name == "water" for s in path[0].solids)
    assert any(s.name == "water" for s in path[-1].solids)
    assert not any(s.name == "feldspar" for s in path[0].solids)   # water delays plagioclase
