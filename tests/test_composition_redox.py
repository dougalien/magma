import math

import pytest

from magma.composition import MOLAR_MASS, PRESETS, Composition
from magma.redox import at_fo2, buffer_log_fo2, fe3_ratio, log_fo2, melt_log_fo2

MORB = PRESETS["MORB"]


def test_molar_masses():
    assert MOLAR_MASS["SiO2"] == pytest.approx(60.0843, abs=1e-3)
    assert MOLAR_MASS["FeO"] == pytest.approx(71.8444, abs=1e-3)
    assert MOLAR_MASS["Fe2O3"] == pytest.approx(159.6882, abs=1e-3)


def test_normalize_and_fractions():
    c = MORB.normalized()
    assert c.total == pytest.approx(100.0)
    x = c.mole_fractions()
    assert sum(x.values()) == pytest.approx(1.0)
    assert sum(c.cation_fractions().values()) == pytest.approx(1.0)
    wet = Composition({**MORB.wt, "H2O": 2.0})
    assert wet.normalized(volatile_free=True).get("H2O") == 0.0


def test_iron_split_keeps_total_iron():
    c = MORB.with_fe3_ratio(0.15)
    assert c.feo_total == pytest.approx(MORB.feo_total)
    assert c.fe3_ratio == pytest.approx(0.15)
    assert c.mg_number > MORB.mg_number          # less Fe2+, higher Mg#


def test_unknown_oxide_rejected():
    with pytest.raises(ValueError):
        Composition({"SiO2": 50, "FeO3": 2})


def test_buffers_at_1200C():
    # Frost (1991): QFM about -8.3 at 1200 C, NNO about 0.7 log units above it, IW about 3.6 below
    assert buffer_log_fo2("QFM", 1200) == pytest.approx(-8.30, abs=0.02)
    assert buffer_log_fo2("NNO", 1200) - buffer_log_fo2("QFM", 1200) == pytest.approx(0.73, abs=0.05)
    assert buffer_log_fo2("QFM", 1200) - buffer_log_fo2("IW", 1200) == pytest.approx(3.6, abs=0.2)
    assert log_fo2("FMQ+1", 1200) == pytest.approx(buffer_log_fo2("QFM", 1200) + 1)
    assert log_fo2("NNO-0.5", 1100, 2000) == pytest.approx(buffer_log_fo2("NNO", 1100, 2000) - 0.5)


def test_kc91_behaviour():
    lo = fe3_ratio(MORB, 1200, 1, "QFM-2")
    mid = fe3_ratio(MORB, 1200, 1, "QFM")
    hi = fe3_ratio(MORB, 1200, 1, "QFM+2")
    assert lo < mid < hi
    # MORB glasses near QFM carry roughly 0.10 to 0.16 Fe3+/total Fe
    assert 0.05 < mid < 0.20
    # slope: a = 0.196 per ln fO2, so X_Fe2O3/X_FeO rises 10**(0.196*2)=2.47x over 2 log units
    r = lambda f: f / (2 * (1 - f))
    assert r(hi) / r(mid) == pytest.approx(10 ** (0.196 * 2), rel=1e-6)


def test_kc91_round_trip():
    target = buffer_log_fo2("QFM", 1150, 1000) + 0.7
    c = at_fo2(MORB, 1150, 1000, target)
    assert melt_log_fo2(c, 1150, 1000) == pytest.approx(target, abs=1e-9)
