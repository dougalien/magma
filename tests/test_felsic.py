"""The late minerals (quartz, tridymite, K-feldspar, muscovite, biotite, hornblende) against
ThermoEngine's MELTS v1.0.2 equilibrium example (Early Bishop Tuff, 760 C, 175 MPa;
Documentation/source/MELTS-v1.0.2-equilibrium.rst)."""
import pytest

from magma.composition import CARDS, PRESETS
from magma.equilibrium import Equilibrator
from magma.thermo.felsic import quartz_gibbs
from magma.thermo.water_fluid import water_gibbs_1bar

BT = PRESETS["Bishop Tuff rhyolite"]


@pytest.fixture(scope="module")
def bt760():
    return Equilibrator(BT, fo2=None).equilibrate(760, 1750)


@pytest.mark.parametrize("phase, melts", [
    ("biotite", -19819.30), ("muscovite", -25188.37), ("hornblende", -54565.30), ("tridymite", -2141.04),
])
def test_undersaturated_phases_match_melts(bt760, phase, melts):
    # MELTS reports A > 0 for undersaturated; Magma the opposite sign
    assert bt760.affinities[phase] == pytest.approx(melts, abs=1500)


def test_quartz_and_sanidine_at_saturation(bt760):
    # MELTS: 1.2 g quartz and 2.3 g sanidine (Or49) out of 105.6 g; Magma is within 0.3 kJ of both
    aff = bt760.affinities
    for ph in ("quartz", "alkali-feldspar"):
        if ph in aff:
            assert abs(aff[ph]) < 300
    assert bt760.liquid.normalized().get("SiO2") == pytest.approx(73.63, abs=0.3)


def test_granite_minimum_below_760():
    r = Equilibrator(BT, fo2=None).equilibrate(740, 1750)
    names = {s.name for s in r.solids}
    assert {"quartz", "alkali-feldspar"} <= names


def test_quartz_is_continuous_at_the_alpha_beta_transition():
    assert quartz_gibbs(848.0 - 1e-6, 1.0) == pytest.approx(quartz_gibbs(848.0 + 1e-6, 1.0), abs=1.0)   # Berman: 0.3 J


def test_water_on_the_minerals_scale():
    assert water_gibbs_1bar(298.15) == pytest.approx(-285830.0 - 298.15 * 69.95, abs=1.0)


def test_cards_name_presets():
    for c in CARDS:
        assert c["preset"] in PRESETS and c["kind"] in ("melt", "rock")
        assert (c["t_start"] > c["t_end"]) == (c["kind"] == "melt")
