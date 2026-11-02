"""Holland & Powell end members (magma.hp.endmembers) against their own definitions and known phase equilibria."""
import math

import pytest
from scipy.optimize import brentq

from magma.hp.endmembers import DATA, gibbs, gibbs_j
from magma.thermo import water


def _s_v(name):
    h = 1e-3
    s = -(gibbs(name, 298.15 + h, 0.001) - gibbs(name, 298.15 - h, 0.001)) / (2 * h)
    v = (gibbs(name, 298.15, 0.001 + h) - gibbs(name, 298.15, 0.001)) / h
    return s, v


@pytest.mark.parametrize("name", ["fo", "fa", "di", "en", "an", "py", "alm", "mu", "phl", "tr", "parg", "q", "ru"])
def test_reference_state(name):
    """At 298.15 K and 1 bar, -dG/dT and dG/dP give the tabulated S and V."""
    s, v = _s_v(name)
    assert s == pytest.approx(DATA[name].s, abs=2e-6)
    assert v == pytest.approx(DATA[name].v, rel=1e-4)


def _boundary(reaction, t_k):
    return brentq(lambda p: sum(c * gibbs(n, t_k, p) for n, c in reaction.items()), 0.001, 80.0)


def test_albite_jadeite_quartz():
    """Albite = jadeite + quartz: 16.2 kbar at 600 C (Holland 1980, P = 0.35 + 0.0265 T)."""
    assert _boundary({"jd": 1, "q": 1, "ab": -1}, 873.15) == pytest.approx(16.24, abs=0.6)


def test_quartz_coesite():
    """Quartz = coesite near 30 kbar at 1000 C (Bose and Ganguly 1995 put it at 30.5)."""
    assert _boundary({"coe": 1, "q": -1}, 1273.15) == pytest.approx(30.5, abs=1.0)


def test_al2sio5_triple_point():
    """Kyanite-andalusite-sillimanite near 550 C and 4.5 kbar (Pattison 1992)."""
    f = lambda t: _boundary({"sill": 1, "ky": -1}, t) - _boundary({"ky": 1, "and": -1}, t)
    t = brentq(f, 750.0, 950.0)
    assert t - 273.15 == pytest.approx(550, abs=25)
    assert _boundary({"ky": 1, "and": -1}, t) == pytest.approx(4.5, abs=0.5)


def test_quartz_lambda_transition():
    """The Landau term puts the alpha-beta heat-capacity spike at 847 K at 1 bar."""
    def cp(t, h=0.01):
        return -t * (gibbs("q", t + h, 0.001) - 2 * gibbs("q", t, 0.001) + gibbs("q", t - h, 0.001)) / h ** 2
    assert cp(846.0) > 1.5 * cp(860.0)
    assert cp(860.0) == pytest.approx(cp(870.0), rel=0.02)


@pytest.mark.parametrize("t_c,p_bar", [(700, 1000), (900, 3000), (1100, 3000)])
def test_water_against_iapws(t_c, p_bar):
    """Pitzer-Sterner H2O against IAPWS-95: G(T, P) - G(T, 1 bar) within 0.25 kJ to 3 kbar."""
    t = t_c + 273.15
    hp = gibbs_j("H2O", t, p_bar) - gibbs_j("H2O", t, 1.0)
    iapws = water.gibbs_iapws(t, p_bar) - water.gibbs_iapws(t, 1.0)
    assert hp == pytest.approx(iapws, abs=250.0)


def test_liquid_end_member_is_finite():
    for name in ("qL", "abL", "foL", "h2oL"):
        assert math.isfinite(gibbs(name, 1473.15, 2.0))


def test_aqueous_species_refused():
    with pytest.raises(ValueError):
        gibbs("Na+", 1000.0, 1.0)


def _oracle_rows():
    from pathlib import Path
    rows = []
    for line in (Path(__file__).parent / "data" / "hp_em_magemin.tsv").read_text().splitlines():
        if line.startswith("#") or line.startswith("name"):
            continue
        n, t, p, g = line.split("\t")
        rows.append((n, float(t), float(p), float(g)))
    return rows


def test_against_magemin():
    """59 end members at three (T, P) points: G matches MAGEMin to 1 J (print precision is 0.001 J)."""
    worst = max(abs(gibbs(n, t, p) - g) for n, t, p, g in _oracle_rows())
    assert worst < 1e-3


def test_gases_without_eos_refused():
    with pytest.raises(ValueError):
        gibbs("CH4", 1000.0, 1.0)
