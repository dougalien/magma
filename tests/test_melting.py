"""Partial melting of KLB-1 peridotite at 1 GPa (spinel field)."""
import pytest

from magma.composition import PRESETS
from magma.melting import iter_melt

ROCK = PRESETS["KLB-1 peridotite"]


def _names(m):
    return {s.name for s in m.result.solids if s.grams > 1e-3}


def test_batch_melting_order_and_mass_balance():
    steps = list(iter_melt(ROCK, 1325, 1425, 25, p_bar=10000, mode="batch"))
    assert [m.result.t_c for m in steps] == sorted(m.result.t_c for m in steps)
    f = [m.result.melt_fraction for m in steps]
    assert all(b >= a - 1e-6 for a, b in zip(f, f[1:]))          # melt grows on heating
    for m in steps:
        r = m.result
        assert r.liquid_grams + sum(s.grams for s in r.solids) == pytest.approx(100.0, rel=1e-9)
    assert {"olivine", "orthopyroxene", "clinopyroxene", "spinel"} <= _names(steps[0])
    # cpx is used up near 20% melt (Hirose and Kushiro 1993 find it near 1350-1400 C at 1 GPa)
    last = steps[-1]
    assert "clinopyroxene" not in _names(last) and 0.15 < last.result.melt_fraction < 0.35


def test_fractional_melting_pools_melt_and_leaves_harzburgite():
    steps = list(iter_melt(ROCK, 1300, 1450, 25, p_bar=10000, mode="fractional"))
    pooled = [m.pooled_grams for m in steps]
    assert all(b >= a for a, b in zip(pooled, pooled[1:])) and pooled[-1] > 5
    for m in steps:                                            # residue keeps 1% melt
        r = m.result
        assert m.extracted_grams == 0 or m.extracted_grams < r.liquid_grams + 1e-9
    assert "clinopyroxene" not in _names(steps[-1]) and "orthopyroxene" in _names(steps[-1])
