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
    # spinel lherzolite: olivine, orthopyroxene, a Ca-pyroxene and spinel.  With the Davis et al. (2009)
    # KLB-1 the Ca-pyroxene (Ca about 0.6 per formula) is found by the low-Ca copy of the clinopyroxene
    # model, which the app labels by composition (server._solid_json)
    first = _names(steps[0])
    assert {"olivine", "orthopyroxene", "spinel"} <= first and first & {"clinopyroxene", "pigeonite"}
    # one pyroxene is used up before about 20% melt, leaving olivine, one pyroxene and spinel
    last = steps[-1]
    pyx = {"orthopyroxene", "clinopyroxene", "pigeonite"}
    assert len(_names(last) & pyx) == 1 and 0.15 < last.result.melt_fraction < 0.35


def test_fractional_melting_pools_melt_and_leaves_harzburgite():
    steps = list(iter_melt(ROCK, 1300, 1450, 25, p_bar=10000, mode="fractional"))
    pooled = [m.pooled_grams for m in steps]
    assert all(b >= a for a, b in zip(pooled, pooled[1:])) and pooled[-1] > 5
    for m in steps:                                            # residue keeps 1% melt
        r = m.result
        assert m.extracted_grams == 0 or m.extracted_grams < r.liquid_grams + 1e-9
    assert "clinopyroxene" not in _names(steps[-1]) and "orthopyroxene" in _names(steps[-1])


def test_ca_pyroxene_is_labelled_by_composition():
    from magma.app.server import _solid_json
    from magma.equilibrium import PhaseResult
    from magma.thermo.pyroxene import CLINOPYROXENE
    names = CLINOPYROXENE.end_members
    augite_like = dict(zip(names, [0.378, 0.328, 0.12, 0.096, -0.08, 0.121, 0.037]))
    assert _solid_json(PhaseResult("pigeonite", 10.0, {}, augite_like))["label"] == "Augite"
    pigeonite_like = dict(zip(names, [0.06, 0.74, 0.14, 0.02, 0.0, 0.03, 0.01]))
    assert _solid_json(PhaseResult("pigeonite", 10.0, {}, pigeonite_like))["label"] == "Pigeonite"
