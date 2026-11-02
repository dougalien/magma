"""The Holland & Powell engine adapter (magma.hp.engine) and the solver's oxygen buffer and warm starts."""
import json
from pathlib import Path

import numpy as np
import pytest

from magma.composition import PRESETS, Composition
from magma.hp import engine
from magma.hp.endmembers import R, gibbs
from magma.hp.minimize import COMPONENTS, equilibrate

POINTS = json.loads((Path(__file__).parent / "data" / "hp_magemin_points.json").read_text())["points"]
WET_BASALT = next(pt for pt in POINTS if pt["P"] == 3.0 and pt["T_C"] == 950 and pt["bulk"][10] > 0)


def test_components_round_trip():
    c = Composition({"SiO2": 50.0, "FeO": 8.0, "Fe2O3": 2.0, "MgO": 8.0, "CaO": 11.0, "Al2O3": 15.0,
                     "Na2O": 2.5, "K2O": 0.5, "TiO2": 1.5, "H2O": 1.5, "MnO": 0.2})
    m = engine.to_components(c)
    assert m[COMPONENTS.index("O")] == pytest.approx(2.0 / 159.69 * 100 / c.total, rel=1e-3)
    back = engine.to_composition(m)
    for ox in ("SiO2", "FeO", "Fe2O3", "MgO", "H2O"):
        assert back.get(ox) == pytest.approx(c.normalized().get(ox) * 100 / (100 - 0.2 * 100 / c.total), rel=1e-3)
    assert "MnO" not in back.wt


def test_buffer_reproduces_closed_system():
    """At the fO2 that a closed system reaches by itself, the buffered system keeps the same Fe3+."""
    t = WET_BASALT["T_C"] + 273.15
    closed = equilibrate(WET_BASALT["bulk"], t, 3.0)
    lf = (2 * closed.gamma[COMPONENTS.index("O")] - gibbs("O2", t, 0.001)) / (R * t * np.log(10))
    buffered = equilibrate(WET_BASALT["bulk"], t, 3.0, log_fo2=lf)
    assert buffered.bulk[COMPONENTS.index("O")] == pytest.approx(WET_BASALT["bulk"][COMPONENTS.index("O")], abs=1e-6)
    assert np.max(np.abs(buffered.gamma - closed.gamma)) < 1e-4
    more = equilibrate(WET_BASALT["bulk"], t, 3.0, log_fo2=lf + 1)
    assert more.bulk[COMPONENTS.index("O")] > buffered.bulk[COMPONENTS.index("O")] + 0.05   # more Fe3+


def test_warm_start_gives_the_same_answer():
    t = WET_BASALT["T_C"] + 273.15
    first = equilibrate(WET_BASALT["bulk"], t, 3.0)
    warm = equilibrate(WET_BASALT["bulk"], t - 10, 3.0, previous=first)
    cold = equilibrate(WET_BASALT["bulk"], t - 10, 3.0)
    assert sorted(p.name for p in warm.phases) == sorted(p.name for p in cold.phases)
    assert warm.g == pytest.approx(cold.g, abs=1e-4)


def test_cooling_path_and_names():
    """Wet andesite at 2 kbar: a few fractional steps below the liquidus give plagioclase and pyroxene,
    a silica-richer melt, and app-ready phases."""
    steps = list(engine.iter_crystallize(PRESETS["Andesite (bulk crust)"], 1000, 950, 25, 2000, "NNO", "fractional"))
    assert len(steps) == 3
    first, last = steps[0], steps[-1]
    names = {s.name for r in steps for s in r.solids}
    assert {"feldspar", "orthopyroxene"} <= names
    assert last.liquid.normalized(volatile_free=True).get("SiO2") > first.liquid.normalized(volatile_free=True).get("SiO2")
    plag = next(s for s in first.solids if s.name == "feldspar")
    assert 30 < plag.numbers["An"] < 90 and plag.formula.endswith("O8")
    assert first.liquid_grams + sum(s.grams for s in first.solids) == pytest.approx(100.0, abs=0.5)


def test_server_routes_engine():
    from magma.app.server import parse_inputs, step_json
    bulk, kw = parse_inputs({"bulk": PRESETS["Andesite (bulk crust)"].wt, "engine": "hp", "t_start": 1000,
                             "t_end": 1000, "step": 5, "p_bar": 2000, "fo2": "NNO"})
    assert kw["engine"] == "hp"
    r = next(engine.iter_crystallize(bulk, 1000, 1000, 5, 2000, "NNO", "equilibrium"))
    j = step_json(r)
    assert j["solids"] and all(s["formula"] for s in j["solids"] if s["name"] != "water")
    with pytest.raises(ValueError):
        parse_inputs({"bulk": {"SiO2": 50}, "engine": "other"})
