"""Mantle to magma: magma-type names, the scarce-oxide spinel fix, and the guided example."""
import pytest

from magma.classify import magma_type, tas_name
from magma.composition import PRESETS, Composition
from magma.differentiation import iter_differentiate, magma_type_stages, melt_source, primary_with_water
from magma.equilibrium import Equilibrator


def _c(si, alk):
    return Composition({"SiO2": si, "Na2O": alk, "Al2O3": 100 - si - alk})


@pytest.mark.parametrize("si, alk, tas, kind", [
    (43, 2, "picrobasalt", "ultramafic melt"), (48, 3, "basalt", "basalt"),
    (54, 3, "basaltic andesite", "basaltic andesite"), (60, 4, "andesite", "andesite"),
    (66, 5, "dacite", "dacite"), (74, 8, "rhyolite", "rhyolite"), (48, 6, "trachybasalt", "basalt"),
    (60, 10.6, "trachyte or trachydacite", "andesite"),
])
def test_tas_and_silica_names(si, alk, tas, kind):
    assert tas_name(_c(si, alk)) == tas
    assert magma_type(_c(si, alk)) == kind


def test_presets_are_basalts():
    assert magma_type(PRESETS["MORB"]) == "basalt" and tas_name(PRESETS["MORB"]) == "basalt"


def test_primary_with_water_keeps_proportions():
    p = primary_with_water(PRESETS["MORB"], 2.0)
    assert p.get("H2O") == pytest.approx(2.0)
    assert p.total == pytest.approx(100.0)
    assert p.get("SiO2") / p.get("MgO") == pytest.approx(PRESETS["MORB"].get("SiO2") / PRESETS["MORB"].get("MgO"))


# an evolved, Cr-poor ferrobasalt liquid from a dry fractionation run (1100 C, 1 kbar, QFM)
FERROBASALT = Composition({"SiO2": 43.16, "TiO2": 1.61, "Al2O3": 13.12, "Fe2O3": 4.06, "Cr2O3": 0.037,
                           "FeO": 20.32, "MgO": 3.71, "CaO": 4.73, "Na2O": 8.36, "K2O": 0.89})


def test_spinel_rich_in_a_scarce_oxide_can_still_be_added():
    # its most stable composition is Cr-rich; 1e-3 mol of it would need more Cr2O3 than there is
    r = Equilibrator(FERROBASALT, fo2="QFM", phases=("spinel",)).equilibrate(1100, 1000)
    assert [s.name for s in r.solids] == ["spinel"]


def test_titanomagnetite_copy_is_dropped_without_chromium():
    names = [p.name for p in Equilibrator(PRESETS["MORB"]).phases]
    assert "spinel" in names and "titanomagnetite" not in names
    names = [p.name for p in Equilibrator(PRESETS["KLB-1 peridotite"]).phases]
    assert "spinel" in names and "titanomagnetite" in names


def test_mantle_to_magma():
    src = melt_source(PRESETS["KLB-1 peridotite"], 1380, 10000)
    assert 0.05 < src.result.melt_fraction < 0.20                       # about 11%
    assert magma_type(src.primary) in ("ultramafic melt", "basalt")
    assert src.primary.get("MgO") > 10                                  # a primitive, picritic melt
    steps = list(iter_differentiate(src.primary, 1340, 940, step=20, p_crust=2000, h2o=2.0, fo2="NNO"))
    first = {}
    for r in steps:
        for s in r.solids:
            first.setdefault(s.name, r.t_c)
    # Bowen's order: olivine first, then pyroxene, then plagioclase; magnetite before plagioclase
    assert min(first, key=lambda n: -first[n]) == "olivine"
    assert first["olivine"] > first["clinopyroxene"] > first["feldspar"]
    assert "titanomagnetite" in first and first["titanomagnetite"] > first["feldspar"]
    kinds = [st.magma_type for st in magma_type_stages(steps)]
    assert "basaltic andesite" in kinds and kinds[-1] in ("andesite", "dacite")
    assert steps[-1].liquid_grams < 40
