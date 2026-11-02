"""Naming a melt: magma type by silica, the total-alkali-silica (TAS) name, and Bowen's series.

TAS: the IUGS classification of volcanic rocks, Le Bas, Le Maitre, Streckeisen and Zanettin (1986)
J. Petrol. 27:745-750 (Fig. 1 there; field corners as given by Le Maitre et al. 2002), on a
volatile-free basis with all iron counted in.  It names a melt as the rock it would quench to.

Magma type by silica: the usual teaching scheme behind "basalt, andesite, dacite, rhyolite" (the
SiO2 limits of the subalkaline TAS fields: 45, 52, 57, 63 and, for rhyolite, 69 wt%; the TAS
dacite-rhyolite boundary slopes, 69 is where it meets 8 wt% alkalis).
"""
from __future__ import annotations

from .composition import Composition

# corners (SiO2, Na2O + K2O) of each TAS field, wt%
TAS_FIELDS: dict[str, list[tuple[float, float]]] = {
    "picrobasalt": [(41, 0), (41, 3), (45, 3), (45, 0)],
    "basalt": [(45, 0), (45, 5), (52, 5), (52, 0)],
    "basaltic andesite": [(52, 0), (52, 5), (57, 5.9), (57, 0)],
    "andesite": [(57, 0), (57, 5.9), (63, 7), (63, 0)],
    "dacite": [(63, 0), (63, 7), (69, 8), (77, 0)],
    "rhyolite": [(77, 0), (69, 8), (69, 13), (77, 13), (87.5, 4.7), (85.9, 0)],
    "trachybasalt": [(45, 5), (49.4, 7.3), (52, 5)],
    "basaltic trachyandesite": [(52, 5), (49.4, 7.3), (53, 9.3), (57, 5.9)],
    "trachyandesite": [(57, 5.9), (53, 9.3), (57.6, 11.7), (63, 7)],
    "trachyte or trachydacite": [(63, 7), (57.6, 11.7), (61, 13.5), (69, 13), (69, 8)],
    "tephrite or basanite": [(41, 3), (41, 7), (45, 9.4), (49.4, 7.3), (45, 5), (45, 3)],
    "phonotephrite": [(49.4, 7.3), (45, 9.4), (48.4, 11.5), (53, 9.3)],
    "tephriphonolite": [(53, 9.3), (48.4, 11.5), (52.5, 14), (57.6, 11.7)],
    "phonolite": [(57.6, 11.7), (52.5, 14), (48.2, 16), (55, 18), (61, 13.5)],
    "foidite": [(41, 3), (37, 3), (35, 9), (37, 14), (48.2, 16), (52.5, 14), (48.4, 11.5), (45, 9.4), (41, 7)],
}

# boundaries drawn on the TAS plot (each a polyline)
TAS_LINES: list[list[tuple[float, float]]] = [
    [(41, 0), (41, 7), (45, 9.4), (48.4, 11.5), (52.5, 14), (48.2, 16)],
    [(45, 0), (45, 5), (49.4, 7.3), (53, 9.3), (57.6, 11.7), (61, 13.5)],
    [(41, 3), (45, 3)],
    [(45, 5), (52, 5), (57, 5.9), (63, 7), (69, 8), (77, 0)],
    [(52, 0), (52, 5), (49.4, 7.3), (45, 9.4)],
    [(57, 0), (57, 5.9), (53, 9.3), (48.4, 11.5)],
    [(63, 0), (63, 7), (57.6, 11.7), (52.5, 14)],
    [(69, 8), (69, 13)],
]

SILICA_TYPES = [(45.0, "ultramafic melt"), (52.0, "basalt"), (57.0, "basaltic andesite"),
                (63.0, "andesite"), (69.0, "dacite"), (101.0, "rhyolite")]

# Bowen's reaction series (Bowen 1928, The Evolution of the Igneous Rocks): the order in which minerals
# crystallize from a cooling basaltic melt.  model: the Magma phases that stand for the mineral.
BOWEN_DISCONTINUOUS = [
    ("olivine", ("olivine",)),
    ("pyroxene", ("orthopyroxene", "pigeonite", "clinopyroxene")),
    ("amphibole", ()),
    ("biotite", ()),
]
BOWEN_CONTINUOUS = [("Ca-rich plagioclase", ("feldspar",)), ("Na-rich plagioclase", ("feldspar",))]
BOWEN_LATE = [("K-feldspar", ()), ("muscovite", ()), ("quartz", ())]


def _tas_point(comp: Composition) -> tuple[float, float]:
    c = comp.with_fe3_ratio(0.0).normalized(volatile_free=True)
    return c.get("SiO2"), c.get("Na2O") + c.get("K2O")


def tas_point(comp: Composition) -> tuple[float, float]:
    """(SiO2, Na2O + K2O) in wt%, volatile-free."""
    return _tas_point(comp)


def _inside(pt, poly) -> bool:
    x, y = pt
    inside = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def tas_name(comp: Composition) -> str:
    """The TAS field of a melt (Le Bas et al. 1986)."""
    pt = _tas_point(comp)
    for name, poly in TAS_FIELDS.items():
        if _inside(pt, poly):
            return name
    s, a = pt
    if s < 41:
        return "ultrabasic (below the TAS fields)" if a < 3 else "foidite"
    return "off the TAS diagram"


def magma_type(comp: Composition) -> str:
    """Basalt, basaltic andesite, andesite, dacite or rhyolite, by SiO2 alone (volatile-free)."""
    s, _ = _tas_point(comp)
    for limit, name in SILICA_TYPES:
        if s < limit:
            return name
    return "rhyolite"
