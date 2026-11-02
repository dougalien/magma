"""From mantle to magma: partially melt a mantle rock, then cool that melt with fractional
crystallization, following Bowen's reaction series and the magma types the melt passes through.

Stage 1 is batch melting of the rock (KLB-1 peridotite by default) at the source pressure and
temperature: the melt there is the primary magma.  Stage 2 takes that melt up to a crustal pressure,
adds water if asked (the mantle rock is dry; the water stands for a hydrous source, as above a
subduction zone), and cools it in steps, taking the crystals out of contact with the melt at every
step.  Masses in stage 2 are grams per 100 g of primary melt.

What limits the felsic end: the MELTS liquid holds sodium too well (ThermoEngine's MELTS does the
same), so a melt this poor in potassium builds up Na2O as it evolves; past about 62 wt% SiO2 (near
800 C with the defaults) its silica falls while Na2O climbs, which is not a real trend.  A crustal
andesite (the "Arc andesite" card) goes on to rhyolite and the granite minimum.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .classify import magma_type, tas_name
from .composition import PRESETS, Composition
from .equilibrium import EquilibriumResult
from .melting import DEFAULT_PHASES, _approach, _scaled
from .path import allowed_phases, iter_crystallize


@dataclass
class Source:
    """The partly molten mantle at the source (grams per 100 g of rock)."""
    result: EquilibriumResult
    primary: Composition                 # the melt, volatile-free, as it leaves the source
    notes: list[str] = field(default_factory=list)


DEFAULTS = dict(t_source=1380.0, p_source=10000.0, p_crust=2000.0, h2o=2.0, fo2="NNO",
                t_end=800.0, step=10.0)


def melt_source(rock: Composition = PRESETS["KLB-1 peridotite"], t_source: float = 1380.0,
                p_source: float = 10000.0, fo2="QFM", phases=DEFAULT_PHASES) -> Source:
    """Batch melting of the rock at (t_source, p_source).  Raises ValueError below the solidus."""
    ph = allowed_phases(phases, p_source)
    r = _approach(rock, t_source, 25.0, p_source, fo2, ph)
    if r is None or r.liquid_grams <= 0 or r.t_c > t_source + 1e-6:
        raise ValueError(f"the rock does not melt at {t_source:.0f} C and {p_source:.0f} bar")
    r = _scaled(r, 100.0)
    if r.t_c != t_source:
        raise ValueError(f"the rock is below its solidus at {t_source:.0f} C and {p_source:.0f} bar")
    return Source(r, r.liquid.normalized(volatile_free=True))


def primary_with_water(primary: Composition, h2o: float) -> Composition:
    """The primary melt with h2o wt% water (of the total)."""
    if not 0.0 <= h2o < 15.0:
        raise ValueError("water must be between 0 and 15 wt%")
    wt = {o: v * (100.0 - h2o) / 100.0 for o, v in primary.normalized(volatile_free=True).wt.items()}
    if h2o > 0:
        wt["H2O"] = h2o
    return Composition(wt, name="primary melt", source="")


def iter_differentiate(primary: Composition, t_start: float, t_end: float = 800.0, step: float = 10.0,
                       p_crust: float = 2000.0, h2o: float = 2.0, fo2="NNO", phases=DEFAULT_PHASES):
    """Stage 2: fractional crystallization of the primary melt.  Yields EquilibriumResult."""
    bulk = primary_with_water(primary, h2o)
    yield from iter_crystallize(bulk, t_start, t_end, step=step, p_bar=p_crust, fo2=fo2,
                                mode="fractional", phases=phases)


@dataclass
class Stage:
    t: float
    melt_grams: float
    sio2: float
    magma_type: str
    tas: str


def magma_type_stages(results) -> list[Stage]:
    """The first step at which the melt reaches each new magma type (by silica)."""
    out: list[Stage] = []
    for r in results:
        mt = magma_type(r.liquid)
        if not out or out[-1].magma_type != mt:
            liq = r.liquid.with_fe3_ratio(0.0).normalized(volatile_free=True)
            out.append(Stage(r.t_c, r.liquid_grams, liq.get("SiO2"), mt, tas_name(r.liquid)))
    return out
