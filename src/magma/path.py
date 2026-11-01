"""Cooling paths: equilibrium or fractional crystallization in temperature steps."""
from __future__ import annotations

from .composition import Composition
from .equilibrium import Equilibrator, EquilibriumResult


def crystallize(bulk: Composition, t_start: float, t_end: float, step: float = 5.0,
                p_bar: float = 1.0, fo2="QFM", mode: str = "equilibrium",
                phases=("olivine", "feldspar", "clinopyroxene")) -> list[EquilibriumResult]:
    """Cool from t_start to t_end (degC).  mode "fractional" removes the solids after each step
    and continues with the remaining liquid.  Masses are grams per 100 g of the starting melt."""
    return list(iter_crystallize(bulk, t_start, t_end, step, p_bar, fo2, mode, phases))


def iter_crystallize(bulk: Composition, t_start: float, t_end: float, step: float = 5.0,
                     p_bar: float = 1.0, fo2="QFM", mode: str = "equilibrium",
                     phases=("olivine", "feldspar", "clinopyroxene")):
    """Same as crystallize, one step at a time (the app shows progress and can stop between steps)."""
    if mode not in ("equilibrium", "fractional"):
        raise ValueError("mode must be 'equilibrium' or 'fractional'")
    if step <= 0 or t_end > t_start:
        raise ValueError("cool from a higher to a lower temperature with a positive step")
    current, scale = bulk, 100.0
    n = int(round((t_start - t_end) / step))
    for k in range(n + 1):
        t = t_start - k * step
        r = Equilibrator(current, fo2=fo2, phases=phases).equilibrate(t, p_bar)
        f = scale / (r.liquid_grams + sum(s.grams for s in r.solids))
        r.liquid_grams *= f
        for s in r.solids:
            s.grams *= f
        yield r
        if mode == "fractional" and r.solids:
            scale = r.liquid_grams
            current = r.liquid
