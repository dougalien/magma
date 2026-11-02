"""Cooling paths: equilibrium or fractional crystallization in temperature steps."""
from __future__ import annotations

import copy

from .composition import Composition
from .equilibrium import Equilibrator, EquilibriumResult


PIGEONITE_MAX_P = 15000.0   # bar; above this pigeonite is left out (Magma's choice: it is a low-pressure
                            # phase, and the second cpx would otherwise take the place of augite + opx)
SOLIDUS_MELT = 0.06     # below this melt fraction a path that loses its assemblage has reached the solidus


def crystallize(bulk: Composition, t_start: float, t_end: float, step: float = 5.0,
                p_bar: float = 1.0, fo2="QFM", mode: str = "equilibrium",
                phases=("olivine", "feldspar", "clinopyroxene", "pigeonite", "orthopyroxene", "spinel", "titanomagnetite", "rhombohedral", "garnet", "water")) -> list[EquilibriumResult]:
    """Cool from t_start to t_end (degC).  mode "fractional" removes the solids after each step
    and continues with the remaining liquid.  Masses are grams per 100 g of the starting melt."""
    return list(iter_crystallize(bulk, t_start, t_end, step, p_bar, fo2, mode, phases))


def iter_crystallize(bulk: Composition, t_start: float, t_end: float, step: float = 5.0,
                     p_bar: float = 1.0, fo2="QFM", mode: str = "equilibrium",
                     phases=("olivine", "feldspar", "clinopyroxene", "pigeonite", "orthopyroxene", "spinel", "titanomagnetite", "rhombohedral", "garnet", "water")):
    """Same as crystallize, one step at a time (the app shows progress and can stop between steps).
    In equilibrium mode the path ends early near the solidus: once the melt is nearly gone and the
    next step cannot keep the assemblage (Magma has no liquid-free equilibrium yet)."""
    if mode not in ("equilibrium", "fractional"):
        raise ValueError("mode must be 'equilibrium' or 'fractional'")
    if step <= 0 or t_end > t_start:
        raise ValueError("cool from a higher to a lower temperature with a positive step")
    phases = allowed_phases(phases, p_bar)
    current, scale, prev = bulk, 100.0, None
    n = int(round((t_start - t_end) / step))
    for k in range(n + 1):
        t = t_start - k * step
        r = Equilibrator(current, fo2=fo2, phases=phases).step_from(t, p_bar, prev)
        if prev is not None and prev.melt_fraction < SOLIDUS_MELT and r.melt_fraction > prev.melt_fraction + 0.002:
            return                               # the previous assemblage could not be followed: solidus
        prev = copy.deepcopy(r) if mode == "equilibrium" else None
        f = scale / (r.liquid_grams + sum(s.grams for s in r.solids))
        r.liquid_grams *= f
        for s in r.solids:
            s.grams *= f
        yield r
        if mode == "fractional" and r.solids:
            scale = r.liquid_grams
            current = r.liquid


def allowed_phases(phases, p_bar):
    return tuple(ph for ph in phases if not (ph == "pigeonite" and p_bar > PIGEONITE_MAX_P))
