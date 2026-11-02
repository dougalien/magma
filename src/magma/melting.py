"""Partial melting of a rock on heating, batch (equilibrium) or fractional.

Batch melting is equilibrium, so it does not depend on the path: Magma computes it by cooling the
rock's own composition from above its liquidus (where the solver has a sure start) down to the
lowest temperature asked, and reports the steps from low to high temperature.  The steps stop at the
solidus, which Magma cannot cross yet.

Fractional melting starts from the batch state at the lowest temperature and heats in steps.  At each
step the melt above a small retained fraction (the porosity of the residue) is removed and saved;
the residue (solids + retained melt) is the next step's rock.  Reported per step: the residue, the
melt removed in that step, and the pooled (accumulated) melt.  Masses are grams per 100 g of rock.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

from .composition import MOLAR_MASS, Composition
from .equilibrium import EquilibriumResult, Equilibrator
from .path import SOLIDUS_MELT, allowed_phases
from .thermo.phases import OXIDES_OF

DEFAULT_PHASES = ("olivine", "feldspar", "clinopyroxene", "pigeonite", "orthopyroxene", "spinel", "titanomagnetite", "rhombohedral", "garnet", "hornblende", "biotite", "alkali-feldspar", "muscovite", "quartz", "tridymite", "water")
T_LIQUIDUS_START = 1900.0       # degC; the cooling run that gives batch melting starts here at most


@dataclass
class MeltStep:
    result: EquilibriumResult                   # the system at this T (grams per 100 g of rock)
    extracted_grams: float = 0.0                # melt removed at this step (fractional)
    extracted: Composition | None = None        # its composition
    pooled_grams: float = 0.0                   # all melt removed so far
    pooled: Composition | None = None           # its composition
    notes: list[str] = field(default_factory=list)


def solid_oxide_grams(r: EquilibriumResult) -> dict[str, float]:
    out: dict[str, float] = {}
    for s in r.solids:
        for e, n in s.moles.items():
            for o, k in OXIDES_OF[e].items():
                out[o] = out.get(o, 0.0) + k * n * MOLAR_MASS[o]
    return out


def _batch_descending(bulk, t_low, t_high, step, p_bar, fo2, phases, on_progress):
    """Cool from above the liquidus to t_low; return the results at t_high, t_high - step, ..."""
    n = int(round((t_high - t_low) / step))
    targets = [t_high - k * step for k in range(n + 1)]
    out, prev = [], None
    # a coarse approach from T_LIQUIDUS_START gives the first fine step a warm start
    t0 = max(T_LIQUIDUS_START, t_high)
    if t0 > t_high:
        n_c = max(1, int(round((t0 - t_high) / 25.0)))
        dt = (t0 - t_high) / n_c
        for k in range(n_c):
            prev = Equilibrator(bulk, fo2=fo2, phases=phases).step_from(t0 - k * dt, p_bar, prev)
    for k, t in enumerate(targets):
        r = Equilibrator(bulk, fo2=fo2, phases=phases).step_from(t, p_bar, prev)
        if prev is not None and prev.melt_fraction < SOLIDUS_MELT and r.melt_fraction > prev.melt_fraction + 0.002:
            break                                 # the solidus: the assemblage cannot be followed
        out.append(r)
        prev = copy.deepcopy(r)
        if on_progress:
            on_progress(k + 1, len(targets))
    return out


def _approach(bulk, t_target, step, p_bar, fo2, phases):
    """Cool in steps of at most 25 degrees from T_LIQUIDUS_START toward t_target; the last state
    before the solidus (or at t_target)."""
    t0 = max(T_LIQUIDUS_START, t_target)
    n = max(1, int(round((t0 - t_target) / min(25.0, step))))
    prev = None
    for k in range(n + 1):
        t = t0 - k * (t0 - t_target) / n
        r = Equilibrator(bulk, fo2=fo2, phases=phases).step_from(t, p_bar, prev)
        if prev is not None and prev.melt_fraction < SOLIDUS_MELT and r.melt_fraction > prev.melt_fraction + 0.002:
            break
        prev = copy.deepcopy(r)
    return prev


def _scaled(r: EquilibriumResult, total: float) -> EquilibriumResult:
    r = copy.deepcopy(r)
    f = total / (r.liquid_grams + sum(s.grams for s in r.solids))
    r.liquid_grams *= f
    for s in r.solids:
        s.grams *= f
        s.moles = {e: v * f for e, v in s.moles.items()}
    return r


def iter_melt(bulk: Composition, t_start: float, t_end: float, step: float = 10.0, p_bar: float = 1.0,
              fo2="QFM", mode: str = "batch", phases=DEFAULT_PHASES, retained: float = 0.01,
              on_progress=None):
    """Heat a rock from t_start up to t_end (degC).  Yields MeltStep, lowest T first.  In batch mode
    nothing is yielded until the cooling run behind it is done (on_progress(k, n) reports it)."""
    if mode not in ("batch", "fractional"):
        raise ValueError("mode must be 'batch' or 'fractional'")
    if step <= 0 or t_end < t_start:
        raise ValueError("heat from a lower to a higher temperature with a positive step")
    if not 0.0 <= retained < 0.5:
        raise ValueError("the retained melt fraction must be between 0 and 0.5")
    phases = allowed_phases(phases, p_bar)
    total = sum(bulk.normalized().wt.values())    # 100 g of rock
    if mode == "batch":
        steps = _batch_descending(bulk, t_start, t_end, step, p_bar, fo2, phases, on_progress)
        for r in reversed(steps):
            yield MeltStep(_scaled(r, total))
        return

    # fractional: the batch state at t_start (or at the solidus, if that is higher), then heat
    start = _approach(bulk, t_start, step, p_bar, fo2, phases)
    t = start.t_c
    rock_total = 100.0
    state = _scaled(start, rock_total)
    pooled_g, pooled_ox = 0.0, {}
    while t <= t_end + 1e-9:
        r = state
        if r.t_c != t:
            r = Equilibrator(rock, fo2=fo2, phases=phases).step_from(t, p_bar, prev)
            r = _scaled(r, rock_total)
        system = r.liquid_grams + sum(s.grams for s in r.solids)
        keep = retained * system
        ms = MeltStep(copy.deepcopy(r))
        if r.liquid_grams > keep:
            out_g = r.liquid_grams - keep
            liq = r.liquid.normalized()
            for o, v in liq.wt.items():
                pooled_ox[o] = pooled_ox.get(o, 0.0) + v / 100.0 * out_g
            pooled_g += out_g
            ms.extracted_grams, ms.extracted = out_g, liq
            r.liquid_grams = keep
        ms.pooled_grams = pooled_g
        ms.pooled = Composition(dict(pooled_ox)).normalized() if pooled_g > 0 else None
        yield ms
        # the residue (solids + retained melt) is the next rock
        res = solid_oxide_grams(r)
        for o, v in r.liquid.normalized().wt.items():
            res[o] = res.get(o, 0.0) + v / 100.0 * r.liquid_grams
        rock_total = sum(res.values())
        if rock_total <= 1e-6 or not r.solids:
            return                                # all melted
        rock = Composition(res)
        prev = r
        t += step
        state = r
