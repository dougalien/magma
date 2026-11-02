"""The Holland & Powell engine behind the app: wt% in, cooling and melting paths out.

This adapter turns a ``Composition`` (oxide wt%) into the eleven components of the solver, runs
``minimize.equilibrate`` step by step (each step starting from the one before), and returns results
with the same fields as the MELTS engine's ``EquilibriumResult``, so the app shows either.

Phase names: the HP phases are reported under the app's mineral names (olivine, feldspar,
clinopyroxene, ...).  One HP model can stand for two app minerals: a feldspar with Or above 30% is
K-feldspar, a clinopyroxene with Wo below 25% is pigeonite, a spinel with more than half magnetite and
ulvospinel is titanomagnetite.  Phases that MELTS does not have keep names of their own (cordierite,
epidote, chlorite, rutile, sphene, sillimanite, ...).

Iron: FeO (total) and O (the oxygen in excess of FeO, i.e. moles of Fe2O3).  With an oxygen buffer the
system is open to oxygen at that fO2; with fo2=None it is closed and keeps the Fe2O3 it started with.
MnO, P2O5, NiO, CoO and CO2 are not in the HP igneous set and are left out.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

import numpy as np

from ..composition import MOLAR_MASS, Composition
from ..redox import log_fo2 as _log_fo2
from .minimize import COMPONENTS, MELT_MIN_T, PURE_PHASES, Result, equilibrate

_M = np.array([MOLAR_MASS["FeO"] if c == "FeO" else 15.9994 if c == "O" else MOLAR_MASS[c] for c in COMPONENTS])
_I = {c: i for i, c in enumerate(COMPONENTS)}

# app mineral (the picker) -> HP phases it switches on
APP_TO_HP = {
    "olivine": ("ol",), "feldspar": ("fsp",), "alkali-feldspar": ("fsp",), "clinopyroxene": ("cpx",),
    "pigeonite": ("cpx",), "orthopyroxene": ("opx",), "spinel": ("spl",), "titanomagnetite": ("spl",),
    "rhombohedral": ("ilm",), "garnet": ("g",), "hornblende": ("amp",), "biotite": ("bi",),
    "muscovite": ("mu",), "quartz": ("q", "coe", "stv"), "tridymite": ("trd", "crst"), "water": ("fl",),
}
# always on: HP phases the picker does not list
ALWAYS = ("liq", "cd", "ep", "chl", "ky", "sill", "and", "ru", "sph", "cor", "ne")

HP_NAME = {"ol": "olivine", "opx": "orthopyroxene", "ilm": "rhombohedral", "g": "garnet", "amp": "hornblende",
           "bi": "biotite", "mu": "muscovite", "fl": "water", "cd": "cordierite", "ep": "epidote",
           "chl": "chlorite", "q": "quartz", "trd": "tridymite", "crst": "cristobalite", "coe": "coesite",
           "stv": "stishovite", "ky": "kyanite", "sill": "sillimanite", "and": "andalusite", "ru": "rutile",
           "sph": "sphene", "cor": "corundum", "ne": "nepheline"}
LABEL = {"cordierite": "Cordierite", "epidote": "Epidote", "chlorite": "Chlorite", "cristobalite": "Cristobalite",
         "coesite": "Coesite", "stishovite": "Stishovite", "kyanite": "Kyanite", "sillimanite": "Sillimanite",
         "andalusite": "Andalusite", "rutile": "Rutile", "sphene": "Sphene (titanite)", "corundum": "Corundum",
         "nepheline": "Nepheline"}
# anhydrous oxygens per formula unit, for the structural formula
OXYGENS = {"ol": 4, "opx": 6, "cpx": 6, "fsp": 8, "spl": 4, "ilm": 3, "g": 12, "amp": 23, "bi": 11, "mu": 11,
           "cd": 18, "ep": 12.5, "chl": 14}
_CATION = {"SiO2": ("Si", 1, 2), "TiO2": ("Ti", 1, 2), "Al2O3": ("Al", 2, 3), "Cr2O3": ("Cr", 2, 3),
           "MgO": ("Mg", 1, 1), "CaO": ("Ca", 1, 1), "Na2O": ("Na", 2, 1), "K2O": ("K", 2, 1)}


def to_components(bulk: Composition) -> np.ndarray:
    """Moles of COMPONENTS in 100 g of the composition (normalized)."""
    n = bulk.normalized().moles()
    out = np.zeros(len(COMPONENTS))
    for c in COMPONENTS:
        if c in ("FeO", "O"):
            continue
        out[_I[c]] = n.get(c, 0.0)
    out[_I["FeO"]] = n.get("FeO", 0.0) + 2.0 * n.get("Fe2O3", 0.0)
    out[_I["O"]] = n.get("Fe2O3", 0.0)
    return out


def to_composition(moles: np.ndarray, name: str = "") -> Composition:
    """Oxide wt% (normalized) of moles of COMPONENTS, with Fe3+ as Fe2O3."""
    wt = {}
    for c in COMPONENTS:
        if c in ("FeO", "O"):
            continue
        if moles[_I[c]] > 1e-12:
            wt[c] = moles[_I[c]] * MOLAR_MASS[c]
    fe2o3 = max(moles[_I["O"]], 0.0)
    feo = max(moles[_I["FeO"]] - 2.0 * fe2o3, 0.0)
    if feo > 1e-12:
        wt["FeO"] = feo * MOLAR_MASS["FeO"]
    if fe2o3 > 1e-12:
        wt["Fe2O3"] = fe2o3 * MOLAR_MASS["Fe2O3"]
    return Composition(wt, name).normalized()


def grams(moles: np.ndarray) -> float:
    return float(moles @ _M)


@dataclass
class HPPhase:
    """One phase, in the shape the app expects (name, grams, moles, x), plus its formula and numbers."""
    name: str
    grams: float
    moles: dict[str, float]
    x: dict[str, float]
    formula: str = ""
    numbers: dict[str, float] = field(default_factory=dict)
    hp: str = ""                         # the HP phase it came from
    oxide_moles: np.ndarray | None = None


@dataclass
class HPResult:
    t_c: float
    p_bar: float
    liquid: Composition | None
    liquid_grams: float
    solids: list = field(default_factory=list)
    affinities: dict = field(default_factory=dict)
    warm_start: str = ""
    raw: Result | None = None            # the solver's answer (moles, gamma), for the next step
    liquid_moles: np.ndarray | None = None
    log_fo2: float | None = None

    @property
    def melt_fraction(self) -> float:
        tot = self.liquid_grams + sum(s.grams for s in self.solids)
        return self.liquid_grams / tot if tot > 0 else 0.0


def hp_phases(app_phases) -> tuple:
    out = list(ALWAYS)
    for p in app_phases:
        for h in APP_TO_HP.get(p, ()):
            if h not in out:
                out.append(h)
    return tuple(out)


def _formula(hp: str, c: np.ndarray) -> str:
    """Cations per OXYGENS[hp] anhydrous oxygens, e.g. Ca0.85Mg0.90Fe0.21Al0.08Si1.96O6."""
    if hp not in OXYGENS:
        return ""
    fe3 = 2.0 * max(c[_I["O"]], 0.0)
    fe2 = max(c[_I["FeO"]] - fe3, 0.0)
    cat, oxy = {}, 0.0
    for ox, (el, k, o) in _CATION.items():
        v = c[_I[ox]]
        if v > 0:
            cat[el] = k * v
            oxy += o * v
    cat["Fe"] = fe2 + fe3
    oxy += fe2 + 1.5 * fe3
    if oxy <= 0:
        return ""
    f = OXYGENS[hp] / oxy
    order = ("K", "Na", "Ca", "Mg", "Fe", "Ti", "Cr", "Al", "Si")
    parts = [f"{el}{cat[el] * f:.2f}" for el in order if cat.get(el, 0) * f >= 0.005]
    o = OXYGENS[hp]
    tail = {"amp": "O22(OH)2", "bi": "O10(OH)2", "mu": "O10(OH)2", "chl": "O10(OH)8", "ep": "O12(OH)"}.get(
        hp, f"O{o:g}")
    return "".join(parts) + tail


def _numbers(hp: str, c: np.ndarray, x: dict) -> dict:
    fe3 = 2.0 * max(c[_I["O"]], 0.0)
    fe2 = max(c[_I["FeO"]] - fe3, 0.0)
    mg, ca = c[_I["MgO"]], c[_I["CaO"]]
    mgn = 100 * mg / (mg + fe2) if mg + fe2 > 0 else 0.0
    if hp == "ol":
        return {"Fo": mgn}
    if hp == "fsp":
        na, k = 2 * c[_I["Na2O"]], 2 * c[_I["K2O"]]
        s = ca + na + k
        return {"An": 100 * ca / s, "Or": 100 * k / s} if s > 0 else {}
    if hp in ("cpx", "opx"):
        fe = c[_I["FeO"]]
        return {"Wo": 100 * ca / (ca + mg + fe), "Mg#": mgn}
    if hp == "spl":
        al, cr = 2 * c[_I["Al2O3"]], 2 * c[_I["Cr2O3"]]
        out = {"Usp": 100 * x.get("usp", 0.0), "Mg#": mgn}
        if cr > 1e-6:
            out["Cr#"] = 100 * cr / (cr + al)
        return out
    if hp == "ilm":
        return {"Ilm": 100 * (1 - x.get("hm", 0.0)), "Hem": 100 * x.get("hm", 0.0)}
    if hp in ("amp", "bi", "g", "cd", "chl"):
        return {"Mg#": mgn}
    return {}


def _app_name(hp: str, nums: dict, x: dict) -> str:
    if hp == "fsp":
        return "alkali-feldspar" if nums.get("Or", 0) > 30 else "feldspar"
    if hp == "cpx":
        return "pigeonite" if nums.get("Wo", 100) < 25 else "clinopyroxene"
    if hp == "spl":
        return "titanomagnetite" if sum(x.get(e, 0) for e in ("nmt", "imt", "usp")) > 0.5 else "spinel"
    return HP_NAME.get(hp, hp)


def _wrap(r: Result, lf: float | None) -> HPResult:
    t_c, p_bar = round(r.t - 273.15, 6), round(r.p * 1000.0, 6)
    liq_m = np.zeros(len(COMPONENTS))
    solids = []
    for ph in r.phases:
        m = ph.oxide_moles
        if ph.name == "liq":
            liq_m = liq_m + m
            continue
        if ph.prop is None:
            x, mol = {ph.name: 1.0}, {ph.name: ph.amount}
        else:
            x = {e: float(v) for e, v in zip(ph.end_members, ph.prop)}
            mol = {e: float(v * ph.amount) for e, v in zip(ph.end_members, ph.prop)}
        nums = _numbers(ph.name, ph.composition, x)
        name = _app_name(ph.name, nums, x)
        solids.append(HPPhase(name, grams(m), mol, x, _formula(ph.name, ph.composition), nums, ph.name, m))
    # one row per app mineral (two amphiboles across a solvus are summed; the larger gives the formula)
    merged: dict[str, HPPhase] = {}
    for s in sorted(solids, key=lambda s: -s.grams):
        if s.name in merged:
            m = merged[s.name]
            m.grams += s.grams
            m.oxide_moles = m.oxide_moles + s.oxide_moles
            for e, v in s.moles.items():
                m.moles[e] = m.moles.get(e, 0.0) + v
        else:
            merged[s.name] = copy.copy(s)
    lg = grams(liq_m)
    liquid = to_composition(liq_m, "liquid") if lg > 1e-9 else None
    return HPResult(t_c, p_bar, liquid, lg, list(merged.values()), raw=r, liquid_moles=liq_m, log_fo2=lf)


def equilibrium(bulk_moles: np.ndarray, t_c: float, p_bar: float, fo2="QFM", phases=None,
                previous: HPResult | None = None) -> HPResult:
    """One equilibrium.  fo2: a buffer ('QFM+1'), a log10 value, or None (closed to oxygen)."""
    lf = None if fo2 is None else _log_fo2(fo2, t_c, p_bar)
    r = equilibrate(bulk_moles, t_c + 273.15, p_bar / 1000.0, phases=phases, log_fo2=lf,
                    previous=previous.raw if previous is not None else None)
    return _wrap(r, lf)


def iter_crystallize(bulk: Composition, t_start: float, t_end: float, step: float = 5.0, p_bar: float = 1.0,
                     fo2="QFM", mode: str = "equilibrium", phases=tuple(APP_TO_HP)):
    """Cool from t_start to t_end (degC), grams per 100 g of starting melt.  Fractional mode takes the
    solids (and any fluid) out after each step.  The path ends when the melt is gone."""
    if mode not in ("equilibrium", "fractional"):
        raise ValueError("mode must be 'equilibrium' or 'fractional'")
    if step <= 0 or t_end > t_start:
        raise ValueError("cool from a higher to a lower temperature with a positive step")
    hp = hp_phases(phases)
    b = to_components(bulk) * 100.0 / grams(to_components(bulk))
    prev = None
    n = int(round((t_start - t_end) / step))
    for k in range(n + 1):
        t = t_start - k * step
        if t + 273.15 < MELT_MIN_T:
            return
        r = equilibrium(b, t, p_bar, fo2, hp, prev)
        if r.liquid is None:
            return                                   # below the solidus
        yield r
        prev = r          # in fractional mode too: the next step re-solves its assemblage for the new bulk
        b = r.liquid_moles.copy() if mode == "fractional" else r.raw.bulk.copy()   # (after any O exchange)


@dataclass
class MeltStep:
    result: HPResult
    extracted_grams: float = 0.0
    extracted: Composition | None = None
    pooled_grams: float = 0.0
    pooled: Composition | None = None
    notes: list = field(default_factory=list)


def iter_melt(bulk: Composition, t_start: float, t_end: float, step: float = 10.0, p_bar: float = 1.0,
              fo2="QFM", mode: str = "batch", phases=tuple(APP_TO_HP), retained: float = 0.01,
              on_progress=None):
    """Heat a rock from t_start to t_end (degC), grams per 100 g of rock; yields MeltStep from the first
    step with melt.  Batch: the rock and its melt stay together.  Fractional: melt above `retained` of
    the system is taken out at each step and pooled."""
    if mode not in ("batch", "fractional"):
        raise ValueError("mode must be 'batch' or 'fractional'")
    if step <= 0 or t_end < t_start:
        raise ValueError("heat from a lower to a higher temperature with a positive step")
    if not 0.0 <= retained < 0.5:
        raise ValueError("the retained melt fraction must be between 0 and 0.5")
    hp = hp_phases(phases)
    b = to_components(bulk) * 100.0 / grams(to_components(bulk))
    prev = None
    pooled = np.zeros(len(COMPONENTS))
    n = int(round((t_end - t_start) / step))
    for k in range(n + 1):
        t = t_start + k * step
        r = equilibrium(b, t, p_bar, fo2, hp, prev)
        if on_progress:
            on_progress(k + 1, n + 1)
        prev = r
        b = r.raw.bulk.copy()
        if r.liquid is None:
            continue                                 # below the solidus: nothing to show yet
        ms = MeltStep(copy.copy(r))
        if mode == "fractional":
            system = r.liquid_grams + sum(s.grams for s in r.solids)
            keep = retained * system
            if r.liquid_grams > keep:
                f = 1.0 - keep / r.liquid_grams
                out = r.liquid_moles * f
                pooled += out
                ms.extracted_grams, ms.extracted = grams(out), r.liquid
                b = b - out
                ms.result.liquid_grams = keep
                ms.result.liquid_moles = r.liquid_moles - out
        if pooled.any():
            ms.pooled_grams, ms.pooled = grams(pooled), to_composition(pooled, "pooled melt")
        yield ms
        if not r.solids:
            return                                   # all molten


def melt_source(rock: Composition, t_source: float, p_source: float, fo2="QFM", phases=tuple(APP_TO_HP),
                step: float = 50.0):
    """Batch melting at the source: heat the rock from 1000 C in steps to t_source.  Returns the result
    (grams per 100 g of rock) and the melt, volatile-free."""
    t0 = min(1000.0, t_source)
    n = max(1, int(np.ceil((t_source - t0) / step)))
    last = None
    for ms in iter_melt(rock, t0, t_source, step=(t_source - t0) / n if t_source > t0 else 1.0,
                        p_bar=p_source, fo2=fo2, mode="batch", phases=phases):
        last = ms.result
    if last is None or abs(last.t_c - t_source) > 1e-6 or last.liquid is None:
        raise ValueError(f"the rock does not melt at {t_source:.0f} C and {p_source:.0f} bar")
    return last, last.liquid.normalized(volatile_free=True)
