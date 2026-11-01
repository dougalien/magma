"""Equilibrium between the liquid and solid solutions at fixed T, P, bulk composition and fO2.

The calculation follows the logic of MELTS (Ghiorso 1985, 1994): add the solid with the largest
positive affinity, solve the equilibrium conditions for the current assemblage, drop any solid that
runs out, repeat until no absent phase is oversaturated.

Equilibrium conditions, one per end member of every present solid phase:

    g0_i + d(G_mix)/d(n_i) = sum_k nu_ik mu_liq_k

with the liquid = bulk - solids.  The fO2 is imposed by giving the liquid's iron the Fe3+/Fe2+ of
Kress and Carmichael (1991) at that T, P and fO2 (the solids in use so far carry no ferric iron).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .composition import FEO_PER_FE2O3, MOLAR_MASS, Composition
from .redox import fe3_ratio
from .thermo.liquid import COMPONENTS, MODEL_OXIDES, mu_mixing, to_components
from .thermo.phases import OXIDES_OF, PHASES, SolutionPhase, mu_in_liquid
from .thermo.standard_state import liquid_gibbs


def _oxide_vector(end_member: str) -> dict[str, float]:
    return OXIDES_OF[end_member]


def _grams(oxides: dict[str, float]) -> float:
    return sum(n * MOLAR_MASS[o] for o, n in oxides.items())


@dataclass
class PhaseResult:
    name: str
    grams: float
    moles: dict[str, float]          # end member -> moles
    x: dict[str, float]              # end member -> mole fraction


@dataclass
class EquilibriumResult:
    t_c: float
    p_bar: float
    liquid: Composition              # wt%, normalized, Fe split for the fO2
    liquid_grams: float
    solids: list[PhaseResult] = field(default_factory=list)
    affinities: dict[str, float] = field(default_factory=dict)   # absent phases, J (A < 0: undersaturated)

    @property
    def melt_fraction(self) -> float:
        tot = self.liquid_grams + sum(s.grams for s in self.solids)
        return self.liquid_grams / tot

    def summary(self) -> str:
        lines = [f"T {self.t_c:.1f} C  P {self.p_bar:.0f} bar  melt {100 * self.melt_fraction:.1f} wt%"]
        for s in self.solids:
            comp = " ".join(f"{k} {v:.3f}" for k, v in s.x.items())
            lines.append(f"  {s.name:<10s} {s.grams:8.3f} g  {comp}")
        return "\n".join(lines)


class Equilibrator:
    """Equilibrium of one bulk composition (oxide wt%, i.e. grams per 100 g) at T, P, fO2."""

    def __init__(self, bulk: Composition, fo2="QFM", phases=("olivine", "feldspar")):
        """fo2: a buffer ("QFM+1"), a log10 value, or None for a system closed to oxygen."""
        b = bulk.with_fe3_ratio(0.0)                     # carry iron as total FeO
        self.bulk_fe2o3 = bulk.moles().get("Fe2O3", 0.0)  # used when fo2 is None (closed to oxygen)
        self.bulk_moles = {o: v for o, v in b.moles().items() if o in MODEL_OXIDES or o == "FeO"}
        self.ignored = {o: v for o, v in bulk.wt.items() if o not in MODEL_OXIDES and v > 0}
        self.fo2 = fo2
        self.phases: list[SolutionPhase] = [PHASES[n] for n in phases]

    # ---- liquid ------------------------------------------------------------------------------
    def _liquid(self, solid_oxides: dict[str, float], t_c, p) -> Composition:
        mol = dict(self.bulk_moles)
        for o, n in solid_oxides.items():
            mol[o] = mol.get(o, 0.0) - n
        wt = {o: n * MOLAR_MASS[o] for o, n in mol.items()}
        if min(wt.values()) < -1e-9:
            raise ValueError("solids exceed the bulk composition")
        liq = Composition({o: max(v, 0.0) for o, v in wt.items()})
        if self.fo2 is None:      # closed to oxygen: all ferric iron stays in the liquid
            fe_tot = mol.get("FeO", 0.0)
            return liq.with_fe3_ratio(min(2 * self.bulk_fe2o3 / fe_tot, 1.0) if fe_tot > 0 else 0.0)
        return liq.with_fe3_ratio(fe3_ratio(liq, t_c, p, self.fo2))

    def _mu_liquid(self, liq: Composition, t_k, p) -> np.ndarray:
        g0 = np.array([liquid_gibbs(c, t_k, p) for c in COMPONENTS])
        return g0 + mu_mixing(to_components(liq).x, t_k)

    # ---- solve one assemblage ----------------------------------------------------------------
    def _solve(self, present, t_c, p, guess):
        t = t_c + 273.15
        layout = [(ph, e) for ph in present for e in ph.end_members]

        def unpack(y):
            n = np.exp(y)
            solid_ox: dict[str, float] = {}
            for (ph, e), ni in zip(layout, n):
                for o, k in _oxide_vector(e).items():
                    solid_ox[o] = solid_ox.get(o, 0.0) + k * ni
            return n, solid_ox

        def residual(y):
            n, solid_ox = unpack(y)
            try:
                liq = self._liquid(solid_ox, t_c, p)
            except ValueError:
                return np.full(len(y), 1e6)
            mu_l = self._mu_liquid(liq, t, p)
            res, i = [], 0
            for ph in present:
                m = len(ph.end_members)
                ni = n[i:i + m]
                mu_s = ph.mu(ni / ni.sum(), t, p)
                for j, e in enumerate(ph.end_members):
                    res.append((mu_s[j] - mu_in_liquid(e, mu_l)) / (8.3143 * t))
                i += m
            return np.array(res)

        sol = least_squares(residual, guess, method="lm", xtol=1e-12, ftol=1e-12, max_nfev=4000)
        n, solid_ox = unpack(sol.x)
        return sol, n, layout, solid_ox

    def equilibrate(self, t_c: float, p_bar: float = 1.0, max_iter: int = 20) -> EquilibriumResult:
        t = t_c + 273.15
        present: list[SolutionPhase] = []
        guess = np.array([])
        n = np.array([])
        layout, solid_ox = [], {}
        for _ in range(max_iter):
            liq = self._liquid(solid_ox, t_c, p_bar)
            mu_l = self._mu_liquid(liq, t, p_bar)
            absent = [ph for ph in self.phases if ph not in present]
            aff = {ph.name: ph.affinity(mu_l, t, p_bar) for ph in absent}
            best = max(aff.items(), key=lambda kv: kv[1][0], default=None)
            if best is None or best[1][0] <= 1e-3:
                break
            ph = next(p for p in absent if p.name == best[0])
            x = best[1][1]
            present.append(ph)
            guess = np.concatenate([guess, np.log(np.maximum(x * 1e-3, 1e-12))])
            sol, n, layout, solid_ox = self._solve(present, t_c, p_bar, guess)
            guess = sol.x
            # drop phases that have essentially vanished
            keep, i, g2 = [], 0, []
            for p_ in present:
                m = len(p_.end_members)
                if n[i:i + m].sum() > 1e-10:
                    keep.append(p_)
                    g2.extend(guess[i:i + m])
                i += m
            if len(keep) != len(present):
                present, guess = keep, np.array(g2)
                if present:
                    sol, n, layout, solid_ox = self._solve(present, t_c, p_bar, guess)
                    guess = sol.x
                else:
                    n, layout, solid_ox = np.array([]), [], {}
        liq = self._liquid(solid_ox, t_c, p_bar)
        mu_l = self._mu_liquid(liq, t, p_bar)
        result = EquilibriumResult(t_c, p_bar, liq.normalized(), liq.total)
        i = 0
        for ph in present:
            m = len(ph.end_members)
            ni = n[i:i + m]
            moles = dict(zip(ph.end_members, ni.tolist()))
            grams = sum(_grams({o: k * v for o, k in _oxide_vector(e).items()}) for e, v in moles.items())
            result.solids.append(PhaseResult(ph.name, grams, moles,
                                             {e: v / ni.sum() for e, v in moles.items()}))
            i += m
        for ph in self.phases:
            if ph not in present:
                result.affinities[ph.name] = ph.affinity(mu_l, t, p_bar)[0]
        return result
