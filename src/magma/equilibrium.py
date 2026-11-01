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
class _Solution:
    x: np.ndarray
    norm: float
    success: bool


def _newton(residual, y, tol=1e-9, max_iter=60, h=1e-7, stall=False) -> _Solution:
    """Damped Gauss-Newton with a finite-difference Jacobian.  `residual` returns None outside the
    feasible region (a solid exceeding the bulk, or an impossible site occupancy); steps are halved
    until they stay feasible and reduce the residual."""
    r = residual(y)
    if r is None:
        return _Solution(y, np.inf, False)
    norm = float(np.linalg.norm(r))
    history = []
    for it in range(max_iter):
        if norm < tol:
            break
        history.append(norm)
        if stall and it >= 15 and norm > 0.5 * history[-11] and norm > 1e-3:
            break                            # stalled: no real progress in ten iterations
        jac = np.empty((len(r), len(y)))
        for j in range(len(y)):
            dy = h * max(1.0, abs(y[j]))
            yp = y.copy()
            yp[j] += dy
            rp = residual(yp)
            if rp is None:
                yp[j] -= 2 * dy
                rp = residual(yp)
                if rp is None:
                    return _Solution(y, norm, False)
                dy = -dy
            jac[:, j] = (rp - r) / dy
        if not np.all(np.isfinite(jac)):
            return _Solution(y, norm, False)
        try:
            step = np.linalg.lstsq(jac, -r, rcond=None)[0]
        except np.linalg.LinAlgError:
            return _Solution(y, norm, False)
        big = np.max(np.abs(step))
        if big > 1.0:                        # at most a factor e in any amount, 1 in any fraction
            step /= big
        lam, accepted = 1.0, False
        while lam > 1e-8:
            y_new = y + lam * step
            r_new = residual(y_new)
            if r_new is not None and np.linalg.norm(r_new) < (1 - 1e-4 * lam) * norm:
                y, r, norm, accepted = y_new, r_new, float(np.linalg.norm(r_new)), True
                break
            lam *= 0.5
        if not accepted:
            break
    return _Solution(y, norm, norm < 1e-4)          # 1e-4 RT, about 1 J per end member


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
    warm_start: str = ""             # "kept": the start's assemblage solved as it was; "dropped": one less

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

    def __init__(self, bulk: Composition, fo2="QFM", phases=("olivine", "feldspar", "clinopyroxene", "pigeonite", "orthopyroxene", "spinel", "rhombohedral", "garnet", "water")):
        """fo2: a buffer ("QFM+1"), a log10 value, or None for a system closed to oxygen."""
        b = bulk.with_fe3_ratio(0.0)                     # carry iron as total FeO
        self.bulk_fe2o3 = bulk.moles().get("Fe2O3", 0.0)  # used when fo2 is None (closed to oxygen)
        self.bulk_moles = {o: v for o, v in b.moles().items() if o in MODEL_OXIDES or o == "FeO"}
        self.bulk_grams = {o: v * MOLAR_MASS[o] for o, v in self.bulk_moles.items()}
        self.ignored = {o: v for o, v in bulk.wt.items() if o not in MODEL_OXIDES and v > 0}
        self.fo2 = fo2
        absent = {o for o in MODEL_OXIDES if o not in ("Fe2O3", "H2O") and self.bulk_moles.get(o, 0.0) <= 0}
        self.phases: list[SolutionPhase] = []
        for name in phases:
            ph = PHASES[name]
            if hasattr(ph, "restricted"):        # end members that need a missing oxide cannot form
                ph = ph.restricted([e for e in ph.end_members if set(OXIDES_OF[e]) & absent])
            self.phases.append(ph)

    # ---- liquid ------------------------------------------------------------------------------
    def _liquid(self, solid_oxides: dict[str, float], t_c, p) -> Composition:
        mol = dict(self.bulk_moles)                 # iron carried as total FeO
        solid_fe2o3 = solid_oxides.get("Fe2O3", 0.0)
        for o, n in solid_oxides.items():
            if o == "Fe2O3":
                mol["FeO"] = mol.get("FeO", 0.0) - 2 * n
            else:
                mol[o] = mol.get(o, 0.0) - n
        wt = {o: n * MOLAR_MASS[o] for o, n in mol.items()}
        if min(wt.values()) < -1e-9:
            raise ValueError("solids exceed the bulk composition")
        liq = Composition({o: max(v, 0.0) for o, v in wt.items()})
        if self.fo2 is None:      # closed to oxygen: the ferric iron not in solids stays in the liquid
            fe_tot = mol.get("FeO", 0.0)
            fe3 = 2 * (self.bulk_fe2o3 - solid_fe2o3)
            if fe3 < -1e-12:
                raise ValueError("solids hold more ferric iron than the system")
            return liq.with_fe3_ratio(min(max(fe3, 0.0) / fe_tot, 1.0) if fe_tot > 0 else 0.0)
        return liq.with_fe3_ratio(fe3_ratio(liq, t_c, p, self.fo2))

    def _mu_liquid(self, liq: Composition, t_k, p) -> np.ndarray:
        g0 = np.array([liquid_gibbs(c, t_k, p) for c in COMPONENTS])
        return g0 + mu_mixing(to_components(liq).x, t_k)

    # ---- solve one assemblage ----------------------------------------------------------------
    def _solve(self, present, t_c, p, guess, stall=False):
        t = t_c + 273.15
        layout = [(ph, e) for ph in present for e in ph.end_members]

        def unpack(y):
            parts, i = [], 0
            for ph in present:
                parts.append(ph.n_from_params(y[i:i + ph.n_params]))
                i += ph.n_params
            n = np.concatenate(parts)
            solid_ox: dict[str, float] = {}
            for (ph, e), ni in zip(layout, n):
                for o, k in _oxide_vector(e).items():
                    solid_ox[o] = solid_ox.get(o, 0.0) + k * ni
            return n, solid_ox

        def residual(y):
            try:
                n, solid_ox = unpack(y)
            except OverflowError:
                return None
            if not np.all(np.isfinite(n)):
                return None
            try:
                liq = self._liquid(solid_ox, t_c, p)
                for ph, ni in self._split(present, n):
                    if hasattr(ph, "feasible_x") and not ph.feasible_x(ni / ni.sum()):
                        raise ValueError("infeasible solid composition")
            except ValueError:
                return None
            mu_l = self._mu_liquid(liq, t, p)
            res, i = [], 0
            for ph in present:
                m = len(ph.end_members)
                ni = n[i:i + m]
                mu_s = ph.mu(ni / ni.sum(), t, p)
                active = getattr(ph, "active", None)
                for j, e in enumerate(ph.end_members):
                    if active is None or active[j]:
                        res.append((mu_s[j] - mu_in_liquid(e, mu_l)) / (8.3143 * t))
                i += m
            res = np.array(res)
            return res if np.all(np.isfinite(res)) else None

        sol = _newton(residual, np.asarray(guess, float), stall=stall)
        n, solid_ox = unpack(sol.x)
        return sol, n, layout, solid_ox

    @staticmethod
    def _split(present, n):
        i = 0
        for ph in present:
            m = len(ph.end_members)
            yield ph, n[i:i + m]
            i += m

    def step_from(self, t_c: float, p_bar: float, start: "EquilibriumResult | None", depth: int = 0
                  ) -> EquilibriumResult:
        """equilibrate(t_c, start=start), going through halfway temperatures when the start's whole
        assemblage cannot be carried over in one step (large changes near phase-out boundaries)."""
        r = self.equilibrate(t_c, p_bar, start=start)
        if start is None or r.warm_start == "kept" or depth >= 2 or abs(t_c - start.t_c) < 0.5:
            return r
        mid = self.step_from(0.5 * (t_c + start.t_c), p_bar, start, depth + 1)
        return self.equilibrate(t_c, p_bar, start=mid)

    def equilibrate(self, t_c: float, p_bar: float = 1.0, max_iter: int = 20,
                    start: "EquilibriumResult | None" = None) -> EquilibriumResult:
        """start: a result at a nearby T (the previous step of a path); its solids, scaled to this
        bulk, are the first guess.  Without it (or if that fails) solids are added one at a time."""
        t = t_c + 273.15
        present: list[SolutionPhase] = []
        guess = np.array([])
        n = np.array([])
        layout, solid_ox = [], {}
        skipped: set[str] = set()
        warm = ""
        if start is not None and start.solids:
            by_name = {ph.name: ph for ph in self.phases}
            scale = sum(self.bulk_grams.values()) / (start.liquid_grams + sum(s.grams for s in start.solids))
            trial, tguess, grams_k = [], [], []
            for s_ in start.solids:
                if s_.name in by_name and s_.grams > 0:
                    ph = by_name[s_.name]
                    x = np.array([s_.x[e] for e in ph.end_members])
                    trial.append(ph)
                    grams_k.append(s_.grams)
                    tguess.append(ph.params_from(x, max(sum(s_.moles.values()) * scale * 0.98, 1e-8)))
            # the previous assemblage, else the same without one phase (one that is running out)
            # smallest first: the phase that is running out is usually the one with least left
            by_size = sorted(range(len(trial)), key=lambda j: grams_k[j])
            options = [list(range(len(trial)))] + [[k for k in range(len(trial)) if k != j]
                                                  for j in by_size] if len(trial) > 1 else [[0]]
            for i_opt, keep in enumerate(options):
                if not keep:
                    continue
                ph_k = [trial[k] for k in keep]
                sol, tn, tlayout, tsolid = self._solve(ph_k, t_c, p_bar, np.concatenate([tguess[k] for k in keep]),
                                                      stall=True)   # carrying a start over: give up early
                if sol.success and all(nk.sum() > 0 for _, nk in self._split(ph_k, tn)):
                    present, guess, n, layout, solid_ox = ph_k, sol.x, tn, tlayout, tsolid
                    warm = "kept" if i_opt == 0 else "dropped"
                    break
        for _ in range(max_iter):
            liq = self._liquid(solid_ox, t_c, p_bar)
            mu_l = self._mu_liquid(liq, t, p_bar)
            absent = [ph for ph in self.phases if ph not in present and ph.name not in skipped]
            aff = {ph.name: ph.affinity(mu_l, t, p_bar) for ph in absent}
            best = max(aff.items(), key=lambda kv: kv[1][0], default=None)
            if best is None or best[1][0] <= 1e-3:
                break
            ph = next(p for p in absent if p.name == best[0])
            x = best[1][1]
            trial = present + [ph]
            for amount in (1e-3, 1e-2, 3e-2):   # moles of the new phase to start from
                tguess = np.concatenate([guess, ph.params_from(x, amount)])
                sol, tn, tlayout, tsolid = self._solve(trial, t_c, p_bar, tguess)
                if sol.success:
                    break
            if not sol.success:
                # a phase driven toward zero while the new one grows is being replaced by it
                gone = [k for k, (_, nk) in enumerate(self._split(trial, tn)) if nk.sum() < 1e-4 * max(tn.sum(), 1e-30)]
                if gone and len(gone) < len(trial):
                    keep = [k for k in range(len(trial)) if k not in gone]
                    starts = np.cumsum([0] + [q.n_params for q in trial])
                    g2 = np.concatenate([sol.x[starts[k]:starts[k + 1]] for k in keep])
                    trial = [trial[k] for k in keep]
                    sol, tn, tlayout, tsolid = self._solve(trial, t_c, p_bar, g2)
            if not sol.success or ph not in trial:
                skipped.add(ph.name)             # cannot be added at this T; try the others
                continue
            present, guess, n, layout, solid_ox = trial, sol.x, tn, tlayout, tsolid
            # drop phases that have essentially vanished
            keep, i, j, g2 = [], 0, 0, []
            for p_ in present:
                m = len(p_.end_members)
                if n[i:i + m].sum() > 1e-10:
                    keep.append(p_)
                    g2.extend(guess[j:j + p_.n_params])
                i += m
                j += p_.n_params
            if len(keep) != len(present):
                present, guess = keep, np.array(g2)
                if present:
                    sol, n, layout, solid_ox = self._solve(present, t_c, p_bar, guess)
                    guess = sol.x
                else:
                    n, layout, solid_ox = np.array([]), [], {}
        liq = self._liquid(solid_ox, t_c, p_bar)
        mu_l = self._mu_liquid(liq, t, p_bar)
        result = EquilibriumResult(t_c, p_bar, liq.normalized(), liq.total, warm_start=warm)
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
