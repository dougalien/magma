"""The late minerals of Bowen's series: quartz, tridymite, alkali feldspar, muscovite, biotite and
hornblende, as in MELTS v1.0.2 (ThermoEngine EquilibrateUsingMELTSv102.m lists all of them).

Parameter values are from ThermoEngine v2 (src/BermanStoichiometricPhases.m, BiotiteBerman.m,
HornblendeBerman.m; AGPL-3), numbers only, with Douglas's approval (2026-11-01 policy, extended to
these phases on 2026-11-02).  The code is our own, written from the equations.

* Quartz: Berman (1988).  Below its alpha-beta transition (848 K at 1 bar, + 0.0237 K/bar) alpha quartz
  with the lambda heat capacity, above it beta quartz; MELTS adds -1291 J to both (QUARTZ_ADJUSTMENT,
  on by default in ThermoEngine).
* Tridymite: beta tridymite (Berman 1988), stable form above 383 K.
* Muscovite: Berman (1988), a pure phase.
* Biotite: annite-phlogopite, ideal on one site plus a regular W = 3.21 kcal (Sack and Ghiorso 1989
  as used in MELTS BiotiteBerman.m).
* Hornblende: pargasite - ferropargasite - magnesiohastingsite, a reciprocal solution with Fe-Mg on
  four M1,2 sites and Fe3+-Al on one M3 site (MELTS HornblendeBerman.m), W(Fe-Mg) = 4 x 1.68 kcal,
  W(Fe3+-Al) = 16.78 kJ.  MELTS's own documentation calls its hornblende poorly constrained.
* Alkali feldspar: a second copy of the plagioclase model (Elkins and Grove 1990), as MELTS has
  Plagioclase and Sanidine; each copy looks for its saturating composition on its own side of the
  feldspar solvus (Or below or above OR_SPLIT).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize

from .minerals import MINERALS
from .phases import R, SolutionPhase, g_ex_feldspar, mu_in_liquid

# ---- pure phases ---------------------------------------------------------------------------------
QUARTZ_ADJUSTMENT = -1291.0
_QZ_L1, _QZ_L2, _QZ_TT, _QZ_DTDP, _QZ_TRL = -9.187e-2, 24.607e-5, 848.0, 0.0237, 373.0


def quartz_gibbs(t: float, p: float) -> float:
    """G of quartz, J/mol (Berman 1988 with the MELTS adjustment)."""
    cp_t = _QZ_TT + _QZ_DTDP * (p - 1.0)
    if t > cp_t:
        return MINERALS["quartz_beta"].gibbs(t, p) + QUARTZ_ADJUSTMENT
    l1, l2 = _QZ_L1, _QZ_L2
    d = _QZ_TT - cp_t
    t0 = _QZ_TRL - d
    x1 = l1 * l1 * d + 2 * l1 * l2 * d * d + l2 * l2 * d ** 3
    x2 = l1 * l1 + 4 * l1 * l2 * d + 3 * l2 * l2 * d * d
    x3 = 2 * l1 * l2 + 3 * l2 * l2 * d
    x4 = l2 * l2
    h = x1 * (t - t0) + x2 * (t * t - t0 * t0) / 2 + x3 * (t ** 3 - t0 ** 3) / 3 + x4 * (t ** 4 - t0 ** 4) / 4
    s = x1 * (math.log(t) - math.log(t0)) + x2 * (t - t0) + x3 * (t * t - t0 * t0) / 2 + x4 * (t ** 3 - t0 ** 3) / 3
    return MINERALS["quartz_alpha"].gibbs(t, p) + h - t * s + QUARTZ_ADJUSTMENT


class PurePhase:
    """A phase of fixed composition (one end member)."""
    n_params = 1

    def __init__(self, name: str, end_member: str, gibbs=None, formula: str = ""):
        self.name, self.end_members = name, (end_member,)
        self._gibbs = gibbs or MINERALS[end_member].gibbs
        self._formula = formula or MINERALS[end_member].formula

    def g0(self, t, p) -> np.ndarray:
        return np.array([self._gibbs(float(t), float(p))])

    def mu(self, x, t, p) -> np.ndarray:
        return self.g0(t, p)

    def gibbs(self, x, t, p) -> float:
        return float(self.g0(t, p)[0])

    def g_mix(self, x, t, p) -> float:
        return 0.0

    def n_from_params(self, prm) -> np.ndarray:
        return np.array([math.exp(prm[0])])

    def params_from(self, x, n_total) -> np.ndarray:
        return np.array([math.log(n_total)])

    def affinity(self, mu_liq, t, p):
        m = mu_in_liquid(self.end_members[0], mu_liq)
        if not np.isfinite(m):
            return -np.inf, np.array([1.0])
        return m - float(self.g0(t, p)[0]), np.array([1.0])

    def formula(self, x) -> str:
        return self._formula


QUARTZ = PurePhase("quartz", "quartz_beta", quartz_gibbs, "SiO2")
TRIDYMITE = PurePhase("tridymite", "tridymite")
MUSCOVITE = PurePhase("muscovite", "muscovite")


# ---- biotite ---------------------------------------------------------------------------------------
W_BIO = 4.184 * 1000.0 * 3.21


def g_ex_biotite(x, t, p):
    return W_BIO * x[0] * x[1]


class Biotite(SolutionPhase):
    def formula(self, x) -> str:
        return f"Phl{100 * x[1]:.0f} Ann{100 * x[0]:.0f}"


BIOTITE = Biotite("biotite", ("annite", "phlogopite"), g_ex_biotite, sites=1.0)


# ---- hornblende ------------------------------------------------------------------------------------
W_FEMG = 4.0 * 1.68 * 4.184 * 1000.0
W_FEAL = 16.78 * 1000.0


def _xl(v):
    return v * math.log(v) if v > 0 else 0.0


class Hornblende(SolutionPhase):
    """Order: pargasite, ferropargasite, magnesiohastingsite; r0 = Fe2+ on M1,2, r1 = Fe3+ on M3."""

    @staticmethod
    def _r(x):
        x = np.asarray(x, float)
        s = x.sum()
        return x[1] / s, x[2] / s

    def g_mix(self, x, t, p) -> float:
        r0, r1 = self._r(x)
        return (W_FEMG * r0 * (1 - r0) + W_FEAL * r1 * (1 - r1)
                + R * t * (4 * (_xl(r0) + _xl(1 - r0)) + _xl(r1) + _xl(1 - r1)))

    def mu(self, x, t, p, h=None) -> np.ndarray:
        r0, r1 = self._r(x)
        e = 1e-300
        g = self.g_mix(x, t, p)
        d0 = W_FEMG * (1 - 2 * r0) + 4 * R * t * (math.log(max(r0, e)) - math.log(max(1 - r0, e)))
        d1 = W_FEAL * (1 - 2 * r1) + R * t * (math.log(max(r1, e)) - math.log(max(1 - r1, e)))
        re = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        return self.g0(t, p) + g + (re - np.array([r0, r1])) @ np.array([d0, d1])

    def formula(self, x) -> str:
        r0, r1 = self._r(x)
        return f"Mg#{100 * (1 - r0):.0f} Fe3+ {r1:.2f}"


HORNBLENDE = Hornblende("hornblende", ("pargasite", "ferropargasite", "magnesiohastingsite"), None)


# ---- the two feldspars -----------------------------------------------------------------------------
OR_SPLIT = 0.30      # mole fraction orthoclase between plagioclase (below) and alkali feldspar (above)


class Feldspar(SolutionPhase):
    """The Elkins and Grove (1990) ternary feldspar, searched on one side of the solvus."""

    def __init__(self, name, or_range):
        super().__init__(name, ("albite", "anorthite", "sanidine"), g_ex_feldspar)
        object.__setattr__(self, "or_range", or_range)

    def affinity(self, mu_liq, t, p):
        a, x = SolutionPhase.affinity(self, mu_liq, t, p)
        lo, hi = self.or_range
        if lo - 1e-9 <= x[2] <= hi + 1e-9 or not np.isfinite(a):
            return a, x
        # the unconstrained minimum is on the other side of the solvus: look on this side only
        d = self.g0(t, p) - np.array([mu_in_liquid(e, mu_liq) for e in self.end_members])
        if not np.all(np.isfinite(d)):
            return -np.inf, x
        eps = 1e-9

        def f(y):
            v = np.clip(np.array([y[0], y[1], 1 - y[0] - y[1]]), eps, 1)
            return float(v @ d) + self.g_mix(v, t, p)

        cons = [{"type": "ineq", "fun": lambda y: 1 - y[0] - y[1] - lo},
                {"type": "ineq", "fun": lambda y: hi - (1 - y[0] - y[1])}]
        best = None
        for orx in (lo + 0.02, 0.5 * (lo + hi), hi - 0.02):
            orx = min(max(orx, lo + 1e-4), hi - 1e-4)
            y0 = np.array([0.8 * (1 - orx), 0.2 * (1 - orx)])
            res = minimize(f, y0, method="SLSQP", bounds=[(eps, 1), (eps, 1)], constraints=cons,
                           options=dict(maxiter=300, ftol=1e-12))
            if np.all(np.isfinite(res.x)) and (best is None or res.fun < best.fun):
                best = res
        if best is None:
            return -np.inf, x
        y = best.x
        v = np.clip(np.array([y[0], y[1], 1 - y[0] - y[1]]), eps, 1)
        return -f(y), v / v.sum()

    def formula(self, x) -> str:
        return f"An{100 * x[1]:.0f} Ab{100 * x[0]:.0f} Or{100 * x[2]:.0f}"


ALKALI_FELDSPAR = Feldspar("alkali-feldspar", (OR_SPLIT, 1.0))
