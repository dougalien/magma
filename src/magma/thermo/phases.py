"""Solid solution phases and their affinity for crystallizing from a liquid.

Each phase is a set of end members (MELTS values, minerals.py) and a molar Gibbs energy of mixing.
The affinity of a phase for crystallizing from a liquid with chemical potentials mu_liq is

    A = -min over x of [ sum_i x_i (g0_i - mu~_i) + G_mix(x) ],   mu~_i = sum_k nu_ik mu_liq_k

where nu_ik writes end member i in liquid components (Table 3 of Ghiorso and Sack 1995).  A > 0:
the phase is oversaturated; the x at the minimum is the composition that would crystallize.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import minimize

from .liquid import MODEL_OXIDES, _TO_COMPONENTS
from .minerals import MINERALS

R = 8.3143

OXIDES_OF = {
    "forsterite": {"MgO": 2, "SiO2": 1},
    "fayalite": {"FeO": 2, "SiO2": 1},
    "albite": {"Na2O": 0.5, "Al2O3": 0.5, "SiO2": 3},
    "anorthite": {"CaO": 1, "Al2O3": 1, "SiO2": 2},
    "sanidine": {"K2O": 0.5, "Al2O3": 0.5, "SiO2": 3},
    "diopside": {"CaO": 1, "MgO": 1, "SiO2": 2},
    "clinoenstatite": {"MgO": 2, "SiO2": 2},
    "hedenbergite": {"CaO": 1, "FeO": 1, "SiO2": 2},
    "alumino-buffonite": {"CaO": 1, "TiO2": 0.5, "MgO": 0.5, "Al2O3": 0.5, "SiO2": 1},
    "buffonite": {"CaO": 1, "TiO2": 0.5, "MgO": 0.5, "Fe2O3": 0.5, "SiO2": 1},
    "essenite": {"CaO": 1, "Fe2O3": 0.5, "Al2O3": 0.5, "SiO2": 1},
    "jadeite": {"Na2O": 0.5, "Al2O3": 0.5, "SiO2": 2},
    "chromite": {"FeO": 1, "Cr2O3": 1},
    "hercynite": {"FeO": 1, "Al2O3": 1},
    "magnetite": {"FeO": 1, "Fe2O3": 1},
    "spinel": {"MgO": 1, "Al2O3": 1},
    "ulvospinel": {"FeO": 2, "TiO2": 1},
    "geikielite": {"MgO": 1, "TiO2": 1},
    "hematite": {"Fe2O3": 1},
    "ilmenite": {"FeO": 1, "TiO2": 1},
    "corundum": {"Al2O3": 1},
    "almandine": {"FeO": 3, "Al2O3": 1, "SiO2": 3},
    "grossular": {"CaO": 3, "Al2O3": 1, "SiO2": 3},
    "pyrope": {"MgO": 3, "Al2O3": 1, "SiO2": 3},
    "water": {"H2O": 1},
}


def liquid_stoichiometry(end_member: str) -> np.ndarray:
    """Coefficients of the liquid components that make one mole of the end member."""
    ox = np.array([OXIDES_OF[end_member].get(o, 0.0) for o in MODEL_OXIDES])
    nu = _TO_COMPONENTS @ ox
    nu[np.abs(nu) < 1e-12] = 0.0
    return nu


def mu_in_liquid(end_member: str, mu_liq: np.ndarray) -> float:
    """sum_k nu_ik mu_liq_k, skipping components the end member does not use (their mu may be -inf)."""
    nu = liquid_stoichiometry(end_member)
    nz = nu != 0
    return float(nu[nz] @ mu_liq[nz])


def _xlogx(x):
    return sum(v * math.log(v) for v in x if v > 0)


# ---- olivine (Sack and Ghiorso 1989; Hirschmann 1991), Fe-Mg binary --------------------------
W_OL = 10150.0 + 5075.0 + 5075.0


def g_ex_olivine(x, t, p):
    fo, fa = x
    return W_OL * fo * fa


# ---- feldspar (Elkins and Grove 1990, Table 4), order: albite, anorthite, sanidine ------------
EG = dict(whabor=18810.0, wsabor=10.3, wvabor=0.4602, whorab=27320.0, wsorab=10.3, wvorab=0.3264,
          whaban=7924.0, whanab=0.0, whoran=40317.0, whanor=38974.0, wvanor=-0.1037,
          whabanor=12545.0, wvabanor=-1.095)


# ---- garnet (Berman 1990; Berman and Koziol 1991, as in MELTS GarnetBerman.m) ----------------
# asymmetric ternary Margules on the three dodecahedral sites, W = WH - T WS + P WV per formula
# key ijk: 1 = grossular, 2 = pyrope, 3 = almandine
W_GT = {"112": (21560.0, 18.79, 0.10), "122": (69200.0, 18.79, 0.10), "113": (20320.0, 5.08, 0.17),
        "133": (2620.0, 5.08, 0.09), "223": (230.0, 0.0, 0.01), "233": (3720.0, 0.0, 0.06)}


def g_ex_garnet(x, t, p):
    al, gr, py = x
    w = {k: h - t * s + p * v for k, (h, s, v) in W_GT.items()}
    return (w["112"] * gr * py * (gr + al / 2) + w["122"] * gr * py * (py + al / 2)
            + w["113"] * gr * al * (gr + py / 2) + w["133"] * gr * al * (al + py / 2)
            + w["223"] * py * al * (py + gr / 2) + w["233"] * py * al * (al + gr / 2))


def g_ex_feldspar(x, t, p):
    ab, an, orr = x
    w = EG
    h = (w["whaban"] * ab * an * (an + orr / 2) + w["whanab"] * ab * an * (ab + orr / 2)
         + w["whabor"] * ab * orr * (orr + an / 2) + w["whorab"] * ab * orr * (ab + an / 2)
         + w["whanor"] * an * orr * (orr + ab / 2) + w["whoran"] * an * orr * (an + ab / 2)
         + w["whabanor"] * ab * an * orr)
    s = (w["wsabor"] * ab * orr * (orr + an / 2)
         + w["wsorab"] * ab * orr * (ab + an / 2))
    v = (w["wvabor"] * ab * orr * (orr + an / 2) + w["wvorab"] * ab * orr * (ab + an / 2)
         + w["wvanor"] * an * orr * (orr + ab / 2) + w["wvabanor"] * ab * an * orr)
    return h - t * s + (p - 1) * v


@dataclass(frozen=True)
class SolutionPhase:
    """G_mix = sites * RT sum x ln x + g_ex(x): ideal mixing of whole end members on `sites` sites."""
    name: str
    end_members: tuple[str, ...]
    g_ex: Callable
    sites: float = 1.0

    def g0(self, t, p) -> np.ndarray:
        return np.array([MINERALS[e].gibbs(t, p) for e in self.end_members])

    def g_mix(self, x, t, p) -> float:
        return self.sites * R * t * _xlogx(x) + self.g_ex(x, t, p)

    def mu(self, x, t, p, h=1e-6) -> np.ndarray:
        """Chemical potentials of the end members (Darken equation on g_ex, ideal part exact)."""
        x = np.asarray(x, float)
        gex = self.g_ex(x, t, p)
        grad = np.empty(len(x))
        for i in range(len(x)):
            up, dn = x.copy(), x.copy()
            up[i] += h
            dn[i] -= h
            grad[i] = (self.g_ex(up, t, p) - self.g_ex(dn, t, p)) / (2 * h)
        ex = gex + grad - x @ grad
        with np.errstate(divide="ignore"):
            ideal = self.sites * R * t * np.log(x)
        return self.g0(t, p) + ideal + ex

    def gibbs(self, x, t, p) -> float:
        x = np.asarray(x, float)
        return float(x @ self.g0(t, p)) + self.g_mix(x, t, p)

    # amounts as solver parameters: log moles of each end member
    @property
    def n_params(self) -> int:
        return len(self.end_members)

    def n_from_params(self, prm) -> np.ndarray:
        return np.exp(prm)

    def params_from(self, x, n_total) -> np.ndarray:
        return np.log(np.maximum(np.asarray(x) * n_total, 1e-12))

    def x_of(self, n) -> np.ndarray:
        return n / n.sum()

    def formula(self, x) -> str:
        if self.name == "olivine":
            return f"Fo{100 * x[0]:.0f}"
        if self.name == "feldspar":
            return f"An{100 * x[1]:.0f} Ab{100 * x[0]:.0f} Or{100 * x[2]:.0f}"
        if self.name == "garnet":
            return f"Py{100 * x[2]:.0f} Alm{100 * x[0]:.0f} Grs{100 * x[1]:.0f}"
        return ""

    def affinity(self, mu_liq: np.ndarray, t: float, p: float) -> tuple[float, np.ndarray]:
        """(A, x): A > 0 means oversaturated; x is the composition at the minimum."""
        d = self.g0(t, p) - np.array([mu_in_liquid(e, mu_liq) for e in self.end_members])
        n = len(d)

        def f(z):
            e = np.exp(z - z.max())
            x = e / e.sum()
            return float(x @ d) + self.g_mix(np.clip(x, 1e-300, 1), t, p)

        starts = [np.zeros(n)] + [np.where(np.arange(n) == i, 6.0, 0.0) for i in range(n)]
        best = min((minimize(f, z0, method="BFGS") for z0 in starts), key=lambda r: r.fun)
        e = np.exp(best.x - best.x.max())
        return -best.fun, e / e.sum()


OLIVINE = SolutionPhase("olivine", ("forsterite", "fayalite"), g_ex_olivine, sites=2.0)
FELDSPAR = SolutionPhase("feldspar", ("albite", "anorthite", "sanidine"), g_ex_feldspar)
GARNET = SolutionPhase("garnet", ("almandine", "grossular", "pyrope"), g_ex_garnet, sites=3.0)
from .pyroxene import CLINOPYROXENE, ORTHOPYROXENE, PIGEONITE   # noqa: E402
from .spinel import SPINEL, TITANOMAGNETITE              # noqa: E402
from .rhombohedral import RHOMBOHEDRAL                   # noqa: E402
from .water_fluid import WATER_FLUID                     # noqa: E402

PHASES = {p.name: p for p in (OLIVINE, FELDSPAR, CLINOPYROXENE, PIGEONITE, ORTHOPYROXENE, SPINEL, TITANOMAGNETITE, RHOMBOHEDRAL, GARNET, WATER_FLUID)}
