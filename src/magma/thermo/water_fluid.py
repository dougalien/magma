"""A pure-water fluid phase (MELTS's "Water"), saturating when the melt holds as much H2O as it can.

G of the fluid is measured from pure water at the same T and 1 bar, like the dissolved-H2O standard
state of the liquid (standard_state.h2o_liquid_gibbs), so G_fluid = G_water(T, P) - G_water(T, 1 bar),
from IAPWS-95 (thermo/water.py).  MELTS uses the Haar et al. (1984) equation of state, and the
Delany and Helgeson (1978) correction above 10 kbar, for the same quantity.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from . import water


@lru_cache(maxsize=4096)
def fluid_gibbs(t_k: float, p_bar: float) -> float:
    """G of the water fluid on the minerals' scale (water_gibbs_1bar plus the pressure term)."""
    if p_bar <= 1.0:
        return water_gibbs_1bar(t_k)
    return water.gibbs_iapws(t_k, p_bar) - water.gibbs_iapws(t_k, 1.0) + water_gibbs_1bar(t_k)


# Water on the Berman-Brown apparent-energy scale of the minerals (H = enthalpy of formation from
# the elements at 298.15 K, S = third-law entropy), so that hydrous minerals (micas, amphibole) can be
# set against the H2O of the melt.  Liquid water at 298.15 K and 1 bar: H = -285830 J, S = 69.95 J/K
# (Robie and Hemingway 1995, the values of the Berman 1988 database).  IAPWS-95 has its own zero
# (liquid at the triple point), so the two differ by a + b T.
H_F_WATER, S_WATER = -285830.0, 69.95


@lru_cache(maxsize=1)
def _iapws_offset() -> tuple[float, float]:
    t, h = 298.15, 1e-3
    g = water.gibbs_iapws(t, 1.0)
    s = -(water.gibbs_iapws(t + h, 1.0) - water.gibbs_iapws(t - h, 1.0)) / (2 * h)
    return H_F_WATER - (g + t * s), S_WATER - s


@lru_cache(maxsize=4096)
def water_gibbs_1bar(t_k: float) -> float:
    """G of pure water (steam above 373 K) at T and 1 bar on the minerals' scale, J/mol."""
    dh, ds = _iapws_offset()
    return water.gibbs_iapws(t_k, 1.0) + dh - t_k * ds


class WaterFluid:
    name = "water"
    end_members = ("water",)
    n_params = 1

    def g0(self, t, p) -> np.ndarray:
        return np.array([fluid_gibbs(float(t), float(p))])

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
        from .phases import mu_in_liquid
        m = mu_in_liquid("water", mu_liq)
        if not np.isfinite(m):
            return -np.inf, np.array([1.0])
        return m - float(self.g0(t, p)[0]), np.array([1.0])

    @staticmethod
    def formula(x) -> str:
        return "H2O"


WATER_FLUID = WaterFluid()
