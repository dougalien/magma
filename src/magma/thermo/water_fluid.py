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
    if p_bar <= 1.0:
        return 0.0
    return water.gibbs_iapws(t_k, p_bar) - water.gibbs_iapws(t_k, 1.0)


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
