"""Saturation of a liquid of fixed composition in one solid phase (no crystallization)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .composition import Composition
from .redox import at_fo2
from .thermo.liquid import COMPONENTS, mu_mixing, to_components
from .thermo.phases import PHASES
from .thermo.standard_state import liquid_gibbs


@dataclass(frozen=True)
class Saturation:
    phase: str
    t_c: float
    affinity: float                 # J; > 0 oversaturated
    x: dict[str, float]             # composition that would crystallize


def saturation(phase: str, comp: Composition, t_c: float, p_bar: float = 1.0, fo2="QFM") -> Saturation:
    t = t_c + 273.15
    melt = at_fo2(comp, t_c, p_bar, fo2) if fo2 is not None else comp
    mu = np.array([liquid_gibbs(c, t, p_bar) for c in COMPONENTS]) + mu_mixing(to_components(melt).x, t)
    ph = PHASES[phase]
    a, x = ph.affinity(mu, t, p_bar)
    return Saturation(phase, t_c, a, dict(zip(ph.end_members, x.tolist())))


def saturation_temperature(phase: str, comp: Composition, p_bar=1.0, fo2="QFM", lo=900.0, hi=1700.0):
    """Temperature (degC) where the liquid, unchanged, just saturates in the phase."""
    t = brentq(lambda tc: saturation(phase, comp, tc, p_bar, fo2).affinity, lo, hi, xtol=0.01)
    return saturation(phase, comp, t, p_bar, fo2)
