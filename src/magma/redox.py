"""Oxygen fugacity: solid buffers and the ferric/ferrous ratio of silicate melts.

Buffers: Frost (1991), Reviews in Mineralogy 25, Table 1,
    log10 fO2 = A/T + B + C (P - 1)/T          T in K, P in bar

Melt Fe3+/Fe2+: Kress and Carmichael (1991), Contrib. Mineral. Petrol. 108, 82-92, eq. 7,
    ln(X_Fe2O3 / X_FeO) = a ln fO2 + b/T + c + sum_i d_i X_i
                          + e [1 - T0/T - ln(T/T0)] + f P/T + g (T - T0) P/T + h P^2/T
with anhydrous oxide mole fractions X_i, total iron counted as FeO in the sum, P in Pa.
"""
from __future__ import annotations

import math

from .composition import Composition

# ---- buffers ------------------------------------------------------------------------------------
# name -> list of (T_min_K, A, B, C); the last row whose T_min <= T applies
_BUFFERS: dict[str, list[tuple[float, float, float, float]]] = {
    # quartz-fayalite-magnetite: alpha-quartz below 573 C, beta-quartz above
    "QFM": [(0.0, -26455.3, 10.344, 0.092), (846.15, -25096.3, 8.735, 0.110)],
    "NNO": [(0.0, -24930.0, 9.36, 0.046)],
    "IW": [(0.0, -27489.0, 6.702, 0.055)],
    "MH": [(0.0, -25700.6, 14.558, 0.019), (950.0, -26452.6, 15.455, 0.019)],
}
BUFFERS = tuple(_BUFFERS)


def buffer_log_fo2(buffer: str, t_c: float, p_bar: float = 1.0) -> float:
    """log10 fO2 of a buffer at T (degC) and P (bar)."""
    key = buffer.upper().replace("FMQ", "QFM")
    if key not in _BUFFERS:
        raise ValueError(f"unknown buffer {buffer!r}; known: {', '.join(BUFFERS)}")
    t = t_c + 273.15
    a, b, c = next((a, b, c) for tmin, a, b, c in reversed(_BUFFERS[key]) if t >= tmin)
    return a / t + b + c * (p_bar - 1.0) / t


def log_fo2(spec: str | float, t_c: float, p_bar: float = 1.0) -> float:
    """log10 fO2 from a number or a buffer with offset: 'QFM', 'QFM+1', 'NNO-0.5'."""
    if isinstance(spec, (int, float)):
        return float(spec)
    s = spec.strip().replace(" ", "").upper()
    for sign in ("+", "-"):
        i = s.find(sign, 1)
        if i > 0:
            return buffer_log_fo2(s[:i], t_c, p_bar) + float(s[i:])
    return buffer_log_fo2(s, t_c, p_bar)


# ---- Kress and Carmichael (1991) ----------------------------------------------------------------
KC91 = {
    "a": 0.196, "b": 1.1492e4, "c": -6.675,
    "d": {"Al2O3": -2.243, "FeO": -1.828, "CaO": 3.201, "Na2O": 5.854, "K2O": 6.215},
    "e": -3.36, "f": -7.01e-7, "g": -1.54e-10, "h": 3.85e-17, "T0": 1673.0,
}


def _kc91_rest(comp: Composition, t_k: float, p_bar: float) -> float:
    """Everything in eq. 7 except a ln fO2."""
    k = KC91
    x = comp.mole_fractions(volatile_free=True, total_iron_as_feo=True)
    p = p_bar * 1e5
    t0 = k["T0"]
    return (k["b"] / t_k + k["c"] + sum(d * x.get(ox, 0.0) for ox, d in k["d"].items())
            + k["e"] * (1.0 - t0 / t_k - math.log(t_k / t0))
            + k["f"] * p / t_k + k["g"] * (t_k - t0) * p / t_k + k["h"] * p * p / t_k)


def fe3_ratio(comp: Composition, t_c: float, p_bar: float, fo2: str | float) -> float:
    """Fe3+ / total Fe of a melt at T (degC), P (bar) and fO2 (log10 value or buffer string)."""
    t = t_c + 273.15
    lf = log_fo2(fo2, t_c, p_bar) * math.log(10.0)
    r = math.exp(KC91["a"] * lf + _kc91_rest(comp, t, p_bar))      # X_Fe2O3 / X_FeO
    return 2.0 * r / (1.0 + 2.0 * r)


def melt_log_fo2(comp: Composition, t_c: float, p_bar: float) -> float:
    """log10 fO2 implied by the melt's own Fe2O3 and FeO (inverse of eq. 7)."""
    n = comp.moles()
    if n.get("Fe2O3", 0.0) <= 0 or n.get("FeO", 0.0) <= 0:
        raise ValueError("need both FeO and Fe2O3 to compute fO2 from the melt")
    r = n["Fe2O3"] / n["FeO"]
    t = t_c + 273.15
    return (math.log(r) - _kc91_rest(comp, t, p_bar)) / KC91["a"] / math.log(10.0)


def at_fo2(comp: Composition, t_c: float, p_bar: float, fo2: str | float) -> Composition:
    """The composition with its total iron split for the given fO2."""
    return comp.with_fe3_ratio(fe3_ratio(comp, t_c, p_bar, fo2))
