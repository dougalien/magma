"""End-member Gibbs energies in the Holland & Powell (2011) dataset (data/hp_ds636.tsv).

Units follow the dataset: kJ, kbar, K.  `gibbs_j` gives J/mol at T in K and P in bar for use beside the
MELTS engine.  G is the apparent Gibbs energy of formation from the elements (H at 298.15 K, 1 bar, minus T S):

    G(T, P) = H - T S + int Cp dT - T int Cp/T dT + int V dP + G_order

* Heat capacity Cp = a + b T + c / T^2 + d / sqrt(T).
* Solids: the modified Tait equation of state with a thermal pressure from one Einstein temperature,
  theta = 10636 / (S / atoms + 6.44) (S in J/K), rounded to a whole kelvin as THERMOCALC does
  (Holland and Powell 2011, eqs 3-4 and the appendix).
* Liquid end members (names ending in L): V(T) = V0 exp(alpha0 (T - 298.15)) and
  kappa(T) = kappa0 + dkappa/dT (T - 298.15) in the same Tait form, with no thermal pressure
  (Holland and Powell 2011, section on melts; Holland, Green and Powell 2018).
* Order-disorder: the Landau model of Holland and Powell (1998) for quartz, hematite, magnetite and
  others; Bragg-Williams for feldspars, spinels, ilmenite, sillimanite and others (Holland and Powell 1996).
  Both are taken at internal equilibrium unless a state is given.
* H2O and CO2 fluids: the Pitzer and Sterner (1994) equation of state.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

from scipy.optimize import brentq

T0, P0 = 298.15, 0.001        # K, kbar (1 bar)
R = 0.0083144                 # kJ/K/mol
EPS = 1e-8

OXIDES = ("SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "MnO", "Cr2O3", "H2O", "CO2",
          "Sel", "Cl")


@dataclass(frozen=True)
class EndMember:
    name: str
    comp: dict
    atoms: float
    h: float
    s: float
    v: float
    cp: tuple
    alpha0: float
    kappa0: float
    kappa0p: float
    kappa0pp: float
    od: float
    odp: tuple

    @property
    def is_liquid(self) -> bool:
        return self.name.endswith("L")


def _load(file: str = "hp_ds636.tsv") -> dict[str, EndMember]:
    text = resources.files("magma.data").joinpath(file).read_text("utf-8")
    rows = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    out = {}
    for r in csv.DictReader(rows, delimiter="\t"):
        f = lambda k: float(r[k])
        comp = {o: f(o) for o in OXIDES if f(o) != 0.0}
        out[r["name"]] = EndMember(
            r["name"], comp, f("atoms"), f("H"), f("S"), f("V"), (f("a"), f("b"), f("c"), f("d")),
            f("alpha0"), f("kappa0"), f("kappa0p"), f("kappa0pp"), f("od"),
            tuple(f(f"od{i}") for i in range(1, 7)))
    return out


DATA = _load()
_AQUEOUS = {"H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3--", "AlOH3", "AlOH4-", "KOH",
            "HCl", "KCl", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl+", "FeCl2", "aqSi", "HS-", "HSO3-", "SO42-",
            "HSO4-"}


def _cp_terms(cp, t: float) -> float:
    """int_T0^T Cp dT - T int_T0^T Cp/T dT, kJ."""
    a, b, c, d = cp
    h = a * (t - T0) + b * (t * t - T0 * T0) / 2 - c * (1 / t - 1 / T0) + 2 * d * (math.sqrt(t) - math.sqrt(T0))
    s = (a * math.log(t / T0) + b * (t - T0) - c / 2 * (t ** -2 - T0 ** -2)
         - 2 * d * (t ** -0.5 - T0 ** -0.5))
    return h - t * s


def _tait(v0: float, kappa0: float, kp: float, kpp: float, p: float, pth: float) -> float:
    """int_P0^P V dP of the modified Tait equation (Holland and Powell 2011, eq 3), kJ."""
    a = (1 + kp) / (1 + kp + kappa0 * kpp)
    b = (kp + kp * kp - kappa0 * kpp) / (kappa0 * (1 + kp))
    c = (1 + kp + kappa0 * kpp) / (kp + kp * kp - kappa0 * kpp)
    num = (p - P0) * (1 - a) + a * (-(1 + b * (p - pth)) ** (1 - c) + (1 + b * (P0 - pth)) ** (1 - c)) / (b * (c - 1))
    return v0 * num / ((1 - a) + a * (1 + b * P0) ** (-c))


def _thermal_pressure(e: EndMember, t: float) -> float:
    theta = math.floor(10636.0 / (e.s * 1e3 / e.atoms + 6.44) + 0.5)
    u0 = theta / T0
    xi0 = math.exp(u0) * u0 * u0 / (math.exp(u0) - 1) ** 2
    return theta * e.alpha0 * e.kappa0 / xi0 * (1 / (math.exp(theta / t) - 1) - 1 / (math.exp(u0) - 1))


# ---- fluids: Pitzer and Sterner (1994) -----------------------------------------------------------------
def _ps_coeffs(name: str, t: float):
    if name == "H2O":
        return (0.24657688e6 / t + 0.51359951e2,
                0.58638965e0 / t - 0.28646939e-2 + 0.31375577e-4 * t,
                -0.62783840e1 / t + 0.14791599e-1 + 0.35779579e-3 * t + 0.15432925e-7 * t * t,
                -0.42719875e0 - 0.16325155e-4 * t,
                0.56654978e4 / t - 0.16580167e2 + 0.76560762e-1 * t,
                0.10917883e0,
                0.38878656e13 / t ** 4 - 0.13494878e9 / t ** 2 + 0.30916564e6 / t + 0.75591105e1,
                -0.65537898e5 / t + 0.18810675e3,
                -0.14182435e14 / t ** 4 + 0.18165390e9 / t ** 2 - 0.19769068e6 / t - 0.23530318e2,
                0.92093375e5 / t + 0.12246777e3)
    return (0.18261340e7 / t + 0.79224365e2,
            0.66560660e-4 + 0.57152798e-5 * t + 0.30222363e-9 * t * t,
            0.59957845e-2 + 0.71669631e-4 * t + 0.62416103e-8 * t * t,
            -0.13270279e1 / t - 0.15210731e0 + 0.53654244e-3 * t - 0.71115142e-7 * t * t,
            0.12456776e0 / t + 0.49045367e1 + 0.98220560e-2 * t + 0.55962121e-5 * t * t,
            0.75522299e0,
            -0.39344644e12 / t ** 4 + 0.90918237e8 / t ** 2 + 0.42776716e6 / t - 0.22347856e2,
            0.40282608e3 / t + 0.11971627e3,
            0.22995650e8 / t ** 2 - 0.78971817e5 / t - 0.63376456e2,
            0.95029765e5 / t + 0.18038071e2)


R_CM3 = 83.144     # cm3 bar / K / mol


def pitzer_sterner_volume(name: str, t: float, p_bar: float) -> float:
    """Molar volume (cm3/mol) of H2O or CO2 at T (K) and P (bar), Pitzer and Sterner (1994) eq 2."""
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 = _ps_coeffs(name, t)

    def pressure(v):
        r = 1 / v
        den = c2 + c3 * r + c4 * r * r + c5 * r ** 3 + c6 * r ** 4
        dden = c3 + 2 * c4 * r + 3 * c5 * r * r + 4 * c6 * r ** 3
        return R_CM3 * t * (r + c1 * r * r - r * r * dden / den ** 2 + c7 * r * r * math.exp(-c8 * r)
                            + c9 * r * r * math.exp(-c10 * r))

    return brentq(lambda v: pressure(v) - p_bar, 3.0, R_CM3 * t / (p_bar / 1000.0), xtol=1e-14, maxiter=500)


def _fluid_term(name: str, t: float, p: float) -> float:
    """Pitzer-Sterner G of the fluid relative to the ideal gas at 1 bar, kJ (P in kbar)."""
    p_bar = 1000.0 * p
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 = _ps_coeffs(name, t)
    v = pitzer_sterner_volume(name, t, p_bar)
    r = 1 / v
    a_res = R_CM3 * t * (c1 * r + (1 / (c2 + c3 * r + c4 * r * r + c5 * r ** 3 + c6 * r ** 4) - 1 / c2)
                         - c7 / c8 * (math.exp(-c8 * r) - 1) - c9 / c10 * (math.exp(-c10 * r) - 1))
    return (a_res + p_bar * v + R_CM3 * t * (math.log(R_CM3 * t / v) - 1)) * 1e-4


# ---- order-disorder ------------------------------------------------------------------------------------
def _landau(e: EndMember, t: float, p: float, state: str) -> float:
    tc0, smax, vmax = e.odp[0], e.odp[1], e.odp[2]
    q20 = math.sqrt(1 - T0 / tc0)
    if state == "ordered":
        return smax * tc0 * (-2 / 3 + q20 * (1 - q20 * q20 / 3)) - t * smax * (q20 - 1) + p * vmax * (q20 - 1)
    if state == "disordered":
        return smax * tc0 * q20 * (1 - q20 * q20 / 3) - t * smax * q20 + p * vmax * q20
    tc = tc0 + (p * vmax / smax if vmax else 0.0)
    q2 = 0.0 if t > tc else math.sqrt((tc - t) / tc0)
    return (smax * (tc0 * (q20 * (1 - q20 * q20 / 3) + q2 ** 3 / 3) - q2 * tc) - t * smax * (q20 - q2)
            + p * vmax * q20)


def _bragg_williams(e: EndMember, t: float, p: float, state: str) -> float:
    dh, dhv, w, wv, n, fac = e.odp
    if n <= 0:
        return 0.0
    if state == "ordered":
        return 0.0
    if state == "disordered":
        sod = fac * R * (math.log(1 / (n + 1)) + n * math.log(n / (n + 1)))
        if fac < 0:
            sod *= (1 / fac - n) / (n + 1)
        return dh + p * dhv + t * sod

    if fac < 0:
        def dg(q):
            return (dh + p * dhv + (w + p * wv) * (2 * q - 1)
                    + n / (n + 1) * R * t * (math.log(n * (1 - q) / (1 + n * q)) - fac * math.log((1 - q) / (n + q))))
        lo, hi = EPS, 1 - EPS
    else:
        def dg(q):
            return (dh + p * dhv + (w + p * wv) * (2 * q - 1)
                    + fac * n / (n + 1) * R * t * math.log(n * (1 - q) ** 2 / ((1 + n * q) * (n + q))))
        lo, hi = (EPS, 1 - EPS) if dg(EPS) * dg(1 - EPS) <= 0 else (0.0, 1 - EPS)
    f_lo, f_hi = dg(lo), dg(hi)
    if f_lo * f_hi > 0:
        q = lo if abs(f_lo) < abs(f_hi) else hi
    else:
        q = brentq(dg, lo, hi, xtol=1e-12, maxiter=500)

    if fac < 0:
        sod = ((1 + n * q) * math.log((1 + n * q) / (n + 1)) + n * (1 - q) * math.log(n * (1 - q) / (n + 1))
               - fac * (n * (1 - q) * math.log((1 - q) / (n + 1)) + n * (n + q) * math.log((n + q) / (n + 1)))) / (n + 1)
    else:
        sod = fac * ((1 + n * q) * math.log((1 + n * q) / (n + 1)) + n * (1 - q) * math.log((1 - q) / (n + 1))
                     + n * (1 - q) * math.log(n * (1 - q) / (n + 1)) + n * (n + q) * math.log((n + q) / (n + 1))) / (n + 1)
    return dh + p * dhv + q * (w - dh + p * (wv - dhv)) - q * q * (w + p * wv) + R * t * sod


# ---- G ---------------------------------------------------------------------------------------------------
@lru_cache(maxsize=65536)
def gibbs(name: str, t: float, p: float, state: str = "equilibrium") -> float:
    """G of an end member, kJ/mol, at T (K) and P (kbar)."""
    if name in _AQUEOUS:
        raise ValueError(f"{name} is an aqueous species (DEW model), not part of this engine")
    e = DATA[name]
    g = e.h - t * e.s + _cp_terms(e.cp, t)

    if name in ("H2O", "CO2"):
        return g + _fluid_term(name, t, p)
    if name == "O2":
        return g
    if e.is_liquid:
        v0 = e.v * math.exp(e.alpha0 * (t - T0))
        kappa0 = e.kappa0 + e.od * (t - T0)
        return g + _tait(v0, kappa0, e.kappa0p, e.kappa0pp, p, 0.0)

    g += _tait(e.v, e.kappa0, e.kappa0p, e.kappa0pp, p, _thermal_pressure(e, t))
    if e.od == 1.0 and e.odp[1] > 0:
        g += _landau(e, t, p, state)
    elif e.od == 2.0:
        g += _bragg_williams(e, t, p, state)
    return g


def gibbs_j(name: str, t_k: float, p_bar: float, state: str = "equilibrium") -> float:
    """G in J/mol at T (K) and P (bar)."""
    return 1000.0 * gibbs(name, float(t_k), float(p_bar) / 1000.0, state)
