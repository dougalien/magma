"""Mineral end members (Berman 1988 form) with the values MELTS uses (data/melts_minerals.tsv).

Gibbs energies are Berman-Brown apparent energies, H - T S with absolute S, J/mol, the same
convention as the liquid standard states.

Two special cases:
* albite: plus the Salje (1985) Landau term for Al-Si (q) and monoclinic-triclinic (q_od) ordering,
  minimized over q and q_od at each T, P.  Above about 1270 K both vanish.
* sanidine: plus the Berman (1988) disorder function, integrated from 298.15 K to min(T, 1436 K).
* lambda:l1,l2,Tt,dH (magnetite): Berman (1988) lambda heat capacity T (l1 + l2 T)^2 up to Tt,
  plus the enthalpy dH at Tt.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

import numpy as np
from scipy.optimize import minimize

TR, PR = 298.15, 1.0
_SANIDINE_DIS = (282.98, -4830.0, 3621000.0, -0.15733, 34.770e-6, 410630.0)
_SALJE = dict(a0=5.479, b=6854.0, aod0=41.620, bod=-9301.0, cod=43600.0, d0=-2.171, d1=-3.043,
              d2=-0.001569, d3=0.000002109, tc=1251.0, tod=824.1, vscale=335282.925)


@dataclass(frozen=True)
class EndMember:
    name: str
    formula: str
    h: float
    s: float
    k: tuple[float, float, float, float]
    v: tuple[float, float, float, float, float]   # v0, v1 (P), v2 (P^2), v3 (T), v4 (T^2)
    special: str = ""

    def gibbs(self, t: float, p: float) -> float:
        k0, k1, k2, k3 = self.k
        h = (self.h + k0 * (t - TR) + 2 * k1 * (math.sqrt(t) - math.sqrt(TR)) - k2 * (1 / t - 1 / TR)
             - 0.5 * k3 * (1 / t ** 2 - 1 / TR ** 2))
        s = (self.s + k0 * math.log(t / TR) - 2 * k1 * (1 / math.sqrt(t) - 1 / math.sqrt(TR))
             - 0.5 * k2 * (1 / t ** 2 - 1 / TR ** 2) - k3 / 3 * (1 / t ** 3 - 1 / TR ** 3))
        v0, v1, v2, v3, v4 = self.v
        dt, dp = t - TR, p - PR
        g = h - t * s + v0 * ((1 + v3 * dt + v4 * dt * dt) * dp + v1 / 2 * dp * dp + v2 / 3 * dp ** 3)
        if self.special.startswith("lambda:"):
            l1, l2, tt, dh = (float(v) for v in self.special[7:].split(","))
            tl = min(t, tt)
            h_l = 0.5 * l1 * l1 * (tl ** 2 - TR ** 2) + 2 / 3 * l1 * l2 * (tl ** 3 - TR ** 3) \
                + 0.25 * l2 * l2 * (tl ** 4 - TR ** 4)
            s_l = l1 * l1 * (tl - TR) + l1 * l2 * (tl ** 2 - TR ** 2) + l2 * l2 / 3 * (tl ** 3 - TR ** 3)
            if t > tt:
                h_l += dh
                s_l += dh / tt
            g += h_l - t * s_l
        if self.special == "disorder":
            g += _disorder_g(t, p)
        elif self.special == "albite_salje":
            g += _salje_g(t, p)
        return g


def _disorder_g(t: float, p: float) -> float:
    d0, d1, d2, d3, d4, d5 = _SANIDINE_DIS
    td = min(t, 1436.0)
    dh = (d0 * (td - TR) + 2 * d1 * (math.sqrt(td) - math.sqrt(TR)) - d2 * (1 / td - 1 / TR)
          + d3 * (td * td - TR * TR) / 2 + d4 * (td ** 3 - TR ** 3) / 3)
    ds = (d0 * math.log(td / TR) - 2 * d1 * (1 / math.sqrt(td) - 1 / math.sqrt(TR))
          - d2 * (1 / td ** 2 - 1 / TR ** 2) / 2 + d3 * (td - TR) + d4 * (td * td - TR * TR) / 2)
    return dh - t * ds + (dh / d5) * (p - PR)


@lru_cache(maxsize=4096)
def _salje_g(t: float, p: float) -> float:
    c = _SALJE
    cross_h = c["d0"] - c["d2"] * t * t - 2 * c["d3"] * t ** 3
    cross_s = c["d1"] + 2 * c["d2"] * t + 3 * c["d3"] * t * t
    pf = 1 + (p - 1) / c["vscale"]

    def g(q):
        q0, q1 = q
        h = (-0.5 * c["a0"] * c["tc"] * q0 ** 2 + 0.25 * c["b"] * q0 ** 4
             - 0.5 * c["aod0"] * c["tod"] * q1 ** 2 + 0.25 * c["bod"] * q1 ** 4 + c["cod"] / 6 * q1 ** 6
             + cross_h * q0 * q1)
        s = -(0.5 * c["a0"] * q0 ** 2 + 0.5 * c["aod0"] * q1 ** 2 + cross_s * q0 * q1)
        return h * pf - t * s

    best = min((minimize(g, x0, method="L-BFGS-B", bounds=[(-1, 1), (-1, 1)]) for x0 in
                ([0.0, 0.0], [0.9, 0.9], [0.5, 0.2])), key=lambda r: r.fun)
    return float(min(best.fun, 0.0))


def _load() -> dict[str, EndMember]:
    text = resources.files("magma.data").joinpath("melts_minerals.tsv").read_text("utf-8")
    rows = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    out = {}
    for r in csv.DictReader(rows, delimiter="\t"):
        f = lambda k: float(r[k])
        out[r["name"]] = EndMember(r["name"], r["formula"], f("H"), f("S"),
                                   (f("k0"), f("k1"), f("k2"), f("k3")),
                                   (f("v0"), f("v1"), f("v2"), f("v3"), f("v4")), r["special"] or "")
    return out


MINERALS = _load()
