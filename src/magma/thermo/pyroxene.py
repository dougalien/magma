"""Clinopyroxene and orthopyroxene: the MELTS model of Sack and Ghiorso (1994c), written for Magma.

Seven end members (values in data/melts_minerals.tsv, from ThermoEngine CpxBerman.m):
    diopside CaMgSi2O6, clinoenstatite Mg2Si2O6, hedenbergite CaFeSi2O6,
    alumino-buffonite CaTi0.5Mg0.5AlSiO6, buffonite CaTi0.5Mg0.5Fe3+SiO6,
    essenite CaFe3+AlSiO6, jadeite NaAlSi2O6.

Composition is described by six independent variables r and two ordering variables s:
    r0 = X_hd, r1 = X_abf + X_jd/2, r2 = X_buf - X_jd/2, r3 = X_ess + X_jd/2, r4 = X_jd, r5 = X_en
    s0: Fe3+-Al exchange between M1 and the tetrahedral site; s1: Fe2+-Mg exchange between M1 and M2.
End-member mole fractions can be negative (an Al-rich, Ti-poor augite needs negative buffonite);
what must stay positive is every site occupancy.

The non-configurational Gibbs energy is a cubic Taylor expansion in (r, s) whose coefficients are
in data/melts_cpx_taylor.tsv (evaluated from the parameters of CpxBerman.m).  The configurational
entropy is that of the site occupancies, with the coupled Ti-Mg (M1) and Al-Fe3+ (T) terms of the
model.  At each T, P and r the ordering s minimizes G.  The Gibbs energy of mixing is relative to
the seven end members in their own reference states (essenite with its own Fe3+-Al ordering), and
the end members carry the Berman-form standard states of minerals.py.

Our own implementation of the published model's equations (Douglas's approval, 2026-11-01).
"""
from __future__ import annotations

import math
from functools import lru_cache
from importlib import resources

import numpy as np
from scipy.optimize import minimize

from .minerals import MINERALS

R = 8.3143
END_MEMBERS = ("diopside", "clinoenstatite", "hedenbergite", "alumino-buffonite", "buffonite",
               "essenite", "jadeite")
_EPS = 1e-300


class Taylor:
    """The cubic Taylor expansion of one pyroxene structure: G = sum (H - T S + (P-1) V) * monomial."""

    def __init__(self, filename: str, names=("X2", "X3", "X4", "X5", "X6", "X7", "S1", "S2")):
        text = resources.files("magma.data").joinpath(filename).read_text("utf-8")
        rows = [ln.split("\t") for ln in text.splitlines() if ln and not ln.startswith("#")][1:]
        idx = {nm: k for k, nm in enumerate(names)}
        nv = len(names)
        self.deg = np.zeros((len(rows), nv), int)
        self.hsv = np.zeros((len(rows), 3))
        for k, row in enumerate(rows):
            name, vals = row[0], [float(v) for v in row[1:]] + [0.0] * (4 - len(row))
            if name != "0":
                for i in range(0, len(name), 2):
                    self.deg[k, idx[name[i:i + 2]]] += 1
            self.hsv[k] = vals[:3]
        # derivative of every monomial with respect to each variable: factor deg, lowered powers
        self.ddeg = np.stack([np.maximum(self.deg - np.eye(nv, dtype=int)[j], 0) for j in range(nv)])
        self.dfac = self.deg.T.astype(float)
        self._cache: dict = {}

    def coef(self, t, p) -> np.ndarray:
        key = (t, p)
        c = self._cache.get(key)
        if c is None:
            if len(self._cache) > 512:
                self._cache.clear()
            c = self._cache[key] = self.hsv[:, 0] - t * self.hsv[:, 1] + (p - 1.0) * self.hsv[:, 2]
        return c

    def poly(self, z, t, p) -> float:
        mono = np.prod(np.power(np.asarray(z, float)[None, :], self.deg), axis=1)
        return float(mono @ self.coef(t, p))

    def grad(self, z, t, p) -> np.ndarray:
        z = np.asarray(z, float)
        mono = np.prod(np.power(z[None, None, :], self.ddeg), axis=2) * self.dfac
        return mono @ self.coef(t, p)


CPX = Taylor("melts_cpx_taylor.tsv")
OPX = Taylor("melts_opx_taylor.tsv")


def poly(z, t, p, tay: Taylor = CPX) -> float:
    """Non-configurational G of the Taylor expansion at z = (r0..r5, s0, s1)."""
    return tay.poly(z, t, p)


# The configurational entropy as -R sum_k c_k v_k ln v_k with every v_k affine in z = (r, s).
#            r0    r1    r2    r3    r4     r5   s0    s1      b
_SITE = {
    "ti_m1":  ([0,   .5,   .5,   0,    0,     0,   0,    0],    0),
    "ca_m2":  ([0,   0,    0,    0,    -1,    -1,  0,    0],    1),
    "na_m2":  ([0,   0,    0,    0,    1,     0,   0,    0],    0),
    "si_t":   ([0,   -.5,  -.5,  -.5,  .25,   0,   0,    0],    1),
    "al_m1":  ([0,   0,    0,    .5,   .25,   0,   -.5,  0],    0),
    "fe3_m1": ([0,   0,    0,    .5,   .25,   0,   .5,   0],    0),
    "fe2_m1": ([1,   0,    0,    0,    0,     -.5, 0,    -.5],  0),
    "mg_m1":  ([-1,  -.5,  -.5,  -1,   -.5,   .5,  0,    .5],   1),
    "fe2_m2": ([0,   0,    0,    0,    0,     .5,  0,    .5],   0),
    "mg_m2":  ([0,   0,    0,    0,    0,     .5,  0,    -.5],  0),
    "al_t":   ([0,   .5,   0,    .25,  -.125, 0,   .25,  0],    0),
    "fe3_t":  ([0,   0,    .5,   .25,  -.125, 0,   -.25, 0],    0),
}


def _aff(name):
    a, b = _SITE[name]
    return np.array(a, float), float(b)


def _comb(*parts):
    a, b = np.zeros(8), 0.0
    for sign, name in parts:
        ai, bi = _aff(name) if isinstance(name, str) else name
        a, b = a + sign * ai, b + sign * bi
    return a, b


_ONE = (np.zeros(8), 1.0)
_MGFE = _comb((1, "mg_m1"), (1, "fe2_m1"))
_ENTROPY_TERMS = [            # (c_k, (a_k, b_k))
    (1, _aff("mg_m1")), (1, _aff("fe2_m1")), (1, _aff("al_m1")), (1, _aff("fe3_m1")), (1, _aff("ti_m1")),
    (1, _comb((1, _MGFE), (-1, "ti_m1"))), (-1, _comb((1, _ONE), (-1, "ti_m1"))), (-1, _MGFE),
    (-2, _comb((1, _ONE), (-1, "si_t"))), (2, _aff("al_t")), (2, _aff("fe3_t")),
    (1, _aff("ca_m2")), (1, _aff("na_m2")), (1, _aff("mg_m2")), (1, _aff("fe2_m2")),
    (1, _comb((1, _ONE), (-1, _MGFE), (-1, "na_m2"))), (-1, _comb((1, _ONE), (-1, _MGFE))),
    (-1, _comb((1, _ONE), (-1, "na_m2"))),
]
_EC = np.array([c for c, _ in _ENTROPY_TERMS], float)
_EA = np.array([a for _, (a, b) in _ENTROPY_TERMS])
_EB = np.array([b for _, (a, b) in _ENTROPY_TERMS])


def entropy_and_grad(z):
    v = _EA @ np.asarray(z, float) + _EB
    pos = v > 1e-300
    lv = np.where(pos, np.log(np.where(pos, v, 1.0)), 0.0)
    s = -R * float(_EC @ (v * lv))
    ds = -R * ((_EC * np.where(pos, lv + 1.0, 0.0)) @ _EA)
    return s, ds


# ---- composition -------------------------------------------------------------------------------
_RX = np.array([[0, 0, 0, 0, 0, 0],      # diopside
                [0, 0, 0, 0, 0, 1],      # clinoenstatite
                [1, 0, 0, 0, 0, 0],      # hedenbergite
                [0, 1, 0, 0, 0, 0],      # alumino-buffonite
                [0, 0, 1, 0, 0, 0],      # buffonite
                [0, 0, 0, 1, 0, 0],      # essenite
                [0, .5, -.5, .5, 1, 0]],  # jadeite
               float)


_XR = np.array([[-1, -1, -1, -1, -.5, -1], [0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, -.5, 0], [0, 0, 1, 0, .5, 0], [0, 0, 0, 1, -.5, 0], [0, 0, 0, 0, 1, 0]], float)


def r_from_x(x) -> np.ndarray:
    return np.asarray(x, float) @ _RX


def x_from_r(r) -> np.ndarray:
    r0, r1, r2, r3, r4, r5 = r
    return np.array([1 - r0 - r1 - r2 - r3 - r4 / 2 - r5, r5, r0, r1 - r4 / 2, r2 + r4 / 2,
                     r3 - r4 / 2, r4])


def sites(r, s) -> dict:
    r0, r1, r2, r3, r4, r5 = r
    s0, s1 = s
    return dict(
        ti_m1=(r1 + r2) / 2, ca_m2=1 - r4 - r5, na_m2=r4, si_t=(4 - 2 * r1 - 2 * r2 - 2 * r3 + r4) / 4,
        al_m1=(2 * r3 + r4 - 2 * s0) / 4, fe3_m1=(2 * r3 + r4 + 2 * s0) / 4,
        fe2_m1=(2 * r0 - r5 - s1) / 2, mg_m1=1 - r0 - r3 - 0.5 * (r1 + r2 + r4 - r5 - s1),
        fe2_m2=(r5 + s1) / 2, mg_m2=(r5 - s1) / 2,
        al_t=(4 * r1 + 2 * r3 - r4 + 2 * s0) / 8, fe3_t=(4 * r2 + 2 * r3 - r4 - 2 * s0) / 8,
    )


def _xl(v):
    return v * math.log(v) if v > _EPS else 0.0


def config_entropy(r, s) -> float:
    """Configurational entropy, J/K per mole of pyroxene."""
    y = sites(r, s)
    mgfe = y["mg_m1"] + y["fe2_m1"]
    return -R * (_xl(y["mg_m1"]) + _xl(y["fe2_m1"]) + _xl(y["al_m1"]) + _xl(y["fe3_m1"]) + _xl(y["ti_m1"])
                 + _xl(mgfe - y["ti_m1"]) - _xl(1 - y["ti_m1"]) - _xl(mgfe)
                 - 2 * _xl(1 - y["si_t"]) + 2 * _xl(y["al_t"]) + 2 * _xl(y["fe3_t"])
                 + _xl(y["ca_m2"]) + _xl(y["na_m2"]) + _xl(y["mg_m2"]) + _xl(y["fe2_m2"])
                 + _xl(1 - mgfe - y["na_m2"]) - _xl(1 - mgfe) - _xl(1 - y["na_m2"]))


def s_bounds(r):
    """Ranges of s0 and s1 that keep every ordered site occupancy positive."""
    r0, r1, r2, r3, r4, r5 = r
    lo0 = max(-(2 * r3 + r4) / 2, -(4 * r1 + 2 * r3 - r4) / 2)
    hi0 = min((2 * r3 + r4) / 2, (4 * r2 + 2 * r3 - r4) / 2)
    lo1 = max(-r5, 2 * (r0 + r3 - 1) + r1 + r2 + r4 - r5)
    hi1 = min(r5, 2 * r0 - r5)
    return (lo0, hi0), (lo1, hi1)


def g_rs(r, s, t, p, tay: Taylor = CPX) -> float:
    z = np.concatenate([r, s])
    return tay.poly(z, t, p) - t * entropy_and_grad(z)[0]


def g_and_grad(z, t, p, tay: Taylor = CPX):
    """G(r, s) and its partial derivatives with respect to the 8 variables."""
    s, ds = entropy_and_grad(z)
    return tay.poly(z, t, p) - t * s, tay.grad(z, t, p) - t * ds


def order(r, t, p, tay: Taylor = CPX, tol=1e-9):
    """(s, G) with s minimizing G at fixed r."""
    bounds = []
    for lo, hi in s_bounds(r):
        if hi - lo > 1e-12:
            pad = 1e-10 * (hi - lo)
            bounds.append((lo + pad, hi - pad))
        else:
            bounds.append((0.5 * (lo + hi),) * 2)
    s_start = np.array([0.5 * (b[0] + b[1]) for b in bounds])
    r = np.asarray(r, float)

    def f(s):
        g, dg = g_and_grad(np.concatenate([r, s]), t, p, tay)
        return g, dg[6:]

    res = minimize(f, s_start, jac=True, method="L-BFGS-B", bounds=bounds,
                   options=dict(ftol=1e-15, gtol=tol))
    return res.x, float(res.fun)


@lru_cache(maxsize=512)
def end_member_g(t: float, p: float) -> np.ndarray:
    """Taylor-expansion G of the seven end members in their reference states."""
    z = lambda r, s: np.array(list(r) + list(s), float)
    di = poly(z((0, 0, 0, 0, 0, 0), (0, 0)), t, p)
    en = poly(z((0, 0, 0, 0, 0, 1), (0, -1)), t, p)          # Mg on both M sites
    hd = poly(z((1, 0, 0, 0, 0, 0), (0, 0)), t, p)
    ca = poly(z((0, 1, 0, 0, 0, 0), (0, 0)), t, p)
    cf = poly(z((0, 0, 1, 0, 0, 0), (0, 0)), t, p)
    jd = poly(z((0, 0.5, -0.5, 0.5, 1, 0), (-1, 0)), t, p)   # Al on M1

    def es_g(v):
        s = float(np.clip(v[0], -1 + 1e-15, 1 - 1e-15))
        conf = R * t * ((1 - s) * math.log(1 - s) + (1 + s) * math.log(1 + s) - 2 * math.log(2))
        return conf + poly(z((0, 0, 0, 1, 0, 0), (s, 0)), t, p)

    es = minimize(es_g, [0.5], method="L-BFGS-B", bounds=[(-1 + 1e-12, 1 - 1e-12)],
                  options=dict(ftol=1e-15)).fun
    return np.array([di, en, hd, ca, cf, float(es), jd])


def _site_margins(z) -> np.ndarray:
    """Every site occupancy, and the coupled-site combinations in the entropy, less a margin."""
    r, s = z[:6], z[6:]
    y = sites(r, s)
    mgfe = y["mg_m1"] + y["fe2_m1"]
    v = [y[k] for k in ("al_m1", "fe3_m1", "fe2_m1", "mg_m1", "fe2_m2", "mg_m2", "al_t", "fe3_t",
                        "ca_m2", "na_m2", "ti_m1")]
    v += [1 - y["ti_m1"], 1 - y["si_t"], y["si_t"], mgfe - y["ti_m1"], 1 - mgfe - y["na_m2"]]
    return np.array(v) - 1e-10


# starting compositions for the affinity search (end-member mole fractions)
_AUGITE = [np.array(v, float) for v in (
    [0.55, 0.15, 0.15, 0.04, 0.02, 0.06, 0.03],     # augite
    [0.20, 0.60, 0.12, 0.02, 0.01, 0.03, 0.02],     # pigeonite-like
    [0.75, 0.05, 0.05, 0.05, 0.02, 0.05, 0.03],     # diopside
)]
_OPX = [np.array(v, float) for v in (
    [0.06, 0.74, 0.14, 0.02, 0.0, 0.03, 0.01],      # bronzite
    [0.10, 0.55, 0.30, 0.02, 0.0, 0.02, 0.01],      # hypersthene
    [0.20, 0.60, 0.12, 0.02, 0.01, 0.03, 0.02],     # pigeonite-like
)]


class Pyroxene:
    """Clinopyroxene or orthopyroxene: the same composition space and end-member reference
    states, different Taylor coefficients."""
    end_members = END_MEMBERS
    n_params = 7               # log N and r

    def __init__(self, name: str, tay: Taylor, starts):
        self.name, self.tay, self.starts = name, tay, starts

    def g0(self, t, p) -> np.ndarray:
        return np.array([MINERALS[e].gibbs(t, p) for e in END_MEMBERS])

    @staticmethod
    def feasible_x(x) -> bool:
        r = r_from_x(x)
        (lo0, hi0), (lo1, hi1) = s_bounds(r)
        if hi0 < lo0 or hi1 < lo1:
            return False
        return bool(_site_margins(np.concatenate([r, [0.5 * (lo0 + hi0), 0.5 * (lo1 + hi1)]])).min() >= -1e-10)

    def g_mix(self, x, t, p) -> float:
        _, g = order(r_from_x(x), t, p, self.tay)
        return g - float(np.asarray(x) @ end_member_g(t, p))

    def gibbs(self, x, t, p) -> float:
        x = np.asarray(x, float)
        return float(x @ self.g0(t, p)) + self.g_mix(x, t, p)

    def mu(self, x, t, p, h=1e-6) -> np.ndarray:
        """Chemical potentials of the end members.  With g(r) the molar G at the equilibrium
        ordering and r linear in x, mu_j = g0_j - E_j + g + sum_k dg/dr_k (r_k(e_j) - r_k)."""
        r = r_from_x(x)
        s, _ = order(r, t, p, self.tay)
        g, dg = g_and_grad(np.concatenate([r, s]), t, p, self.tay)
        grad = dg[:6]                    # at the equilibrium ordering, dG/dr = partial dG/dr
        return self.g0(t, p) - end_member_g(t, p) + g + (_RX - r) @ grad

    # ---- amounts as solver parameters ------------------------------------------------------
    def n_from_params(self, prm) -> np.ndarray:
        return math.exp(prm[0]) * x_from_r(prm[1:])

    def params_from(self, x, n_total) -> np.ndarray:
        return np.concatenate([[math.log(n_total)], r_from_x(x)])

    def x_of(self, n) -> np.ndarray:
        return n / n.sum()

    def affinity(self, mu_liq, t, p):
        """(A, x): A > 0 means oversaturated; x is the composition at the minimum."""
        from .phases import mu_in_liquid
        d = self.g0(t, p) - end_member_g(t, p) - np.array([mu_in_liquid(e, mu_liq) for e in END_MEMBERS])

        dz = np.concatenate([d @ _XR, [0.0, 0.0]])          # d(x . d)/dz

        def f(z):
            g, dg = g_and_grad(z, t, p, self.tay)
            return float(x_from_r(z[:6]) @ d) + g, dg + dz

        cons = [{"type": "ineq", "fun": _site_margins}]
        best = None
        for x0 in self.starts:
            r0 = r_from_x(x0)
            z0 = np.concatenate([r0, [0.5 * (lo + hi) for lo, hi in s_bounds(r0)]])
            res = minimize(f, z0, jac=True, method="SLSQP", constraints=cons,
                           options=dict(maxiter=500, ftol=1e-12))
            if (np.all(np.isfinite(res.x)) and _site_margins(res.x).min() > -1e-8
                    and (best is None or res.fun < best.fun)):
                best = res
        if best is None:
            return -np.inf, self.starts[0]
        r = best.x[:6]
        return -float(f(np.concatenate([r, order(r, t, p, self.tay)[0]]))[0]), x_from_r(r)

    @staticmethod
    def formula(x) -> str:
        r = r_from_x(x)
        ca = 1 - r[4] - r[5]
        fe = r[0]
        mg = 1 - r[0] - r[3] - 0.5 * (r[1] + r[2] + r[4]) + r[5]          # Mg on M1 + M2
        tot = ca + mg + fe
        return f"Wo{100 * ca / tot:.0f} En{100 * mg / tot:.0f} Fs{100 * fe / tot:.0f}"


CLINOPYROXENE = Pyroxene("clinopyroxene", CPX, _AUGITE)
ORTHOPYROXENE = Pyroxene("orthopyroxene", OPX, _OPX)
