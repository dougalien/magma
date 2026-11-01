"""Rhombohedral oxide (ilmenite-hematite): the MELTS model of Ghiorso and Evans (2008), written for Magma.

End members (MELTS order): geikielite MgTiO3, hematite Fe2O3, ilmenite FeTiO3, pyrophanite MnTiO3,
corundum Al2O3.  Magma's liquid has no MnO, so pyrophanite is left out (its fraction is zero) and the
phase has four end members: geikielite, hematite, ilmenite, corundum.

Composition: r0 = X_ilmenite, r1 = X_geikielite, r2 = X_pyrophanite, r3 = X_corundum; hematite is
the rest, q = 1 - r0 - r1 - r2 - r3.  Ordering s_i of Fe2+ (i = 0), Mg (1), Mn (2) against Ti over
the two layers A and B (|s_i| <= r_i):

    A: M_i (r_i + s_i)/2, Ti sum(r_i - s_i)/2;   B: M_i (r_i - s_i)/2, Ti sum(r_i + s_i)/2
    disordered ("ID"): M_i r_i/2, Ti sum(r_i)/2, Al r3, Fe3+ q

S = -R [sum_A x ln x + sum_B x ln x - 2 sum_ID(M, Ti) x ln x]
    - (1 - c) 2R sum_ID(all) x ln x + c 2R ln 2,  with c = 0.0730205 (short-range order)
G = H(r, s, P) - T S, with s minimizing G from MELTS's ordered start (s = 0.9 r).  End-member
reference states are the pure phases with their own ordering.  Parameter values from ThermoEngine
RhombohedralBerman.m (published numbers, cited; Douglas's approval 2026-11-01); our own code.
The ternary (W_ilm-ilm-gei ...) terms of the model are all zero in MELTS and are left out here.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.optimize import minimize

from .minerals import MINERALS

R = 8.3143
C_SRO = 0.0730205
END_MEMBERS = ("geikielite", "hematite", "ilmenite", "corundum")
MELTS_ORDER = ("geikielite", "hematite", "ilmenite", "pyrophanite", "corundum")

# per ordering cation i = ilmenite (Fe), geikielite (Mg), pyrophanite (Mn)
DV, WV = 0.010758, 0.035089                    # ordering volume terms, J/bar
DWV_HM, WV_HM2, DWV_CRN, WV_HM = 0.013701, 0.0, 0.013701, -0.11764
DH, WH = 17477.0, 3189.0
DWH_HM, WH_HM2, DWH_CRN = -5626.63, -833.14, 0.0
WH_HM_CRN = 69000.0
WH_HM, WCRN = 22535.6, 22535.6
# cation pairs (0,1) ilm-gei, (0,2) ilm-pyr, (1,2) gei-pyr: W, W_V, W_T
PAIRS = {(0, 1): (2600.0, 0.0, 88099.9), (0, 2): (2200.0, 0.0, 30244.0), (1, 2): (2600.0, 0.0, 2600.0)}


def _site_matrix():
    """Rows: site fractions as A z + b for z = (r0, r1, r2, r3, s0, s1, s2); weights for sum w x ln x."""
    rows, b, w = [], [], []

    def add(coef, const, weight):
        rows.append(coef)
        b.append(const)
        w.append(weight)

    for sign in (+1, -1):                      # layer A, then B
        for i in range(3):
            c = np.zeros(7)
            c[i], c[4 + i] = 0.5, 0.5 * sign
            add(c, 0.0, 1.0)
        c = np.zeros(7)
        c[:3], c[4:] = 0.5, -0.5 * sign
        add(c, 0.0, 1.0)
    for i in range(3):                         # disordered M_i and Ti: -2 and (1 - c) 2
        c = np.zeros(7)
        c[i] = 0.5
        add(c, 0.0, -2.0 + 2.0 * (1 - C_SRO))
    c = np.zeros(7)
    c[:3] = 0.5
    add(c, 0.0, -2.0 + 2.0 * (1 - C_SRO))
    c = np.zeros(7)
    c[3] = 1.0
    add(c, 0.0, 2.0 * (1 - C_SRO))             # Al
    c = np.zeros(7)
    c[:4] = -1.0
    add(c, 1.0, 2.0 * (1 - C_SRO))             # Fe3+
    return np.array(rows), np.array(b), np.array(w)


_A, _B, _W = _site_matrix()


def site_fractions(z):
    return _A @ np.asarray(z, float) + _B


def entropy_and_grad(z):
    v = site_fractions(z)
    pos = v > 1e-300
    lv = np.where(pos, np.log(np.where(pos, v, 1.0)), 0.0)
    s = -R * float(_W @ (v * lv)) + C_SRO * 2 * R * math.log(2.0)
    return s, -R * ((_W * np.where(pos, lv + 1.0, 0.0)) @ _A)


def h_and_grad(z, p):
    """The enthalpy-like part H(r, s) of the model (it includes the P dV terms) and its gradient."""
    z = np.asarray(z, float)
    r, s = z[:4], z[4:]
    q = 1.0 - r.sum()
    dp = p - 1.0
    whm = WH_HM + dp * WV_HM
    dwhm = DWH_HM + dp * DWV_HM
    wcrn, dwcrn = WCRN, DWH_CRN + dp * DWV_CRN
    d = DH + dp * DV
    e = dwhm / 2.0 + (WH_HM2 + dp * WV_HM2) / 4.0
    wo = WH + dp * WV
    h = 0.0
    g = np.zeros(7)
    dq = -np.ones(4)                           # dq/dr
    for i in range(3):
        ri, si = r[i], s[i]
        # (W_hm + dW_hm (q - r_i)) r_i q
        a = whm + dwhm * (q - ri)
        h += a * ri * q
        da = dwhm * dq.copy()
        da[i] -= dwhm
        gi = da * ri * q + a * ri * dq
        gi[i] += a * q
        g[:4] += gi
        # (W_crn + dW_crn (r_i - r3)) r_i r3
        b = wcrn + dwcrn * (ri - r[3])
        h += b * ri * r[3]
        g[i] += dwcrn * ri * r[3] + b * r[3]
        g[3] += -dwcrn * ri * r[3] + b * ri
        # (D + E q - dW_crn r3/2)(r_i^2 - s_i^2)
        c = d + e * q - dwcrn * r[3] / 2.0
        u = ri * ri - si * si
        h += c * u
        g[:4] += e * dq * u
        g[3] += -dwcrn / 2.0 * u
        g[i] += c * 2 * ri
        g[4 + i] += -c * 2 * si
        # W_o s_i^2 (r_i^2 - s_i^2)
        h += wo * si * si * u
        g[i] += wo * si * si * 2 * ri
        g[4 + i] += wo * (2 * si * u - 2 * si ** 3)
    h += WH_HM_CRN * r[3] * q
    g[:4] += WH_HM_CRN * r[3] * dq
    g[3] += WH_HM_CRN * q
    for (i, j), (w, wv, wt) in PAIRS.items():
        a, b = w + dp * wv + wt, w - wt
        h += a * r[i] * r[j] / 2.0 + b * s[i] * s[j] / 2.0
        g[i] += a * r[j] / 2.0
        g[j] += a * r[i] / 2.0
        g[4 + i] += b * s[j] / 2.0
        g[4 + j] += b * s[i] / 2.0
    return h, g


def g_and_grad(z, t, p):
    h, dh = h_and_grad(z, p)
    s, ds = entropy_and_grad(z)
    return h - t * s, dh - t * ds


def _bounds(r):
    return [(-ri + 1e-15, ri - 1e-15) if ri > 1e-14 else (0.0, 0.0) for ri in r[:3]]


def order(r, t, p, tol=1e-10):
    """s minimizing G at fixed r, from MELTS's ordered first guess s = 0.9 r, then Newton-polished
    so that s (and the chemical potentials) are smooth in r."""
    r = np.asarray(r, float)
    free = [i for i in range(3) if r[i] > 1e-14]
    s = np.zeros(3)
    if not free:
        return s, g_and_grad(np.concatenate([r, s]), t, p)[0]

    def f(sf):
        s_ = np.zeros(3)
        s_[free] = sf
        g, dg = g_and_grad(np.concatenate([r, s_]), t, p)
        return g / (R * t), dg[4:][free] / (R * t)

    bnds = [(-r[i] * (1 - 1e-12), r[i] * (1 - 1e-12)) for i in free]
    res = minimize(f, 0.9 * r[free], jac=True, method="L-BFGS-B", bounds=bnds,
                   options=dict(ftol=1e-15, gtol=1e-12, maxiter=500))
    sf = res.x
    idx = [4 + i for i in free]
    a = _A[:, idx]
    for _ in range(30):                        # Newton polish
        s[free] = sf
        z = np.concatenate([r, s])
        g = g_and_grad(z, t, p)[1][idx]
        if np.max(np.abs(g)) < tol:
            break
        hh = np.empty((len(idx), len(idx)))
        for j, k in enumerate(idx):
            zp = z.copy()
            hstep = 1e-6 * max(r[k - 4], 1e-6)
            zp[k] += hstep
            hh[:, j] = (h_and_grad(zp, p)[1][idx] - h_and_grad(z, p)[1][idx]) / hstep
        v = site_fractions(z)
        pos = v > 1e-300
        w = np.where(pos, _W / np.where(pos, v, 1.0), 0.0)
        hh = hh + R * t * (a.T * w) @ a
        try:
            step = -np.linalg.solve(hh, g)
        except np.linalg.LinAlgError:
            break
        lam = 1.0
        for k, i in enumerate(free):           # stay inside |s_i| < r_i
            new = sf[k] + step[k]
            lim = r[i] * (1 - 1e-12)
            if abs(new) > lim:
                lam = min(lam, 0.5 * (lim - abs(sf[k])) / max(abs(step[k]), 1e-300))
        sf = sf + lam * step
    s[free] = sf
    return s, g_and_grad(np.concatenate([r, s]), t, p)[0]


def _pure_order_g(t, p):
    """G of pure ilmenite-type ordering (MELTS's IL_G, GK_G, PY_G): Newton from s = 0.98 on [0, 1);
    the disordered s = 0 if it does not converge (above the first-order transition)."""
    d, wo = DH + (p - 1) * DV, WH + (p - 1) * WV

    def g(s):
        conf = -R * ((1 + s) * math.log(1 + s) + (1 - s) * math.log(1 - s) - 2 * math.log(2)) if s < 1 else 0.0
        return d * (1 - s * s) + wo * s * s * (1 - s * s) - t * conf

    s, old = 0.98, 2.0
    for it in range(1000):
        if abs(s - old) <= 10 * 2.2e-16:
            break
        old = s
        dg = R * t * (math.log(1 + s) - math.log(1 - s)) - 2 * d * s + wo * (2 * s - 4 * s ** 3)
        d2 = R * t * (1 / (1 + s) + 1 / (1 - s)) - 2 * d + wo * (2 - 12 * s * s)
        s = min(max(s - dg / d2, 0.0), 1 - 10 * 2.2e-16)
    else:
        s = 0.0
    return g(s)


@lru_cache(maxsize=512)
def end_member_g(t: float, p: float) -> np.ndarray:
    """Reference G of the five MELTS end members (geikielite, hematite, ilmenite, pyrophanite, corundum)."""
    il = _pure_order_g(t, p)                   # the three titanates share their ordering parameters
    sro = -t * C_SRO * 2 * R * math.log(2.0)
    return np.array([il, sro, il, il, sro])


def r_from_x(x):
    gk, hm, il, cn = x
    return np.array([il, gk, 0.0, cn])


def x_from_r(r):
    return np.array([r[1], 1 - r[0] - r[1] - r[3], r[0], r[3]])


_XR = np.array([[0, 1, 0], [-1, -1, -1], [1, 0, 0], [0, 0, 1]], float)   # x = _XR @ (r0, r1, r3) + e_hm
_RX = np.array([r_from_x(e) for e in np.eye(4)])
_IDX = [MELTS_ORDER.index(e) for e in END_MEMBERS]


class Rhombohedral:
    name = "rhombohedral"
    end_members = END_MEMBERS
    n_params = 4                             # log N, r0 (ilm), r1 (gei), r3 (crn)

    def g0(self, t, p) -> np.ndarray:
        return np.array([MINERALS[e].gibbs(t, p) for e in END_MEMBERS])

    def ends(self, t, p):
        return end_member_g(t, p)[_IDX]

    def g_mix(self, x, t, p) -> float:
        _, g = order(r_from_x(x), t, p)
        return g - float(np.asarray(x) @ self.ends(t, p))

    def gibbs(self, x, t, p) -> float:
        x = np.asarray(x, float)
        return float(x @ self.g0(t, p)) + self.g_mix(x, t, p)

    def mu(self, x, t, p) -> np.ndarray:
        r = r_from_x(x)
        s, _ = order(r, t, p)
        g, dg = g_and_grad(np.concatenate([r, s]), t, p)
        return self.g0(t, p) - self.ends(t, p) + g + (_RX - r) @ dg[:4]

    @staticmethod
    def feasible_x(x) -> bool:
        return bool(np.min(x) > -1e-12)

    def n_from_params(self, prm) -> np.ndarray:
        r = np.array([prm[1], prm[2], 0.0, prm[3]])
        return math.exp(prm[0]) * x_from_r(r)

    def params_from(self, x, n_total) -> np.ndarray:
        r = r_from_x(x)
        return np.array([math.log(n_total), r[0], r[1], r[3]])

    def affinity(self, mu_liq, t, p):
        """(A, x): the largest affinity over compositions; A > 0 means the liquid is supersaturated."""
        from .phases import mu_in_liquid
        mut = np.array([mu_in_liquid(e, mu_liq) for e in END_MEMBERS])
        dd = self.g0(t, p) - self.ends(t, p) - mut

        def f(y):                            # y = (r0, r1, r3), s at its minimum
            r = np.array([y[0], y[1], 0.0, y[2]])
            s, g = order(r, t, p)
            dg = g_and_grad(np.concatenate([r, s]), t, p)[1]
            val = float(x_from_r(r) @ dd) + g
            grad = _XR.T @ dd + dg[[0, 1, 3]]
            return val / (R * t), grad / (R * t)

        cons = [{"type": "ineq", "fun": lambda y: np.array([1 - y.sum() - 1e-10]), "jac": lambda y: -np.ones((1, 3))}]
        best = None
        for y0 in ([0.85, 0.08, 0.01], [0.5, 0.05, 0.01], [0.15, 0.02, 0.01]):
            res = minimize(f, y0, jac=True, method="SLSQP", constraints=cons,
                           bounds=[(1e-10, 1.0)] * 3, options=dict(maxiter=300, ftol=1e-12))
            if np.all(np.isfinite(res.x)) and (best is None or res.fun < best.fun):
                best = res
        if best is None:
            return -np.inf, x_from_r([0.85, 0.08, 0.0, 0.01])
        r = np.array([best.x[0], best.x[1], 0.0, best.x[2]])
        return -best.fun * R * t, x_from_r(r)

    @staticmethod
    def formula(x) -> str:
        """Ilmenite, geikielite, hematite and corundum mol%."""
        gk, hm, il, cn = x
        return f"Ilm{100 * il:.0f} Gk{100 * gk:.0f} Hem{100 * hm:.0f} Crn{100 * cn:.0f}"


RHOMBOHEDRAL = Rhombohedral()
