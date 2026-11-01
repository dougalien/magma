"""Spinel: the MELTS model of Sack and Ghiorso (1991, with the 1991 erratum), written for Magma.

End members (MELTS order; standard states in data/melts_minerals.tsv from SpinelBerman.m):
    chromite FeCr2O4, hercynite FeAl2O4, magnetite Fe3O4, spinel MgAl2O4, ulvospinel Fe2TiO4.

Composition: r0 = X_spinel, r1 = X_chromite, r2 = X_ulvospinel, r3 = X_magnetite (hercynite is
the rest).  Ordering: s0 (Mg between tetrahedral and octahedral sites), s1 (Al), s2 (Fe3+).
Site occupancies, per formula unit with one tetrahedral and two octahedral sites:

    tet: Mg (r0+s0)/2, Fe2+ r2 - r0/2 - s0/2 + s1 + r1 + s2, Al 1 - r1 - r2 - r3 - s1, Fe3+ r3 - s2
    oct: Mg (r0-s0)/4, Fe2+ (2 - r0 + s0 - 2 s1 - 2 r1 - 2 s2)/4, Al (1 - r1 - r2 - r3 + s1)/2,
         Fe3+ (r3+s2)/2, Cr r1, Ti r2/2

G = quadratic Taylor expansion (data/melts_spinel_taylor.tsv) - T S_config
    + r2 r3 (WV1 r3 + WV2 r2)(P - 1), with s minimizing G.  The end members' reference states are
the same expression at each pure composition with its own ordering.  Our own implementation of
the published model's equations (Douglas's approval, 2026-11-01).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .minerals import MINERALS
from .pyroxene import Taylor

R = 8.3143
END_MEMBERS = ("chromite", "hercynite", "magnetite", "spinel", "ulvospinel")
WV1, WV2 = -0.1250, 0.1018           # J/bar, the x_usp x_mt^2 and x_usp^2 x_mt volume terms

TAY = Taylor("melts_spinel_taylor.tsv", names=("X2", "X3", "X4", "X5", "S1", "S2", "S3"))

#             r0     r1    r2    r3    s0     s1    s2       b
_TET = np.array([
    [.5,    0,    0,    0,    .5,    0,    0,     0],     # Mg
    [-.5,   1,    1,    0,    -.5,   1,    1,     0],     # Fe2+
    [0,     -1,   -1,   -1,   0,     -1,   0,     1],     # Al
    [0,     0,    0,    1,    0,     0,    -1,    0],     # Fe3+
])
_OCT = np.array([
    [.25,   0,    0,    0,    -.25,  0,    0,     0],     # Mg
    [-.25,  -.5,  0,    0,    .25,   -.5,  -.5,   .5],    # Fe2+
    [0,     -.5,  -.5,  -.5,  0,     .5,   0,     .5],    # Al
    [0,     0,    0,    .5,   0,     0,    .5,    0],     # Fe3+
    [0,     1,    0,    0,    0,     0,    0,     0],     # Cr
    [0,     0,    .5,   0,    0,     0,    0,     0],     # Ti
])
_A = np.vstack([_TET[:, :7], _OCT[:, :7]])
_B = np.concatenate([_TET[:, 7], _OCT[:, 7]])
_W = np.concatenate([np.ones(4), 2 * np.ones(6)])          # one tetrahedral, two octahedral sites


def site_fractions(z) -> np.ndarray:
    return _A @ np.asarray(z, float) + _B


def entropy_and_grad(z):
    v = site_fractions(z)
    pos = v > 1e-300
    lv = np.where(pos, np.log(np.where(pos, v, 1.0)), 0.0)
    return -R * float(_W @ (v * lv)), -R * ((_W * np.where(pos, lv + 1.0, 0.0)) @ _A)


def g_and_grad(z, t, p):
    z = np.asarray(z, float)
    s, ds = entropy_and_grad(z)
    r2, r3 = z[2], z[3]
    v = r2 * r3 * (WV1 * r3 + WV2 * r2) * (p - 1.0)
    dv = np.zeros(7)
    dv[2] = (r3 * (WV1 * r3 + WV2 * r2) + r2 * r3 * WV2) * (p - 1.0)
    dv[3] = (r2 * (WV1 * r3 + WV2 * r2) + r2 * r3 * WV1) * (p - 1.0)
    return TAY.poly(z, t, p) - t * s + v, TAY.grad(z, t, p) - t * ds + dv


def x_from_r(r) -> np.ndarray:
    r0, r1, r2, r3 = r
    return np.array([r1, 1 - r0 - r1 - r2 - r3, r3, r0, r2])


def r_from_x(x) -> np.ndarray:
    cr, hc, mt, sp, usp = x
    return np.array([sp, cr, usp, mt])


_XR = np.array([[0, 1, 0, 0], [-1, -1, -1, -1], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]], float)
_RX = np.array([r_from_x(e) for e in np.eye(5)])           # r of each pure end member


def _random_s(r):
    """The ordering of a random distribution over the available sites (MELTS's first guess)."""
    r0, r1, r2, r3 = r
    tot = dict(al=2 * (1 - r1 - r2 - r3), fe3=2 * r3, mg=r0, fe2=1 - r0 + r2)
    ratio = 2.0 - 2 * r1 - r2
    f = ratio / (1 + ratio)
    oct_ = {k: v * f for k, v in tot.items()}
    tet = {k: tot[k] - oct_[k] for k in tot}
    oct_ = {k: v / 2 for k, v in oct_.items()}
    return np.array([tet["mg"] - 2 * oct_["mg"], oct_["al"] - tet["al"] / 2, oct_["fe3"] - tet["fe3"] / 2])


_CONS_A = _A[:, 4:]


def order(r, t, p):
    """(s, G) with s minimizing G at fixed r, keeping every site occupancy positive.  Site rows that
    no ordering can change (Cr, Ti) are left out, and without Fe3+ its ordering s2 stays at zero."""
    r = np.asarray(r, float)
    base = _A[:, :4] @ r + _B
    free = [0, 1, 2] if r[3] > 1e-12 else [0, 1]
    rows = np.any(_CONS_A[:, free] != 0, axis=1)
    ca = _CONS_A[np.ix_(rows, free)]

    def full(sf):
        s = np.zeros(3)
        s[free] = sf
        return s

    def f(sf):
        g, dg = g_and_grad(np.concatenate([r, full(sf)]), t, p)
        return g, dg[4:][free]

    key = (t, p, len(free))
    for r_old, sf_old in _CACHE.get(key, ()):   # a nearby composition: Newton from its ordering
        if np.max(np.abs(r_old - r)) < 0.02:
            sf = _polish(r, sf_old.copy(), free, t, p)
            z = np.concatenate([r, full(sf)])
            if (site_fractions(z)[rows].min() > 0
                    and np.max(np.abs(g_and_grad(z, t, p)[1][4:][free])) < 1e-8):
                return full(sf), g_and_grad(z, t, p)[0]
            break
    cons = [{"type": "ineq", "fun": lambda sf: base[rows] + ca @ sf - 1e-12, "jac": lambda sf: ca}]
    s0 = _random_s(r)[free]
    res = minimize(f, s0, jac=True, method="SLSQP", constraints=cons, options=dict(ftol=1e-14, maxiter=200))
    sf = res.x if np.all(base[rows] + ca @ res.x > -1e-10) else s0
    sf = _polish(r, sf, free, t, p)
    entries = _CACHE.setdefault(key, [])
    entries.insert(0, (r.copy(), sf.copy()))
    del entries[8:]
    if len(_CACHE) > 64:
        _CACHE.clear()
    s = full(sf)
    return s, g_and_grad(np.concatenate([r, s]), t, p)[0]


_CACHE: dict = {}       # (T, P, number of free s) -> recent (r, s): orderings to start Newton from


def _polish(r, sf, free, t, p, tol=1e-10):
    """Newton steps on dG/ds = 0 after SLSQP, so that s, and the chemical potentials that depend on
    it, are smooth functions of r (the solver differentiates them numerically)."""
    idx = [4 + k for k in free]
    a = _A[:, idx]
    z = np.concatenate([r, np.zeros(3)])
    # the Taylor expansion is quadratic, so its Hessian is exact by differences of its gradient
    z[idx] = sf
    g_tay = TAY.grad(z, t, p)[idx]
    h_tay = np.empty((len(idx), len(idx)))
    for j, k in enumerate(idx):
        zp = z.copy()
        zp[k] += 1e-4
        h_tay[:, j] = (TAY.grad(zp, t, p)[idx] - g_tay) / 1e-4
    for _ in range(30):
        z[idx] = sf
        g = g_and_grad(z, t, p)[1][idx]
        if np.max(np.abs(g)) < tol:
            break
        v = site_fractions(z)
        pos = v > 1e-300
        w = np.where(pos, _W / np.where(pos, v, 1.0), 0.0)
        h = h_tay + R * t * (a.T * w) @ a
        try:
            step = -np.linalg.solve(h, g)
        except np.linalg.LinAlgError:
            break
        dv = a @ step
        shrink = [-0.5 * vi / di for vi, di in zip(v, dv) if di < 0 and vi > 0]
        lam = min([1.0] + shrink)
        sf = sf + lam * step
    return sf


def end_member_g(t, p) -> np.ndarray:
    def g(r, s):
        return g_and_grad(np.concatenate([r, s]), t, p)[0]

    eps = 1e-12
    cr = g([0, 1, 0, 0], [0, 0, 0])
    usp = g([0, 0, 1, 0], [0, 0, 0])
    hc = minimize_scalar(lambda v: g([0, 0, 0, 0], [0, v, 0]), bounds=(eps, 1 - eps), method="bounded",
                         options=dict(xatol=1e-12)).fun
    sp = minimize_scalar(lambda v: g([1, 0, 0, 0], [v, (1 + v) / 2, 0]), bounds=(-1 + eps, 1 - eps),
                         method="bounded", options=dict(xatol=1e-12)).fun
    mt = minimize_scalar(lambda v: g([0, 0, 0, 1], [0, 0, v]), bounds=(eps, 1 - eps), method="bounded",
                         options=dict(xatol=1e-12)).fun
    return np.array([cr, hc, mt, sp, usp])


class Spinel:
    name = "spinel"
    end_members = END_MEMBERS

    def __init__(self, inactive=()):
        """inactive: end members that cannot form (an oxide missing from the system); their
        fraction is held at zero and their equilibrium condition dropped."""
        self.inactive = tuple(inactive)
        self.active = np.array([e not in self.inactive for e in END_MEMBERS])
        self._free = [k for k in range(4) if self.active[END_MEMBERS.index(("spinel", "chromite",
                                                                            "ulvospinel", "magnetite")[k])]]
        self.n_params = 1 + len(self._free)

    def restricted(self, inactive):
        inactive = [e for e in inactive if e in END_MEMBERS and e != "hercynite"]
        return Spinel(inactive) if inactive else self

    def g0(self, t, p) -> np.ndarray:
        return np.array([MINERALS[e].gibbs(t, p) for e in END_MEMBERS])

    def g_mix(self, x, t, p) -> float:
        _, g = order(r_from_x(x), t, p)
        return g - float(np.asarray(x) @ end_member_g(t, p))

    def gibbs(self, x, t, p) -> float:
        x = np.asarray(x, float)
        return float(x @ self.g0(t, p)) + self.g_mix(x, t, p)

    def mu(self, x, t, p) -> np.ndarray:
        r = r_from_x(x)
        s, _ = order(r, t, p)
        g, dg = g_and_grad(np.concatenate([r, s]), t, p)
        return self.g0(t, p) - end_member_g(t, p) + g + (_RX - r) @ dg[:4]

    @staticmethod
    def feasible_x(x) -> bool:
        r = r_from_x(x)
        base = _A[:, :4] @ r + _B
        s, _ = order(r, 1400.0, 1.0)
        return bool(np.min(base + _CONS_A @ s) > -1e-9)

    # ---- amounts as solver parameters ------------------------------------------------------
    def n_from_params(self, prm) -> np.ndarray:
        r = np.zeros(4)
        r[self._free] = prm[1:]
        return math.exp(prm[0]) * x_from_r(r)

    def params_from(self, x, n_total) -> np.ndarray:
        return np.concatenate([[math.log(n_total)], r_from_x(x)[self._free]])

    def affinity(self, mu_liq, t, p):
        from .phases import mu_in_liquid
        mut = np.array([mu_in_liquid(e, mu_liq) if a else 0.0 for e, a in zip(END_MEMBERS, self.active)])
        d = self.g0(t, p) - end_member_g(t, p) - mut
        d[~self.active] = 0.0
        dz = np.concatenate([d @ _XR, np.zeros(3)])
        fixed = [k for k in range(4) if k not in self._free]

        def f(z):
            g, dg = g_and_grad(z, t, p)
            return float(x_from_r(z[:4]) @ d) + g, dg + dz

        def fs(z):                               # in units of RT, which SLSQP handles far better
            g, dg = f(z)
            return g / (R * t), dg / (R * t)

        if not self.active[END_MEMBERS.index("magnetite")]:
            fixed = fixed + [6]                  # no Fe3+, so no Fe3+ ordering
        # site rows that stay at zero once the fixed variables are zero cannot be held positive
        rows = np.any(np.delete(_A, fixed, axis=1) != 0, axis=1)
        cons = [{"type": "ineq", "fun": lambda z: site_fractions(z)[rows] - 1e-10, "jac": lambda z: _A[rows]}]
        if fixed:
            cons.append({"type": "eq", "fun": lambda z: z[fixed], "jac": lambda z: np.eye(7)[fixed]})
        best = None
        for x0 in _STARTS:
            r0 = r_from_x(x0)
            r0[[k for k in fixed if k < 4]] = 0.0
            z0 = np.concatenate([r0, _random_s(r0)])
            res = minimize(fs, z0, jac=True, method="SLSQP", constraints=cons, options=dict(maxiter=500, ftol=1e-12))
            if (np.all(np.isfinite(res.x)) and site_fractions(res.x).min() > -1e-8
                    and (best is None or res.fun < best.fun)):
                best = res
        if best is None:
            return -np.inf, _STARTS[0]
        r = best.x[:4]
        return -float(f(np.concatenate([r, order(r, t, p)[0]]))[0]), x_from_r(r)

    @staticmethod
    def formula(x) -> str:
        """Cations per four oxygens."""
        cr, hc, mt, sp, usp = x
        cat = [("Fe2+", hc + cr + mt + 2 * usp), ("Mg", sp), ("Fe3+", 2 * mt), ("Al", 2 * (hc + sp)),
               ("Cr", 2 * cr), ("Ti", usp)]
        return " ".join(f"{k}{v:.2f}" for k, v in cat if v > 0.005) + " O4"


_STARTS = [np.array(v, float) for v in (
    [0.20, 0.25, 0.10, 0.40, 0.05],      # Cr-Al spinel of a basalt
    [0.02, 0.15, 0.45, 0.08, 0.30],      # titanomagnetite
    [0.05, 0.45, 0.05, 0.40, 0.05],      # aluminous spinel
)]

SPINEL = Spinel()
