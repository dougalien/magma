"""Holland & Powell solid solutions (the igneous set of Holland, Green and Powell 2018 and its sources).

Each model is a list of end members, each made from dataset end members plus a DQF term
(g = sum c_k G_k + a + b T + c P), a set of site fractions that are linear in the end-member
proportions p (sf = A p, the end members' site occupancies), the ideal activity of each end member as
a product of site fractions, and a non-ideal term:

    mu_i = g_i + RT ln a_i(sf) + mu_ex_i(p),    G = sum p_i mu_i

mu_ex is the symmetric formalism (Powell and Holland 1993) or, when volumes v_i are given, the
asymmetric van Laar form (Holland and Powell 2003), with pair parameters W_ij = a + b T + c P (kJ).
Proportions may be negative for order-disorder end members, as long as every site fraction is >= 0.

Parameter values are those of the igneous set as distributed in MAGEMin (numbers only, Douglas's
approval 2026-11-02).  Units: kJ, K, kbar.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .endmembers import DATA, OXIDES, R, gibbs

# a recipe term: (coefficient, dataset end member, order state)
Term = tuple


@dataclass(frozen=True)
class EndMemberDef:
    name: str
    recipe: tuple                 # ((coef, dataset name, state), ...)
    dqf: tuple = (0.0, 0.0, 0.0)  # a + b T + c P, kJ

    def g(self, t: float, p: float) -> float:
        a, b, c = self.dqf
        return sum(k * gibbs(n, t, p, s) for k, n, s in self.recipe) + a + b * t + c * p

    def composition(self) -> np.ndarray:
        """Moles of each oxide in OXIDES."""
        out = np.zeros(len(OXIDES))
        for k, n, _ in self.recipe:
            out += k * np.array([DATA[n].comp.get(o, 0.0) for o in OXIDES])
        return out


def em(name, *recipe, dqf=(0.0, 0.0, 0.0)):
    """em("obi", (1/3, "ann"), (2/3, "phl"), dqf=(-6, 0, 0)); a term may carry a state as third item."""
    terms = tuple((t[0], t[1], t[2] if len(t) > 2 else "equilibrium") for t in recipe)
    return EndMemberDef(name, terms, dqf)


@dataclass(frozen=True)
class Activity:
    """a = k * prod sf[j]^e_j"""
    k: float
    powers: tuple                 # ((site-fraction name, exponent), ...)


@dataclass
class HPSolution:
    name: str
    end_members: tuple            # EndMemberDef
    site_names: tuple
    occupancy: dict               # end-member name -> {site name: fraction}
    activities: dict              # end-member name -> Activity
    w: dict                       # (i, j) -> (a, b, c) for W = a + b T + c P
    v: tuple | None = None        # van Laar volumes (asymmetric), or None (symmetric)
    _a: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        n, m = len(self.end_members), len(self.site_names)
        idx = {s: j for j, s in enumerate(self.site_names)}
        self._a = np.zeros((m, n))
        for i, e in enumerate(self.end_members):
            for s, f in self.occupancy[e.name].items():
                self._a[idx[s], i] = f
        self._act = []
        for e in self.end_members:
            act = self.activities[e.name]
            self._act.append((act.k, [(idx[s], x) for s, x in act.powers]))
        self._w = np.zeros((n, n, 3))
        for (i, j), abc in self.w.items():
            self._w[i, j] = self._w[j, i] = abc

    @property
    def names(self) -> tuple:
        return tuple(e.name for e in self.end_members)

    @property
    def n(self) -> int:
        return len(self.end_members)

    def g0(self, t: float, p: float) -> np.ndarray:
        return np.array([e.g(t, p) for e in self.end_members])

    def site_fractions(self, prop) -> np.ndarray:
        return self._a @ np.asarray(prop, float)

    def ln_activity_ideal(self, prop) -> np.ndarray:
        sf = self.site_fractions(prop)
        out = np.empty(self.n)
        for i, (k, pw) in enumerate(self._act):
            v = math.log(k)
            for j, x in pw:
                v += x * math.log(sf[j]) if sf[j] > 0 else -math.inf
            out[i] = v
        return out

    def w_matrix(self, t: float, p: float) -> np.ndarray:
        return self._w[:, :, 0] + self._w[:, :, 1] * t + self._w[:, :, 2] * p

    def mu_excess(self, prop, t: float, p: float) -> np.ndarray:
        prop = np.asarray(prop, float)
        w = self.w_matrix(t, p)
        n = self.n
        if self.v is None:
            phi, scale = prop, np.ones(n)
            wstar = w
        else:
            v = np.asarray(self.v, float)
            phi = prop * v / (prop @ v)
            wstar = w * 2.0 / (v[:, None] + v[None, :])
            scale = v
        out = np.empty(n)
        for k in range(n):
            q = -phi.copy()
            q[k] += 1.0
            out[k] = -scale[k] * sum(q[i] * q[j] * wstar[i, j] for i in range(n) for j in range(i + 1, n))
        return out

    def g_excess(self, prop, t: float, p: float) -> float:
        prop = np.asarray(prop, float)
        w = self.w_matrix(t, p)
        if self.v is None:
            return float(sum(prop[i] * prop[j] * w[i, j] for i in range(self.n) for j in range(i + 1, self.n)))
        v = np.asarray(self.v, float)
        vt = prop @ v
        phi = prop * v / vt
        return float(vt * sum(phi[i] * phi[j] * w[i, j] * 2.0 / (v[i] + v[j])
                              for i in range(self.n) for j in range(i + 1, self.n)))

    def mu(self, prop, t: float, p: float) -> np.ndarray:
        """Chemical potentials of the end members, kJ/mol, at proportions prop (sum 1)."""
        return self.g0(t, p) + R * t * self.ln_activity_ideal(prop) + self.mu_excess(prop, t, p)

    def gibbs(self, prop, t: float, p: float) -> float:
        prop = np.asarray(prop, float)
        return float(prop @ self.mu(prop, t, p))

    def compositions(self) -> np.ndarray:
        """(n_em, n_oxides) moles of oxides per end member."""
        return np.array([e.composition() for e in self.end_members])


def A(k=1.0, **powers):
    return Activity(k, tuple(powers.items()))


def _pairs(n, values):
    """Map the flat list W[0..] (pairs (0,1), (0,2), ..., (n-2, n-1)) to {(i, j): (a, b, c)}."""
    keys = [(i, j) for i in range(n) for j in range(i + 1, n)]
    assert len(keys) == len(values)
    return {k: (v if isinstance(v, tuple) else (v, 0.0, 0.0)) for k, v in zip(keys, values)}


# ---- olivine (Holland, Green and Powell 2018) ------------------------------------------------------------
OLIVINE = HPSolution(
    "ol",
    (em("mont", (1, "mont")), em("fa", (1, "fa")), em("fo", (1, "fo")), em("cfm", (0.5, "fa"), (0.5, "fo"))),
    ("xMgM1", "xFeM1", "xMgM2", "xFeM2", "xCaM2"),
    {"mont": {"xMgM1": 1, "xCaM2": 1}, "fa": {"xFeM1": 1, "xFeM2": 1}, "fo": {"xMgM1": 1, "xMgM2": 1},
     "cfm": {"xMgM1": 1, "xFeM2": 1}},
    {"mont": A(xMgM1=1, xCaM2=1), "fa": A(xFeM1=1, xFeM2=1), "fo": A(xMgM1=1, xMgM2=1), "cfm": A(xMgM1=1, xFeM2=1)},
    _pairs(4, [24.0, 38.0, 24.0, 9.0, 4.5, 4.5]),
)

# ---- feldspar (Holland, Green and Powell 2018, after Holland and Powell 2003) -------------------------------
FELDSPAR = HPSolution(
    "fsp",
    (em("ab", (1, "ab")), em("an", (1, "an")), em("san", (1, "san"))),
    ("xNaA", "xCaA", "xKA", "xAlTB", "xSiTB"),
    {"ab": {"xNaA": 1, "xAlTB": 0.25, "xSiTB": 0.75}, "an": {"xCaA": 1, "xAlTB": 0.5, "xSiTB": 0.5},
     "san": {"xKA": 1, "xAlTB": 0.25, "xSiTB": 0.75}},
    {"ab": A(1.7548, xNaA=1, xAlTB=0.25, xSiTB=0.75), "an": A(2.0, xCaA=1, xAlTB=0.5, xSiTB=0.5),
     "san": A(1.7548, xKA=1, xAlTB=0.25, xSiTB=0.75)},
    _pairs(3, [(14.6, -0.00935, -0.04), (24.1, -0.00957, 0.338), (48.5, 0.0, -0.13)]),
    v=(0.674, 0.55, 1.0),
)

# ---- spinel (Holland, Green and Powell 2018) ---------------------------------------------------------------
_SPL_T = "xMgT xFeT xAlT xFe3T".split()
_SPL_M = "xMgM xFeM xAlM xFe3M xCrM xTiM".split()
SPINEL = HPSolution(
    "spl",
    (em("nsp", (1, "sp", "ordered")),
     em("isp", (1, "sp", "ordered"), dqf=(23.5, -0.005763, 0.0)),
     em("nhc", (1, "herc", "ordered")),
     em("ihc", (1, "herc", "ordered"), dqf=(23.6, -0.005763, 0.0)),
     em("nmt", (1, "mt"), dqf=(0.0, 0.005763, 0.0)),
     em("imt", (1, "mt"), dqf=(0.3, 0.0, 0.0)),
     em("pcr", (1, "picr")),
     em("usp", (1, "usp"))),          # MAGEMin lists this end member as "qndm"; it is made from usp
    tuple(_SPL_T + _SPL_M),
    {"nsp": {"xMgT": 1, "xAlM": 1}, "isp": {"xAlT": 1, "xMgM": 0.5, "xAlM": 0.5},
     "nhc": {"xFeT": 1, "xAlM": 1}, "ihc": {"xAlT": 1, "xFeM": 0.5, "xAlM": 0.5},
     "nmt": {"xFeT": 1, "xFe3M": 1}, "imt": {"xFe3T": 1, "xFeM": 0.5, "xFe3M": 0.5},
     "pcr": {"xMgT": 1, "xCrM": 1}, "usp": {"xFeT": 1, "xFeM": 0.5, "xTiM": 0.5}},
    {"nsp": A(xMgT=1, xAlM=1), "isp": A(2.0, xAlT=1, xMgM=0.5, xAlM=0.5),
     "nhc": A(xFeT=1, xAlM=1), "ihc": A(2.0, xAlT=1, xFeM=0.5, xAlM=0.5),
     "nmt": A(xFeT=1, xFe3M=1), "imt": A(2.0, xFe3T=1, xFeM=0.5, xFe3M=0.5),
     "pcr": A(xMgT=1, xCrM=1), "usp": A(2.0, xFeT=1, xFeM=0.5, xTiM=0.5)},
    _pairs(8, [-6.7, 3.6, -9.8, 43.2, 49.5, (-38.4, 0.0, -0.08), 40.0,
               2.7, -3.5, 36.8, 20.7, (-21.6, 0.0, -0.08), 38.2,
               -6.0, 17.5, 51.6, -53.8, 25.7,
               -4.1, 10.0, -38.8, 21.0,
               18.1, 12.1, 5.2,
               -8.7, 21.5,
               15.0]),
    v=(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9),
)

# ---- ilmenite-hematite-geikielite (Holland, Green and Powell 2018) ----------------------------------------
ILMENITE = HPSolution(
    "ilm",
    (em("oilm", (1, "ilm", "ordered")), em("dilm", (1, "ilm", "disordered")), em("hm", (1, "hem")),
     em("ogk", (1, "geik", "ordered")), em("dgk", (1, "geik", "disordered"))),
    ("xFeA", "xTiA", "xFe3A", "xMgA", "xFeB", "xTiB", "xFe3B", "xMgB"),
    {"oilm": {"xFeA": 1, "xTiB": 1},
     "dilm": {"xFeA": 0.5, "xTiA": 0.5, "xFeB": 0.5, "xTiB": 0.5},
     "hm": {"xFe3A": 1, "xFe3B": 1},
     "ogk": {"xMgA": 1, "xTiB": 1},
     "dgk": {"xMgA": 0.5, "xTiA": 0.5, "xMgB": 0.5, "xTiB": 0.5}},
    {"oilm": A(xFeA=0.5, xTiB=0.5), "dilm": A(2.0, xFeA=0.25, xTiA=0.25, xFeB=0.25, xTiB=0.25),
     "hm": A(xFe3A=0.5, xFe3B=0.5), "ogk": A(xMgA=0.5, xTiB=0.5),
     "dgk": A(2.0, xTiA=0.25, xMgA=0.25, xTiB=0.25, xMgB=0.25)},
    _pairs(5, [(7.05, 0.0, 0.13), 14.3, -7.6, 0.6, (7.25, 0.0, -0.13), -5.5, -2.2, 12.5, 2.7, 8.3]),
)

# ---- garnet (Holland, Green and Powell 2018) ----------------------------------------------------------------
GARNET = HPSolution(
    "g",
    (em("py", (1, "py")), em("alm", (1, "alm")), em("gr", (1, "gr")), em("andr", (1, "andr")),
     em("knom", (1, "knor")),
     em("tig", (1, "py"), (-0.5, "cor"), (0.5, "per"), (0.5, "ru"), dqf=(42.3, -0.0173, 0.0))),
    ("xMgM1", "xFeM1", "xCaM1", "xAlM2", "xCrM2", "xFe3M2", "xMgM2", "xTiM2"),
    {"py": {"xMgM1": 1, "xAlM2": 1}, "alm": {"xFeM1": 1, "xAlM2": 1}, "gr": {"xCaM1": 1, "xAlM2": 1},
     "andr": {"xCaM1": 1, "xFe3M2": 1}, "knom": {"xMgM1": 1, "xCrM2": 1},
     "tig": {"xMgM1": 1, "xAlM2": 0.5, "xMgM2": 0.25, "xTiM2": 0.25}},
    {"py": A(xMgM1=3, xAlM2=2), "alm": A(xFeM1=3, xAlM2=2), "gr": A(xCaM1=3, xAlM2=2),
     "andr": A(xCaM1=3, xFe3M2=2), "knom": A(xMgM1=3, xCrM2=2), "tig": A(8.0, xMgM1=3, xAlM2=1, xMgM2=0.5, xTiM2=0.5)},
    _pairs(6, [(4.0, 0.0, 0.1), (45.4, -0.01, 0.04), (107.0, -0.01, -0.036), 2.0, 1.0,
               (17.0, -0.01, 0.1), (65.0, -0.01, 0.039), (8.2, 0.0, 0.01), 0.0,
               2.0, (5.0, -0.01, 0.18), -3.0,
               (63.0, -0.01, 0.1), -1.0,
               0.0]),
    v=(1.0, 1.0, 2.5, 2.5, 1.0, 1.0),
)

# ---- biotite (Holland, Green and Powell 2018, after White et al. 2014) -------------------------------------
BIOTITE = HPSolution(
    "bi",
    (em("phl", (1, "phl")), em("annm", (1, "ann"), dqf=(-3.0, 0.0, 0.0)),
     em("obi", (1 / 3, "ann"), (2 / 3, "phl"), dqf=(-6.0, 0.0, 0.0)),
     em("eas", (1, "east"), dqf=(2.2, 0.0, 0.0)),
     em("tbi", (1, "phl"), (-1, "br"), (1, "ru"), dqf=(55.0, 0.0, 0.0)),
     em("fbi", (1, "acm"), (1, "east"), (-1, "jd"), dqf=(-7.2, 0.0, 0.0))),
    ("xMgM3", "xFeM3", "xFe3M3", "xTiM3", "xAlM3", "xMgM12", "xFeM12", "xSiT", "xAlT", "xOHV", "xOV"),
    {"phl": {"xMgM3": 1, "xMgM12": 1, "xSiT": 0.5, "xAlT": 0.5, "xOHV": 1},
     "annm": {"xFeM3": 1, "xFeM12": 1, "xSiT": 0.5, "xAlT": 0.5, "xOHV": 1},
     "obi": {"xFeM3": 1, "xMgM12": 1, "xSiT": 0.5, "xAlT": 0.5, "xOHV": 1},
     "eas": {"xAlM3": 1, "xMgM12": 1, "xAlT": 1, "xOHV": 1},
     "tbi": {"xTiM3": 1, "xMgM12": 1, "xSiT": 0.5, "xAlT": 0.5, "xOV": 1},
     "fbi": {"xFe3M3": 1, "xMgM12": 1, "xAlT": 1, "xOHV": 1}},
    {"phl": A(4.0, xMgM3=1, xMgM12=2, xSiT=1, xAlT=1, xOHV=2),
     "annm": A(4.0, xFeM3=1, xFeM12=2, xSiT=1, xAlT=1, xOHV=2),
     "obi": A(4.0, xFeM3=1, xMgM12=2, xSiT=1, xAlT=1, xOHV=2),
     "eas": A(xAlM3=1, xMgM12=2, xAlT=2, xOHV=2),
     "tbi": A(4.0, xOV=2, xTiM3=1, xMgM12=2, xSiT=1, xAlT=1),
     "fbi": A(xFe3M3=1, xMgM12=2, xAlT=2, xOHV=2)},
    _pairs(6, [12.0, 4.0, 10.0, 30.0, 14.0, 8.0, 0.0, 32.0, 4.0, 7.0, 24.0, 7.0, 40.0, 1.0, 40.0]),
)

# ---- white mica (Holland, Green and Powell 2018, after White et al. 2014) ----------------------------------
MUSCOVITE = HPSolution(
    "mu",
    (em("mu", (1, "mu")), em("cel", (1, "cel")), em("fcel", (1, "fcel")), em("pa", (1, "pa")),
     em("mam", (1, "ma"), dqf=(5.0, 0.0, 0.0)),
     em("fmu", (0.5, "andr"), (-0.5, "gr"), (1, "mu"), dqf=(25.0, 0.0, 0.0))),
    ("xKA", "xNaA", "xCaA", "xMgM2A", "xFeM2A", "xAlM2A", "xAlM2B", "xFe3M2B", "xSiT1", "xAlT1"),
    {"mu": {"xKA": 1, "xAlM2A": 1, "xAlM2B": 1, "xSiT1": 0.5, "xAlT1": 0.5},
     "cel": {"xKA": 1, "xMgM2A": 1, "xAlM2B": 1, "xSiT1": 1},
     "fcel": {"xKA": 1, "xFeM2A": 1, "xAlM2B": 1, "xSiT1": 1},
     "pa": {"xNaA": 1, "xAlM2A": 1, "xAlM2B": 1, "xSiT1": 0.5, "xAlT1": 0.5},
     "mam": {"xCaA": 1, "xAlM2A": 1, "xAlM2B": 1, "xAlT1": 1},
     "fmu": {"xKA": 1, "xAlM2A": 1, "xFe3M2B": 1, "xSiT1": 0.5, "xAlT1": 0.5}},
    {"mu": A(4.0, xKA=1, xAlM2A=1, xAlM2B=1, xSiT1=1, xAlT1=1),
     "cel": A(xKA=1, xMgM2A=1, xAlM2B=1, xSiT1=2),
     "fcel": A(xKA=1, xFeM2A=1, xAlM2B=1, xSiT1=2),
     "pa": A(4.0, xNaA=1, xAlM2A=1, xAlM2B=1, xSiT1=1, xAlT1=1),
     "mam": A(xCaA=1, xAlM2A=1, xAlM2B=1, xAlT1=2),
     "fmu": A(4.0, xKA=1, xAlM2A=1, xFe3M2B=1, xSiT1=1, xAlT1=1)},
    _pairs(6, [(0.0, 0.0, 0.2), (0.0, 0.0, 0.2), (10.12, 0.0034, 0.353), 35.0, 0.0,
               0.0, (45.0, 0.0, 0.25), 50.0, 0.0,
               (45.0, 0.0, 0.25), 50.0, 0.0,
               15.0, 30.0,
               35.0]),
    v=(0.63, 0.63, 0.63, 0.37, 0.63, 0.63),
)

SOLUTIONS = {s.name: s for s in (OLIVINE, FELDSPAR, SPINEL, ILMENITE, GARNET, BIOTITE, MUSCOVITE)}
