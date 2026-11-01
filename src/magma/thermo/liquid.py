"""The silicate liquid of Ghiorso and Sack (1995): components and the Gibbs energy of mixing.

    G_mix = RT sum X_i ln X_i + 1/2 sum_i sum_j W_ij X_i X_j
            + RT [X_w ln X_w + (1 - X_w) ln(1 - X_w)]                      (eq. 1)

    mu_k - mu_k0 = RT ln X_k + sum_i X_i W_ik - 1/2 sum sum W_ij X_i X_j + RT ln(1 - X_w)   (2a)
    mu_w - mu_w0 = 2 RT ln X_w + sum_i X_i W_iw - 1/2 sum sum W_ij X_i X_j                   (2b)

The W_ij are read from data/gs1995_liquid_w.tsv (Table 4).  The standard-state chemical
potentials mu_k0 (Table A1, anchored to Berman 1988 solids) come in a later module.

MnO, NiO and CoO are not components of the 1995 model (rhyolite-MELTS adds them later); they
are dropped with a note when a composition is converted.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from importlib import resources

import numpy as np

from ..composition import OXIDES, Composition

R = 8.3143  # J/(mol K), as in MELTS

COMPONENTS = ("SiO2", "TiO2", "Al2O3", "Fe2O3", "MgCr2O4", "Fe2SiO4", "Mg2SiO4", "CaSiO3",
              "Na2SiO3", "KAlSiO4", "Ca3(PO4)2", "H2O")
WATER = COMPONENTS.index("H2O")

# oxide moles in one mole of each component
_STOICH = {
    "SiO2": {"SiO2": 1}, "TiO2": {"TiO2": 1}, "Al2O3": {"Al2O3": 1}, "Fe2O3": {"Fe2O3": 1},
    "MgCr2O4": {"MgO": 1, "Cr2O3": 1}, "Fe2SiO4": {"FeO": 2, "SiO2": 1},
    "Mg2SiO4": {"MgO": 2, "SiO2": 1}, "CaSiO3": {"CaO": 1, "SiO2": 1},
    "Na2SiO3": {"Na2O": 1, "SiO2": 1}, "KAlSiO4": {"K2O": 0.5, "Al2O3": 0.5, "SiO2": 1},
    "Ca3(PO4)2": {"CaO": 3, "P2O5": 1}, "H2O": {"H2O": 1},
}
MODEL_OXIDES = ("SiO2", "TiO2", "Al2O3", "Fe2O3", "Cr2O3", "FeO", "MgO", "CaO", "Na2O", "K2O",
                "P2O5", "H2O")
IGNORED_OXIDES = tuple(ox for ox in OXIDES if ox not in MODEL_OXIDES)

# rows: oxides, columns: components
STOICH = np.array([[_STOICH[c].get(ox, 0.0) for c in COMPONENTS] for ox in MODEL_OXIDES])
_TO_COMPONENTS = np.linalg.inv(STOICH)          # square 12 x 12: the conversion is exact


def _load_w() -> np.ndarray:
    w = np.zeros((len(COMPONENTS), len(COMPONENTS)))
    idx = {c: k for k, c in enumerate(COMPONENTS)}
    text = resources.files("magma.data").joinpath("gs1995_liquid_w.tsv").read_text("utf-8")
    rows = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    for row in csv.DictReader(rows, delimiter="\t"):
        i, j = idx[row["i"]], idx[row["j"]]
        w[i, j] = w[j, i] = float(row["W_J"])
    return w


W = _load_w()


@dataclass(frozen=True)
class LiquidComposition:
    """Moles of the liquid components, and the oxides left out of the model."""
    moles: np.ndarray
    ignored: dict[str, float]

    @property
    def total(self) -> float:
        return float(self.moles.sum())

    @property
    def x(self) -> np.ndarray:
        return self.moles / self.moles.sum()

    def as_dict(self) -> dict[str, float]:
        return dict(zip(COMPONENTS, self.moles.tolist()))


def to_components(comp: Composition) -> LiquidComposition:
    """Oxide wt% (per 100 g as given) to moles of liquid components."""
    n = comp.moles()
    oxide_moles = np.array([n.get(ox, 0.0) for ox in MODEL_OXIDES])
    moles = _TO_COMPONENTS @ oxide_moles
    moles[np.abs(moles) < 1e-14] = 0.0
    ignored = {ox: comp.get(ox) for ox in IGNORED_OXIDES if comp.get(ox) > 0}
    return LiquidComposition(moles, ignored)


def from_components(moles: np.ndarray) -> dict[str, float]:
    """Moles of liquid components back to moles of oxides."""
    return dict(zip(MODEL_OXIDES, (STOICH @ np.asarray(moles)).tolist()))


def gibbs_mixing(x: np.ndarray, t_k: float) -> float:
    """Molar Gibbs energy of mixing (ideal + regular solution), J/mol; eq. 1."""
    x = np.asarray(x, dtype=float)
    pos = x > 0
    ideal = R * t_k * float(np.sum(x[pos] * np.log(x[pos])))
    xw = x[WATER]
    if 0 < xw < 1:
        ideal += R * t_k * (xw * np.log(xw) + (1 - xw) * np.log(1 - xw))
    return ideal + 0.5 * float(x @ W @ x)


def mu_mixing(x: np.ndarray, t_k: float) -> np.ndarray:
    """mu_k - mu_k0 for every component, J/mol; eqs. 2a and 2b.  -inf where X_k = 0."""
    x = np.asarray(x, dtype=float)
    excess = W @ x - 0.5 * float(x @ W @ x)
    with np.errstate(divide="ignore"):
        lnx = np.log(x)
    xw = x[WATER]
    mu = R * t_k * lnx + excess + R * t_k * np.log1p(-xw)
    mu[WATER] = 2 * R * t_k * lnx[WATER] + excess[WATER]
    return mu


def activities(x: np.ndarray, t_k: float) -> np.ndarray:
    """a_k = exp((mu_k - mu_k0) / RT)."""
    return np.exp(mu_mixing(x, t_k) / (R * t_k))
