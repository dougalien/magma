"""Silicate melt and rock compositions as oxides.

A ``Composition`` holds oxide weight percent and converts to moles, mole fractions and cation
fractions.  Iron can be given as FeO and Fe2O3 separately or as total iron; ``redox.py`` splits
total iron at a given oxygen fugacity.

Atomic weights: IUPAC 2007 standard values (as used by most petrology codes).
"""
from __future__ import annotations

from dataclasses import dataclass, field

ATOMIC_WEIGHT = {
    "Si": 28.0855, "Ti": 47.867, "Al": 26.9815386, "Fe": 55.845, "Cr": 51.9961,
    "Mn": 54.938045, "Mg": 24.305, "Ni": 58.6934, "Co": 58.933195, "Ca": 40.078,
    "Na": 22.98976928, "K": 39.0983, "P": 30.973762, "H": 1.00794, "C": 12.0107,
    "O": 15.9994,
}

# oxide -> (cation, cations per formula, oxygens per formula), in the usual MELTS order
OXIDES: dict[str, tuple[str, int, int]] = {
    "SiO2": ("Si", 1, 2), "TiO2": ("Ti", 1, 2), "Al2O3": ("Al", 2, 3), "Fe2O3": ("Fe", 2, 3),
    "Cr2O3": ("Cr", 2, 3), "FeO": ("Fe", 1, 1), "MnO": ("Mn", 1, 1), "MgO": ("Mg", 1, 1),
    "NiO": ("Ni", 1, 1), "CoO": ("Co", 1, 1), "CaO": ("Ca", 1, 1), "Na2O": ("Na", 2, 1),
    "K2O": ("K", 2, 1), "P2O5": ("P", 2, 5), "H2O": ("H", 2, 1), "CO2": ("C", 1, 2),
}
VOLATILES = ("H2O", "CO2")


def molar_mass(oxide: str) -> float:
    cation, n, o = OXIDES[oxide]
    return n * ATOMIC_WEIGHT[cation] + o * ATOMIC_WEIGHT["O"]


MOLAR_MASS = {ox: molar_mass(ox) for ox in OXIDES}
FEO_PER_FE2O3 = 2 * MOLAR_MASS["FeO"] / MOLAR_MASS["Fe2O3"]     # 0.8998


@dataclass(frozen=True)
class Composition:
    """Oxide weight percent.  Missing oxides are zero.  Values need not sum to 100."""
    wt: dict[str, float] = field(default_factory=dict)
    name: str = ""
    source: str = ""

    def __post_init__(self):
        unknown = set(self.wt) - set(OXIDES)
        if unknown:
            raise ValueError(f"unknown oxide(s): {', '.join(sorted(unknown))}")
        bad = [ox for ox, v in self.wt.items() if v < 0]
        if bad:
            raise ValueError(f"negative amount for {', '.join(bad)}")

    # ---- basic views -----------------------------------------------------------------------
    def get(self, oxide: str) -> float:
        return self.wt.get(oxide, 0.0)

    @property
    def total(self) -> float:
        return sum(self.wt.values())

    def normalized(self, volatile_free: bool = False) -> "Composition":
        """Rescaled to 100 wt% (optionally without H2O and CO2)."""
        wt = {ox: v for ox, v in self.wt.items() if not (volatile_free and ox in VOLATILES)}
        s = sum(wt.values())
        if s <= 0:
            raise ValueError("composition is empty")
        return Composition({ox: 100.0 * v / s for ox, v in wt.items()}, self.name, self.source)

    def moles(self) -> dict[str, float]:
        """Moles of each oxide per 100 g of the composition as given."""
        return {ox: v / MOLAR_MASS[ox] for ox, v in self.wt.items() if v > 0}

    def mole_fractions(self, volatile_free: bool = False, total_iron_as_feo: bool = False
                       ) -> dict[str, float]:
        n = self.moles()
        if volatile_free:
            n = {ox: v for ox, v in n.items() if ox not in VOLATILES}
        if total_iron_as_feo:
            n = dict(n)
            n["FeO"] = n.get("FeO", 0.0) + 2.0 * n.pop("Fe2O3", 0.0)
        s = sum(n.values())
        return {ox: v / s for ox, v in n.items()}

    def cation_fractions(self, volatile_free: bool = True) -> dict[str, float]:
        """Fractions on a single-cation basis (SiO2, AlO1.5, ...); Fe3+ kept apart as 'Fe3'."""
        c: dict[str, float] = {}
        for ox, mol in self.moles().items():
            if volatile_free and ox in VOLATILES:
                continue
            cation, k, _ = OXIDES[ox]
            key = "Fe3" if ox == "Fe2O3" else cation
            c[key] = c.get(key, 0.0) + k * mol
        s = sum(c.values())
        return {k: v / s for k, v in c.items()}

    # ---- iron --------------------------------------------------------------------------------
    @property
    def feo_total(self) -> float:
        """Total iron expressed as FeO, wt%."""
        return self.get("FeO") + FEO_PER_FE2O3 * self.get("Fe2O3")

    @property
    def fe3_ratio(self) -> float:
        """Fe3+ / total Fe (molar)."""
        n = self.moles()
        fe3 = 2.0 * n.get("Fe2O3", 0.0)
        tot = fe3 + n.get("FeO", 0.0)
        return fe3 / tot if tot > 0 else 0.0

    def with_fe3_ratio(self, ratio: float) -> "Composition":
        """Same total iron, split so that Fe3+/total Fe = ratio."""
        if not 0.0 <= ratio <= 1.0:
            raise ValueError("Fe3+/total Fe must be between 0 and 1")
        feo_t = self.feo_total
        wt = dict(self.wt)
        wt["FeO"] = feo_t * (1.0 - ratio)
        wt["Fe2O3"] = feo_t * ratio / FEO_PER_FE2O3
        return Composition(wt, self.name, self.source)

    @property
    def mg_number(self) -> float:
        """100 Mg / (Mg + Fe2+), molar."""
        n = self.moles()
        mg, fe2 = n.get("MgO", 0.0), n.get("FeO", 0.0)
        return 100.0 * mg / (mg + fe2) if mg + fe2 > 0 else 0.0

    def table(self) -> str:
        rows = [f"{ox:<6} {self.get(ox):8.3f}" for ox in OXIDES if self.get(ox) > 0]
        return "\n".join(rows + [f"{'Total':<6} {self.total:8.3f}"])


# Starting compositions.  Each carries its source; check values against the paper before
# relying on them for validation.
PRESETS: dict[str, Composition] = {
    "MORB": Composition(
        {"SiO2": 50.47, "TiO2": 1.68, "Al2O3": 14.70, "FeO": 10.43, "MnO": 0.184, "MgO": 7.58,
         "CaO": 11.39, "Na2O": 2.79, "K2O": 0.160, "P2O5": 0.184},
        name="ALL MORB (global mean)",
        source="Gale et al. (2013) G-cubed 14:489-518, Table 1, ALL MORB mean; FeO is total Fe",
    ),
    "MELTS example basalt": Composition(
        {"SiO2": 48.68, "TiO2": 1.01, "Al2O3": 17.64, "Fe2O3": 0.89, "Cr2O3": 0.0425, "FeO": 7.59,
         "MgO": 9.10, "CaO": 12.45, "Na2O": 2.65, "K2O": 0.03, "P2O5": 0.08, "H2O": 0.2},
        name="MELTS example basalt",
        source="ThermoEngine documentation, MELTS v1.0.2 fractionation example (high-alumina basalt)",
    ),
    "KLB-1 peridotite": Composition(
        {"SiO2": 44.48, "TiO2": 0.16, "Al2O3": 3.59, "Cr2O3": 0.31, "FeO": 8.10, "MnO": 0.12,
         "MgO": 39.22, "CaO": 3.44, "Na2O": 0.30, "K2O": 0.02, "NiO": 0.25},
        name="KLB-1 spinel lherzolite",
        source="Davis, Tangeman and Hirschmann (2009) Am. Mineral. 94:176-180, Table 1 (KLB-1 of "
               "Takahashi 1986); FeO is total Fe. Entered from memory: check against the paper",
    ),
}
