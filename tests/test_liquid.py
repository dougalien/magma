import numpy as np
import pytest

from magma.composition import PRESETS, Composition
from magma.redox import at_fo2
from magma.thermo.liquid import (COMPONENTS, W, from_components, gibbs_mixing, mu_mixing,
                                 to_components)

MORB = at_fo2(PRESETS["MORB"], 1200, 1, "QFM")


def test_w_table_symmetric_and_spot_values():
    i = {c: k for k, c in enumerate(COMPONENTS)}
    assert np.allclose(W, W.T) and np.all(np.diag(W) == 0)
    # spot checks against Table 4 of Ghiorso and Sack (1995)
    assert W[i["SiO2"], i["Na2SiO3"]] == -99039
    assert W[i["Fe2O3"], i["Fe2SiO4"]] == -179065
    assert W[i["Al2O3"], i["Na2SiO3"]] == -130785
    assert W[i["Na2SiO3"], i["H2O"]] == -96938
    assert np.count_nonzero(np.triu(W)) == 65        # 66 pairs, H2O-MgCr2O4 blank


def test_component_conversion_round_trip():
    lc = to_components(MORB)
    back = from_components(lc.moles)
    for ox, mol in MORB.moles().items():
        if ox in back:
            assert back[ox] == pytest.approx(mol, rel=1e-12, abs=1e-15)
    assert "MnO" in lc.ignored
    assert np.all(lc.moles >= 0)                     # MORB lies inside the component space


def _numeric_mu(n, t, k, h=1e-7):
    f = lambda m: m.sum() * gibbs_mixing(m / m.sum(), t)
    up, dn = n.copy(), n.copy()
    up[k] += h
    dn[k] -= h
    return (f(up) - f(dn)) / (2 * h)


@pytest.mark.parametrize("water", [0.0, 3.0])
def test_chemical_potentials_are_derivatives_of_G(water):
    comp = Composition({**MORB.wt, "H2O": water})
    n = to_components(comp).moles
    t = 1473.15
    mu = mu_mixing(n / n.sum(), t)
    for k in range(len(COMPONENTS)):
        if n[k] > 0:
            assert mu[k] == pytest.approx(_numeric_mu(n, t, k), rel=1e-5, abs=1e-2)
    x = n / n.sum()
    pos = x > 0
    assert float(x[pos] @ mu[pos]) == pytest.approx(gibbs_mixing(x, t), rel=1e-10)
