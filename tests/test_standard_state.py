import math
import os
from pathlib import Path

import pytest

from magma.thermo.standard_state import END_MEMBERS, _glass_h_s



def test_silica_glass_heat_capacity():
    cp = lambda t: 127.3 - 10.777e-3 * t + 4.3127e5 / t ** 2 - 1463.8 / math.sqrt(t)
    assert cp(298.15) == pytest.approx(44.4, abs=0.5)        # measured silica glass
    assert cp(1480.0) < 81.373                               # glass below the liquid
    h, s = _glass_h_s(1000.0)                                # integrals agree with quadrature
    n = 20000
    ts = [298.15 + (1000 - 298.15) * (k + 0.5) / n for k in range(n)]
    dt = (1000 - 298.15) / n
    assert h == pytest.approx(sum(cp(t) for t in ts) * dt, rel=1e-6)
    assert s == pytest.approx(sum(cp(t) / t for t in ts) * dt, rel=1e-6)


def test_liquid_equals_reference_solid_at_fusion():
    from magma.thermo.standard_state import liquid_gibbs
    for name, e in END_MEMBERS.items():
        h, s_ = e.solid_h_s(e.t_fusion)
        assert liquid_gibbs(name, e.t_fusion, 1.0) == pytest.approx(h - e.t_fusion * s_, abs=1e-6)
        # liquid is the stable phase above fusion, solid below (constant-Cp check near Tf)
        for dt, sign in ((50, -1), (-50, 1)):
            t = e.t_fusion + dt
            hs, ss_ = e.solid_h_s(t)
            assert sign * (liquid_gibbs(name, t, 1.0) - (hs - t * ss_)) > 0, name


def test_morb_saturates_in_olivine_with_experimental_kd():
    from magma.composition import PRESETS
    from magma.redox import at_fo2
    from magma.saturation import saturation_temperature
    sat = saturation_temperature("olivine", PRESETS["MORB"])
    assert 1100 < sat.t_c < 1250
    melt = at_fo2(PRESETS["MORB"], sat.t_c, 1.0, "QFM").moles()
    kd = (sat.x["fayalite"] / sat.x["forsterite"]) / (melt["FeO"] / melt["MgO"])
    assert 0.28 < kd < 0.36                    # experiments on MORB give about 0.30
