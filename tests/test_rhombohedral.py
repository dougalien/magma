"""Ilmenite-hematite solution (Ghiorso and Evans 2008 as in MELTS)."""
import numpy as np
import pytest

from magma.composition import PRESETS
from magma.path import crystallize
from magma.thermo import rhombohedral as rh

T, P = 1373.15, 1000.0


def test_solution_gibbs_matches_melts_definition():
    # RhombohedralBerman.m's H - T S macros evaluated at these r, s (1100 C, 1000 bar)
    z = [0.55, 0.12, 0.05, 0.03, 0.3, 0.05, 0.02]
    assert rh.g_and_grad(z, T, P)[0] == pytest.approx(-17429.8810919111, abs=1e-6)


def test_gradient_and_pure_end_member_reference():
    z = np.array([0.55, 0.12, 0.05, 0.03, 0.3, 0.05, 0.02])
    g, dg = rh.g_and_grad(z, T, P)
    num = [(rh.g_and_grad(z + 1e-6 * e, T, P)[0] - rh.g_and_grad(z - 1e-6 * e, T, P)[0]) / 2e-6 for e in np.eye(7)]
    assert dg == pytest.approx(num, abs=1e-4)
    # the solution at pure ilmenite, with its own ordering, is the ilmenite reference state
    assert rh.order([1, 0, 0, 0], T, P)[1] == pytest.approx(rh.end_member_g(T, P)[2], abs=1e-6)


def test_chemical_potentials_are_partial_derivatives():
    ph = rh.RHOMBOHEDRAL
    x = np.array([0.15, 0.1, 0.7, 0.05])
    mu = ph.mu(x, T, P)
    g = lambda n: n.sum() * ph.gibbs(n / n.sum(), T, P)
    num = [(g(x + 1e-6 * e) - g(x - 1e-6 * e)) / 2e-6 for e in np.eye(4)]
    assert mu == pytest.approx(num, abs=0.05)


def test_morb_grows_ilmenite_late():
    path = crystallize(PRESETS["MORB"], 1130, 1100, 10)
    last = next(s for s in path[-1].solids if s.name == "rhombohedral")
    assert 0.5 < last.x["ilmenite"] < 0.9 and last.x["hematite"] < 0.1
    assert not any(s.name == "rhombohedral" for s in path[0].solids)
