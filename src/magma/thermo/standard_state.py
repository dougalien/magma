"""Standard-state Gibbs energies of the liquid components (Ghiorso and Sack 1995, appendix).

For a component with a crystalline reference phase, the solid's enthalpy and entropy are carried
to its fusion temperature at 1 bar (Berman 1988 heat capacity, with a lambda term
Cp = T (l1 + l2 T)^2 below Tt and a first-order step deltaH at Tt), fused (Tf dSf), then carried
to T with the constant liquid heat capacity:

    H_l(T) = H_s(Tf) + Tf dSf + Cp_l (T - Tf)
    S_l(T) = S_s(Tf) + dSf + Cp_l ln(T / Tf)
    G_l(T, P) = H_l - T S_l + [V + dV/dT (T - 1673.15)](P - Pr)
               + [dV/dP + d2V/dTdP (T - 1673.15)](P^2/2 - Pr P + Pr^2/2)
               + d2V/dP2 (P^3/6 - Pr P^2/2 + Pr^2 P/2 - Pr^3/6)

Energies are Berman-Brown apparent Gibbs energies (H - T S, absolute S), J/mol.  Parameters are in
data/melts_liquid_endmembers.tsv (values as in ThermoEngine, which follow Table A1).

SiO2 is referenced to amorphous silica (footnote q of Table A1, Richet et al. 1982):
H = -901554 J, S = 48.475 J/K at 298.15 K; below 1480 K (the glass transition)
Cp = 127.3 - 10.777e-3 T + 4.3127e5 / T^2 - 1463.8 / sqrt(T); above it Cp = 81.373.
The paper prints the last coefficient as 463.8; that gives Cp = 102 J/K at 298 K and 99.5 J/K at
1480 K (above the liquid), while 1463.8 gives 44.2 J/K at 298 K, matching measured silica glass
(44.4 J/K).  We take 1463.8, which is also the value in ThermoEngine (LiquidMeltsSiO2.m).

H2O: the standard state MELTS v1.0.2 uses (ThermoEngine LiquidMeltsH2O.m, numbers taken with
Douglas's approval): G = RT (a/T + b + phi(T, P)) + G_water(T, 1 bar) - G_Robie(T), with a = -33676,
b = 18.3527, phi a quartic in P (the dissolved-water volume) and G_Robie a fit in T.  Magma measures
every water Gibbs energy from pure water at the same T and 1 bar, so liquid_gibbs("H2O") returns
RT (a/T + b + phi) - G_Robie, and the water fluid (thermo/water_fluid.py) returns
G_water(T, P) - G_water(T, 1 bar) from IAPWS-95.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from importlib import resources


TR, PR, T_VOL = 298.15, 1.0, 1673.15
h2o_standard_state_done = True


@dataclass(frozen=True)
class LiquidEndMember:
    name: str
    h: float
    s: float
    k: tuple[float, float, float, float]
    l1: float
    l2: float
    tt: float
    dh_t: float
    v: float
    dvdt: float
    dvdp: float
    d2vdpdt: float
    d2vdp2: float
    t_fusion: float
    ds_fusion: float
    cp: float

    def solid_h_s(self, t: float) -> tuple[float, float]:
        """Reference solid H and S at t and 1 bar."""
        k0, k1, k2, k3 = self.k
        h = (self.h + k0 * (t - TR) + 2 * k1 * (math.sqrt(t) - math.sqrt(TR)) - k2 * (1 / t - 1 / TR)
             - 0.5 * k3 * (1 / t ** 2 - 1 / TR ** 2))
        s = (self.s + k0 * math.log(t / TR) - 2 * k1 * (1 / math.sqrt(t) - 1 / math.sqrt(TR))
             - 0.5 * k2 * (1 / t ** 2 - 1 / TR ** 2) - k3 / 3 * (1 / t ** 3 - 1 / TR ** 3))
        if self.tt > 0:
            tl, l1, l2 = min(t, self.tt), self.l1, self.l2
            h += 0.5 * l1 * l1 * (tl ** 2 - TR ** 2) + 2 / 3 * l1 * l2 * (tl ** 3 - TR ** 3) \
                + 0.25 * l2 * l2 * (tl ** 4 - TR ** 4)
            s += l1 * l1 * (tl - TR) + l1 * l2 * (tl ** 2 - TR ** 2) + l2 * l2 / 3 * (tl ** 3 - TR ** 3)
            if t > self.tt:
                h += self.dh_t
                s += self.dh_t / self.tt
        return h, s


def _load_table() -> dict[str, LiquidEndMember]:
    text = resources.files("magma.data").joinpath("melts_liquid_endmembers.tsv").read_text("utf-8")
    rows = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    out = {}
    for r in csv.DictReader(rows, delimiter="\t"):
        f = lambda k: float(r[k])
        out[r["component"]] = LiquidEndMember(
            r["component"], f("H"), f("S"), (f("k0"), f("k1"), f("k2"), f("k3")), f("l1"), f("l2"),
            f("Tt"), f("deltaH"), f("vLiq"), f("dvdtLiq"), f("dvdpLiq"), f("d2vdtdpLiq"),
            f("d2vdp2Liq"), f("tFusion"), f("sFusion"), f("cpLiq"))
    return out


END_MEMBERS = _load_table()
# the SiO2 liquid's volume terms (Table A1); its energy is handled separately
_SIO2_V = (2.690, 0.0, -1.89e-5, 1.3e-8, 3.6e-10)


def _volume_term(v, dvdt, dvdp, d2vdpdt, d2vdp2, t: float, p: float) -> float:
    dt = t - T_VOL
    return ((v + dvdt * dt) * (p - PR) + (dvdp + d2vdpdt * dt) * (p * p / 2 - PR * p + PR * PR / 2)
            + d2vdp2 * (p ** 3 / 6 - PR * p * p / 2 + PR * PR * p / 2 - PR ** 3 / 6))


# ---- amorphous SiO2 ---------------------------------------------------------------------------
_SIO2_H0, _SIO2_S0, _T_GLASS, _SIO2_CP_LIQ = -901554.0, 48.475, 1480.0, 81.373


def _glass_h_s(t: float) -> tuple[float, float]:
    """Integrals of the glass Cp (and Cp/T) from 298.15 K to t."""
    a, b, c, d = 127.3, -10.777e-3, 4.3127e5, -1463.8
    h = (a * (t - TR) + b / 2 * (t * t - TR * TR) - c * (1 / t - 1 / TR)
         + 2 * d * (math.sqrt(t) - math.sqrt(TR)))
    s = (a * math.log(t / TR) + b * (t - TR) - c / 2 * (1 / t ** 2 - 1 / TR ** 2)
         - 2 * d * (1 / math.sqrt(t) - 1 / math.sqrt(TR)))
    return h, s


def _sio2_gibbs(t: float, p: float) -> float:
    if t <= _T_GLASS:
        dh, ds = _glass_h_s(t)
    else:
        dh, ds = _glass_h_s(_T_GLASS)
        dh += _SIO2_CP_LIQ * (t - _T_GLASS)
        ds += _SIO2_CP_LIQ * math.log(t / _T_GLASS)
    return (_SIO2_H0 + dh) - t * (_SIO2_S0 + ds) + _volume_term(*_SIO2_V, t, p)


def liquid_gibbs(component: str, t_k: float, p_bar: float) -> float:
    """mu0 of a liquid component, J/mol."""
    if component == "SiO2":
        return _sio2_gibbs(t_k, p_bar)
    if component == "H2O":
        return h2o_liquid_gibbs(t_k, p_bar)
    e = END_MEMBERS[component]
    h_s, s_s = e.solid_h_s(e.t_fusion)
    h = h_s + e.t_fusion * e.ds_fusion + e.cp * (t_k - e.t_fusion)
    s = s_s + e.ds_fusion + e.cp * math.log(t_k / e.t_fusion)
    return h - t_k * s + _volume_term(e.v, e.dvdt, e.dvdp, e.d2vdpdt, e.d2vdp2, t_k, p_bar)


def h2o_liquid_gibbs(t: float, p: float) -> float:
    """Dissolved-H2O standard state relative to pure water at T and 1 bar, J/mol (see module notes)."""
    r, a, b = 8.3143, -33676.0, 18.3527
    phi = ((0.110 / t + 4.432e-5 + 1.405e-7 * t - 2.394e-11 * t * t) * p
           + (7.337e-8 / t - 1.170e-8 - 9.502e-13 * t) * p * p
           + (1.876e-10 / t + 4.586e-13) * p ** 3 - 1.191e-14 * p ** 4 / t)
    g_robie = r * t * (2.9147 * math.log(t) - 9.6863e-4 * t + 6.8593e-8 * t * t + 77.8899 / math.sqrt(t)
                       - 28954.8 / t - 2263.27 / (t * t) - 15.8997)
    return r * t * (a / t + b + phi) - g_robie
