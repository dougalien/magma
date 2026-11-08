# Magma

A Dougalien Systems app ([dougalien.com](https://dougalien.com)).  Open magma modeling in Python:
cool a silicate melt and follow what crystallizes (equilibrium or fractional), or heat a rock and
follow the melt that forms.  The high-temperature companion to Reaction.

Magma finds the assemblage of melt and minerals with the lowest total Gibbs energy at a given
temperature, pressure, bulk composition and oxygen fugacity, in the manner of MELTS (Ghiorso and
Sack 1995).  It is an original implementation from the published equations and parameters; the
open-source ThermoEngine is used only to check numbers.  See [docs/DESIGN.md](docs/DESIGN.md).

## On the web (Streamlit)

Magma also runs as a hosted Streamlit app with the same two windows: Input (the app's address) and
Results (`/results`, opened from the "Open results window" button, so it can sit on a second screen).
The two windows share one run through a private `?w=` id in the address; each visitor has their own.
The science code is the same; only the front end (`src/magma/web/`) is new.

    pip install -r requirements.txt
    streamlit run streamlit_app.py

On Streamlit Community Cloud: New app, this repository, entry point `streamlit_app.py`.

## Status

Milestones 1 and 2 are under way.  Working now, for anhydrous to slightly hydrous basaltic melts:

- `magma.composition`, `magma.redox`: oxide compositions (MORB preset from Gale et al. 2013),
  fO2 buffers (Frost 1991) and melt Fe3+/Fe2+ (Kress and Carmichael 1991).
- `magma.thermo.liquid`, `magma.thermo.standard_state`: the Ghiorso and Sack (1995) liquid,
  complete except for the standard state of H2O (fine while water stays in the melt).
- `magma.thermo.minerals`, `magma.thermo.phases`: olivine (Sack and Ghiorso 1989, Hirschmann
  1991) and ternary feldspar (Elkins and Grove 1990) on MELTS end-member data.
- `magma.thermo.pyroxene`: clinopyroxene and orthopyroxene, the seven-end-member model of Sack and Ghiorso (1994c)
  with Fe-Mg and Fe3+-Al ordering, our own code for the model's equations (its Taylor coefficients
  are in `data/melts_cpx_taylor.tsv` and `data/melts_opx_taylor.tsv`, evaluated from ThermoEngine's parameters with Douglas's
  approval on 2026-11-01); checked term by term against the definitions in ThermoEngine.
- `magma.thermo.spinel`: spinel (chromite, hercynite, magnetite, spinel, ulvospinel; Sack and
  Ghiorso 1991 as in MELTS) with Mg, Al and Fe3+ ordering between the tetrahedral and octahedral
  sites; same approach as the pyroxenes (`data/melts_spinel_taylor.tsv`).  It covers Cr-spinel and
  titanomagnetite.  End members whose oxide is absent (chromite without Cr2O3) are switched off.
- `magma.thermo.rhombohedral`: the ilmenite-hematite solution of Ghiorso and Evans (2008) as in
  MELTS (geikielite, hematite, ilmenite, corundum; pyrophanite is left out because the liquid has no
  MnO), with Fe, Mg-Ti ordering between the layers; its G matches ThermoEngine's definition exactly.
- Pigeonite: a second copy of the clinopyroxene model (as MELTS has Augite and Pigeonite), searched
  below 0.30 Ca per formula; left out above 15 kbar.  MORB at 1 bar now grows augite (Wo37) and
  pigeonite (Wo14) together from 1150 C, with pigeonite replacing olivine.
- Garnet: almandine-grossular-pyrope, Berman (1990) as in MELTS GarnetBerman.m.  KLB-1 at 3 GPa is a
  garnet lherzolite (garnet Py71 Alm12 Grs17) that begins to melt near 1500 C.
- Water: the dissolved-H2O standard state of MELTS v1.0.2 (LiquidMeltsH2O.m) and a pure-water fluid
  from IAPWS-95 (`thermo/water.py`, copied from Douglas's Reaction app).  MORB at 1200 C holds about
  3 wt% H2O at 1 kbar; with 2.5 wt% H2O it saturates near 1040 C, and plagioclase waits until 1020 C.
- `magma.equilibrium`: liquid + solids at one T, P and fO2 (or closed to oxygen), adding and
  dropping phases by their affinities, as MELTS does.
- `magma.path`: equilibrium or fractional crystallization along a cooling path.
- `magma.differentiation`, `magma.classify`: mantle to magma (melting, then fractional
  crystallization), with magma-type and TAS names; see below.
- `magma.melting`: partial melting of a rock on heating, batch (equilibrium; computed by cooling the
  rock's composition from above its liquidus and reported from low to high T) or fractional (melt
  above 1% of the residue removed at each step and pooled).

```python
from magma.composition import PRESETS
from magma.path import crystallize

for r in crystallize(PRESETS["MORB"], 1200, 1150, step=10, fo2="QFM", mode="fractional"):
    print(r.summary())
```

Data: liquid and mineral parameters are in `src/magma/data/*.tsv`, each with its source.  Values
that are hard to read in the scanned papers were taken from ThermoEngine (the open-source MELTS
code, AGPL-3), numbers only, with Douglas's approval on 2026-11-01.

Checks against MELTS (ThermoEngine's MELTS v1.0.2 MORB example, 1220 C, 100 MPa):

| Quantity | MELTS | Magma |
| --- | --- | --- |
| Forsterite affinity of the example liquid | 3495 J | 3491 J |
| Plagioclase on that liquid | present, An81 | 146 J short of saturation, An82 |
| Olivine on that liquid | 119 J undersaturated | 277 J undersaturated |
| Augite on that liquid | 1334 J undersaturated | 1394 J undersaturated |
| Orthopyroxene on that liquid | 2295 J undersaturated | 2397 J undersaturated |
| Spinel on that liquid | saturated, Mg 0.79 Cr 0.45 Al 1.32 per 4 O | 160 J oversaturated, Mg 0.79 Cr 0.44 |
| Equilibrium of the bulk composition | plagioclase 1.47 g An81 + 0.02 g spinel | plagioclase 0.43 g An82 + 0.03 g spinel |

ALL MORB at 1 bar and QFM: plagioclase appears between 1190 and 1200 C, then olivine (Fo~83,
olivine-melt Kd 0.32) and augite (Wo~40, Mg# 85) together near 1180 C; by 1150 C the system is
56% melt.  Olivine runs out near 1130 C; ilmenite (Ilm69 Gk23 Hem2) joins at 1110 C and the melt
fraction then drops quickly (25% to 7% by 1100 C), with olivine back as Fo~50.  The path ends at
1090 C with 3.5% melt: Magma cannot follow a system past its solidus yet.  Without ilmenite,
titanomagnetite (about 90% ulvospinel) takes its place at 1110 C.  Each mineral appears only once
(no augite + pigeonite pair).

The cooled MELTS example basalt (QFM, 1 kbar) grows Cr-spinel from 1245 C, then plagioclase at 1220 C,
olivine at 1210 C and augite at 1200 C.

KLB-1 peridotite at 1 GPa and QFM (batch): the rock is 1.7% melt at 1300 C, cpx is used up near
1390 C at about 20% melt, then the residue is harzburgite (olivine + opx + Cr-spinel), 25% melt at
1450 C; Hirose and Kushiro (1993) find cpx out at a similar degree of melting.  Fractional melting of
the same rock exhausts cpx near 1415 C, after which melting nearly stops (about 21 g of pooled melt
by 1450 C).  The lowest-degree melts are far too rich in Na2O (9-11 wt%), the same failing as MELTS
(the reason pMELTS exists).  Garnet is not in the model, so keep peridotite below about 2 GPa.
The KLB-1 composition was entered from memory and should be checked against Davis et al. (2009).

## From mantle to magma: Bowen's series

`magma.differentiation` (and "Mantle to magma" in the app) runs the textbook story in two stages.
Stage 1 batch-melts KLB-1 at the source (default 1380 C, 10 kbar: about 10% melt, a picritic basalt
with 46 wt% SiO2 and 14 wt% MgO).  Stage 2 takes that melt to the crust (default 2 kbar), adds water
(default 2 wt%, standing for a wet source above a subduction zone; KLB-1 itself is dry), and cools it
at the NNO buffer with fractional crystallization.  `magma.classify` names the melt at every step:
magma type by SiO2 (basalt, basaltic andesite, andesite, dacite, rhyolite) and its TAS field
(Le Bas et al. 1986).

```python
from magma.differentiation import melt_source, iter_differentiate, magma_type_stages
src = melt_source()                                   # KLB-1, 1380 C, 10 kbar
steps = list(iter_differentiate(src.primary, 1380, 900, step=10, p_crust=2000, h2o=2.0, fo2="NNO"))
for s in magma_type_stages(steps):
    print(s.t, s.magma_type, round(s.melt_grams, 1), s.tas)
```

With the defaults the minerals come in Bowen's order: olivine from about 1290 C, augite with Cr-spinel
near 1170 C, titanomagnetite at 1130 C, plagioclase near 1050 C (An falling as it cools) and
orthopyroxene near 930 C.  The melt is a basalt until about 1030 C (37 g left of 100 g), a basaltic
andesite to 950 C, and an andesite by silica (about 60 wt% SiO2) at 900 C with 20 g of melt left.
Dry and at QFM the same melt instead follows the iron-enrichment (Fenner, tholeiitic) trend: plagioclase
comes in at 1220 C, no magnetite forms, and the melt gets poorer in silica.

Where the model runs out: the MELTS liquid holds sodium too well, which ThermoEngine's own MELTS does
too (its MORB fractionation example reaches 7 wt% Na2O at 1025 C; Magma's run of the same example
tracks its plotted melt to within about 2 wt% for each oxide down to 1075 C).  The mantle melt is poor in potassium, so its evolved melt
builds up Na2O: biotite comes in at 820 C at about 62 wt% SiO2, and below about 800 C the melt's silica
falls while Na2O + K2O climbs past 12 wt%.  The app flags that tail and says to read no further.  The
example stops at 800 C by default.

## The late minerals and the granite minimum

Quartz, tridymite, K-feldspar (a second copy of the feldspar model, Or above 30%, as MELTS has
Plagioclase and Sanidine), muscovite, biotite (annite-phlogopite) and hornblende (pargasite,
ferropargasite, magnesiohastingsite) are in `thermo/felsic.py`, with the MELTS v1.0.2 (rhyolite-MELTS)
parameters from ThermoEngine (numbers only; Douglas's approval).  For them, water in the melt and in the
fluid is now on the minerals' energy scale: G_water(T, 1 bar) from IAPWS-95, shifted to H = -285830 J and
S = 69.95 J/K at 298.15 K.  Checked against ThermoEngine's MELTS v1.0.2 example (Early Bishop Tuff,
760 C, 175 MPa):

| Phase | MELTS affinity | Magma |
| --- | --- | --- |
| Biotite | 19819 J undersaturated | 20257 J |
| Muscovite | 25188 J | 25310 J |
| Hornblende | 54565 J | 55573 J |
| Tridymite | 2141 J | 2198 J |
| Quartz, sanidine | present (1.2 g, 2.3 g Or49) | 57 J and 114 J short of saturation |

The "Arc andesite" card (bulk continental crust with 3 wt% H2O, 2 kbar, NNO, fractional) runs the end
of Bowen's series: orthopyroxene and augite, then plagioclase with magnetite, dacite by 990 C, rhyolite
by 870 C, and quartz with biotite at 730 C and K-feldspar at 720 C, where the last melt (77 wt% SiO2)
freezes at the granite minimum.  MELTS's hornblende is a simple model and seldom saturates; it does not
form in these runs.

## Cards

The control window starts with cards: melts to cool (MORB, high-alumina basalt, arc andesite, dacite,
Bishop Tuff rhyolite), rocks to melt (KLB-1 peridotite, amphibolite, lower crust, upper crust) and the
mantle-to-magma story.  Each sets the composition and conditions; `composition.CARDS` lists them.
The crustal compositions (Rudnick and Gao 2003) were entered from memory and need checking against the
paper, like KLB-1.

Two fixes came with this: a spinel whose most stable composition is rich in a scarce oxide (Cr-spinel
from a melt with a few hundred ppm Cr2O3) could not be added even when oversaturated, because the
solver's first trial amount needed more Cr2O3 than the system had; and the spinel model now has a
second, chromite-free copy ("titanomagnetite"), as MELTS allows two spinels, since the Cr-spinel copy
stays on its Cr-rich side of the solvus.  Without Cr2O3 in the system only one copy is used.

Next: speed (each step now takes 5-15 s with sixteen phases), and a better-behaved liquid for sodium.

## A second engine: Holland & Powell

MELTS is weak where Bowen's series ends: its hornblende never forms, and its liquid keeps too much
sodium below about 800 C.  Magma is gaining a second, self-consistent system for cool, wet magmas:
the Holland & Powell end-member dataset with the igneous solution models (Holland and Powell 2011;
Green et al. 2016; Holland, Green and Powell 2018).  Numbers are taken from MAGEMin's distribution
of the dataset, with Douglas's approval (2026-11-02); the code is our own.  The plan has five steps:
dataset, solid solutions, pyroxenes and amphibole, liquid and solver, then the app's engine switch.

Step 1 is in: `data/hp_ds636.tsv` (291 end members) and `magma/hp/endmembers.py`, which gives G(T, P)
from the Holland and Powell (2011) heat capacity, modified Tait equation of state and thermal pressure,
Landau and Bragg-Williams order-disorder, the liquid end-member form, and Pitzer and Sterner (1994) for
H2O and CO2.  Checks (tests/test_hp_endmembers.py):

| check | Magma (HP engine) | reference |
|---|---|---|
| albite = jadeite + quartz, 600 C | 16.5 kbar | 16.2 kbar (Holland 1980) |
| quartz = coesite, 1000 C | 30.4 kbar | about 30.5 kbar (Bose and Ganguly 1995) |
| Al2SiO5 triple point | 538 C, 4.3 kbar | about 550 C, 4.5 kbar (Pattison 1992) |
| H2O, G(P) - G(1 bar), 900 C, 3 kbar | 77.15 kJ | 77.02 kJ (IAPWS-95) |
| quartz alpha-beta | Cp spike at 847 K | 847 K |

Against MAGEMin itself (built and run only as a reference, with Douglas's approval): all 259 end members
with an equation of state agree to within 0.001 J at five (T, P) points from 300 C and 1 kbar to 1500 C and
25 kbar.  A table of 59 of them at three points is kept in tests/data/hp_em_magemin.tsv.

Steps 2 and 3 are in: `magma/hp/solutions.py`, one framework for the site-fraction models (ideal
activities as products of site fractions, symmetric and asymmetric van Laar excess terms, DQF end members)
and thirteen igneous-set solids: olivine, feldspar, spinel, ilmenite, garnet, biotite, white mica,
orthopyroxene, clinopyroxene, amphibole (11 end members, Green et al. 2016), cordierite, epidote and
chlorite.  At seventeen MAGEMin equilibria (basalt, wet basalt, KLB-1, granite, andesite and pelite,
2-20 kbar, 550-1350 C; amphibole at seven of them) every end member's
chemical potential matches MAGEMin's oxide chemical potentials to 0.0001 kJ (tests/test_hp_solutions.py,
points in tests/data/hp_magemin_points.json).  The fluid moves to step 4, since it uses the melt's
formulation.

Step 4 is in: the silicate liquid and aqueous fluid of Holland, Green and Powell (2018), and
`magma/hp/minimize.py`, a Gibbs energy minimizer of our own (linear programming with column generation,
then a Newton solve of the exact equilibrium on the assemblage found).  At the seventeen MAGEMin points
it finds the same assemblage (two feldspars and two amphiboles included) with oxide chemical potentials
within 0.001 kJ (tests/test_hp_minimize.py).  The system can be open to oxygen at a buffer (an O2
reservoir at the buffer's fO2), and a path step starts from the step before: if the assemblage holds, a
step takes about a second instead of 5-30.

Step 5 is in: `magma/hp/engine.py` runs cooling and melting paths and the mantle-to-magma example on this
engine, and the control window has an engine switch (MELTS or Holland & Powell).  The HP phases are shown
under the app's mineral names; a feldspar with Or above 30% is K-feldspar, a clinopyroxene with Wo below
25% is pigeonite, and a spinel that is mostly magnetite and ulvospinel is titanomagnetite.  Phases MELTS
lacks keep their own names (cordierite, epidote, chlorite, rutile, sphene, the Al2SiO5 minerals,
corundum, nepheline).  MnO, P2O5, NiO, CoO and CO2 are left out.  The engine also works below the solidus,
so a cooling path runs until the melt is used up, and batch melting starts from the solid rock.

## The app

`python -m magma.app` (or `Magma.exe`) starts a small local server and opens two browser windows,
like Reaction: a control window (starting melt, temperatures, pressure, oxygen buffer, equilibrium or
fractional, minerals allowed) and a results window (plot of crystals and melt, melt composition and
mineral compositions against temperature, a step table, CSV and SVG downloads). It listens on
127.0.0.1 port 8766 only.

### Building Magma.exe

On Windows, from the Magma folder: `powershell -ExecutionPolicy Bypass -File packaging\build_exe.ps1`.
The exe lands in `dist\Magma.exe` (one file, about 60 MB, with Python, numpy and scipy inside).
Once Magma is on GitHub, `.github/workflows/windows-exe.yml` builds it on every push to main and
publishes it as the "Latest Windows build" release, the same way Reaction does.

## Tests

```
python -m pytest
```
