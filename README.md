# Magma

Magma models silicate melts and rocks: cool a melt and follow what crystallizes (equilibrium or
fractional crystallization), or heat a rock and follow the melt that forms (batch or fractional
melting). It is a Dougalien Systems app ([dougalien.com](https://dougalien.com)), free to use and
open source.

At each temperature, pressure, bulk composition and oxygen fugacity, Magma finds the assemblage of
melt and minerals with the lowest total Gibbs energy. It has two thermodynamic engines, both written
from the published equations and parameters:

- **MELTS**: the liquid of Ghiorso and Sack (1995) with the MELTS minerals, as in MELTS v1.0.2.
- **Holland and Powell**: the end-member dataset of Holland and Powell (2011) with the igneous
  solution models, liquid and aqueous fluid of Holland, Green and Powell (2018).

Every paper and data set is listed in [docs/REFERENCES.md](docs/REFERENCES.md) and on the app's
Sources page. The design is described in [docs/DESIGN.md](docs/DESIGN.md).

## Running Magma

On the web (Streamlit):

    pip install -r requirements.txt
    streamlit run streamlit_app.py

The app has two windows: Input (the app's address) and Results (`/results`, opened from the "Open
results window" button, so it can sit on a second screen). The two share one run through a private
`?w=` id in the address; each visitor has their own.

On a desktop, `python -m magma.app` starts a small local server (127.0.0.1, port 8766) and opens the
same two windows in the browser. On Windows, `packaging\build_exe.ps1` builds a single `Magma.exe`.

From Python:

```python
from magma.composition import PRESETS
from magma.path import crystallize

for r in crystallize(PRESETS["MORB"], 1200, 1150, step=10, fo2="QFM", mode="fractional"):
    print(r.summary())
```

## What is in the models

MELTS engine (`magma.thermo`, `magma.equilibrium`, `magma.path`, `magma.melting`):

- Melt: the Ghiorso and Sack (1995) liquid, with the dissolved-H2O standard state of MELTS v1.0.2;
  melt Fe3+/Fe2+ from Kress and Carmichael (1991); oxygen buffers from Frost (1991).
- Olivine (Sack and Ghiorso 1989; Hirschmann 1991); feldspar, as plagioclase and K-feldspar
  (Elkins and Grove 1990); augite, pigeonite and orthopyroxene (Sack and Ghiorso 1994); spinel, as
  Cr-spinel and titanomagnetite (Sack and Ghiorso 1991); ilmenite-hematite (Ghiorso and Evans 2008);
  garnet (Berman 1990); quartz, tridymite, muscovite, biotite and hornblende (MELTS v1.0.2); and a
  pure-water fluid from IAPWS-95.
- The engine cannot follow a system below its solidus; a cooling path ends when the melt is nearly
  gone.

Holland and Powell engine (`magma.hp`):

- Olivine, feldspar, spinel, ilmenite, garnet, biotite, white mica, orthopyroxene, clinopyroxene,
  amphibole (Green et al. 2016), cordierite, epidote and chlorite; the silicate liquid and aqueous
  fluid of Holland, Green and Powell (2018); H2O and CO2 from Pitzer and Sterner (1994).
- A Gibbs energy minimizer of Magma's own (linear programming with column generation, then a Newton
  solve of the equilibrium on the assemblage found). The system can be closed or held at an oxygen
  buffer. This engine works below the solidus.
- MnO, P2O5, NiO, CoO and CO2 are left out of this engine.

Parameter values that are hard to read in the scanned papers were taken from ThermoEngine (the
open-source MELTS code) and, for the Holland and Powell engine, from MAGEMin (Riel et al. 2022),
numbers only. The code is Magma's own.

`magma.differentiation` and `magma.classify` run the "mantle to magma" example: batch melting of
KLB-1 peridotite at the source, then fractional crystallization of that melt in the crust, with each
melt named by its silica content and its total alkali-silica field (Le Bas et al. 1986).

## Starting cards

The app starts with cards: melts to cool (mid-ocean ridge basalt, high-alumina basalt, arc andesite,
dacite, Bishop Tuff rhyolite), rocks to melt (KLB-1 peridotite, amphibolite, lower crust, upper
crust) and the mantle-to-magma example. Each sets a composition and conditions; `composition.CARDS`
lists them and `composition.PRESETS` gives each composition with its source.

## Checks

Against ThermoEngine's MELTS v1.0.2 MORB example liquid at 1220 C and 100 MPa (affinities; positive
means undersaturated):

| Quantity | MELTS | Magma |
| --- | --- | --- |
| Forsterite | 3495 J | 3491 J |
| Olivine | 119 J | 277 J |
| Augite | 1334 J | 1394 J |
| Orthopyroxene | 2295 J | 2397 J |
| Plagioclase | present, An81 | 146 J short of saturation, An82 |
| Spinel | saturated, Mg 0.79 Cr 0.45 per 4 O | 160 J oversaturated, Mg 0.79 Cr 0.46 |
| Equilibrium of the bulk composition | plagioclase 1.47 g An81 + 0.02 g spinel | plagioclase 0.43 g An82 + 0.03 g spinel |

Against ThermoEngine's MELTS v1.0.2 Early Bishop Tuff example (760 C, 175 MPa):

| Phase | MELTS affinity | Magma |
| --- | --- | --- |
| Biotite | 19819 J | 20257 J |
| Muscovite | 25188 J | 25310 J |
| Hornblende | 54565 J | 55573 J |
| Tridymite | 2141 J | 2198 J |
| Quartz, sanidine | present | 57 J and 114 J short of saturation |

Holland and Powell end members (`tests/test_hp_endmembers.py`):

| Check | Magma | Reference |
| --- | --- | --- |
| albite = jadeite + quartz, 600 C | 16.5 kbar | 16.25 kbar (Holland 1980) |
| quartz = coesite, 1000 C | 30.4 kbar | 30.7 kbar (Bose and Ganguly 1995) |
| Al2SiO5 triple point | 538 C, 4.3 kbar | about 550 C, 4.5 kbar (Pattison 1992) |
| quartz alpha-beta, 1 bar | heat-capacity spike at 847 K | 847 K |

Against MAGEMin 2.1.0 (igneous database, dataset 636): the Gibbs energies of 59 end members at three
(T, P) points agree to 0.001 J, and at seventeen equilibria (basalt, wet basalt, KLB-1, granite,
andesite and pelite, 2-20 kbar, 550-1350 C) the solver finds the same assemblage with oxide chemical
potentials within 0.001 kJ (`tests/data/`).

## Known limits

- The MELTS liquid holds sodium too well at low temperature, as ThermoEngine's MELTS does too: below
  about 800 C an evolved melt builds up Na2O, and the app flags that part of a path. The
  lowest-degree mantle melts are also too rich in Na2O, the reason pMELTS exists.
- MELTS's hornblende is a simple model and seldom saturates.

## Tests

    python -m pytest

## How to cite

If you use Magma in a study or a class, please cite it, and cite the models and data it used for
your result (the app's Sources page lists them). GitHub's "Cite this repository" button gives the
reference from `CITATION.cff`.

## License and credits

Magma is free software under the GNU General Public License, version 3 or later (`LICENSE`). Anyone
may use, copy, change and share it at no cost; copies and changed versions that are passed on must
stay under the same license. It comes with no warranty.

Magma's parameter values come from the published papers listed in docs/REFERENCES.md, read where
needed from ThermoEngine (ENKI project, AGPL-3) and MAGEMin (Riel et al. 2022, GPL-3). Its IAPWS-95
water code is shared with Reaction, another Dougalien Systems app.
