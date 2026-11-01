# Magma

A Dougalien Systems app ([dougalien.com](https://dougalien.com)).  Open magma modelling in Python:
cool a silicate melt and follow what crystallizes (equilibrium or fractional), or heat a rock and
follow the melt that forms.  The high-temperature companion to Reaction.

Magma finds the assemblage of melt and minerals with the lowest total Gibbs energy at a given
temperature, pressure, bulk composition and oxygen fugacity, in the manner of MELTS (Ghiorso and
Sack 1995).  It is an original implementation from the published equations and parameters; the
open-source ThermoEngine is used only to check numbers.  See [docs/DESIGN.md](docs/DESIGN.md).

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
- `magma.equilibrium`: liquid + solids at one T, P and fO2 (or closed to oxygen), adding and
  dropping phases by their affinities, as MELTS does.
- `magma.path`: equilibrium or fractional crystallization along a cooling path.

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
56% melt.  Olivine runs out near 1130 C; titanomagnetite (about 90% ulvospinel) joins at 1110 C and
the melt fraction then drops quickly (18% to 7% by 1100 C), with olivine back as Fo~60.  Ilmenite is
still missing, and each mineral appears only once (no augite + pigeonite pair).  Below about 5% melt
the results are not reliable.

The cooled MELTS example basalt (QFM, 1 kbar) grows Cr-spinel from 1245 C, then plagioclase at 1220 C,
olivine at 1210 C and augite at 1200 C.

Next: the Fe-Ti oxides (ilmenite-hematite), a second pyroxene (pigeonite), then partial melting of a rock.

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
