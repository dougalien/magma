# Magma user manual

Magma 0.1.0. For students and instructors who want to cool a magma and watch it crystallize, or heat
a rock and watch it melt. This manual covers the web app (Streamlit); the desktop app has the same two
windows and the same controls.

Contents

1. [What Magma does](#1-what-magma-does)
2. [Starting Magma](#2-starting-magma)
3. [The Input window](#3-the-input-window)
4. [Running and stopping](#4-running-and-stopping)
5. [The Results window](#5-the-results-window)
6. [The two engines](#6-the-two-engines)
7. [Starting compositions](#7-starting-compositions)
8. [Starting cards](#8-starting-cards)
9. [Limits](#9-limits)
10. [How Magma is checked](#10-how-magma-is-checked)
11. [How to cite](#11-how-to-cite)
12. [Terms used on screen](#12-terms-used-on-screen)

## 1. What Magma does

At each temperature, pressure, bulk composition and oxygen fugacity, Magma finds the mix of melt and
minerals with the lowest total Gibbs energy: the assemblage at equilibrium. It does this one
temperature step at a time along a path you choose:

- **Cool a melt**: start above the liquidus and cool. Crystals either stay with the melt (equilibrium
  crystallization) or are taken out at each step (fractional crystallization).
- **Heat a rock**: start near or below the solidus and heat. The melt either stays with the rock
  (batch melting) or is drained off at each step (fractional melting).
- **Mantle to magma**: melt a mantle rock at depth, then cool that melt in the crust with fractional
  crystallization.

Every number on screen is computed by Magma's own solver from published thermodynamic models. Magma
does not use empirical mineral-melt partitioning models.

## 2. Starting Magma

**On the web.** Open the app's address. The Input window opens first. Press **Open results window** to
open the Results window in a second browser window or tab; put it on a second screen if you have one.
The two windows share one run through a private `?w=` code in the address, so each visitor has their
own run. Copy the address with its `?w=` code to come back to the same run; a run left untouched for
12 hours is dropped.

**On your own computer.**

    pip install -r requirements.txt
    streamlit run streamlit_app.py

or, without Streamlit, `python -m magma.app`, which starts a small local server (127.0.0.1, port 8766)
and opens the Input window (called the control window there) in your browser.

**From Python.**

```python
from magma.composition import PRESETS
from magma.path import crystallize

for r in crystallize(PRESETS["MORB"], 1200, 1150, step=10, fo2="QFM", mode="fractional"):
    print(r.summary())
```

## 3. The Input window

### Starting cards

The top panel holds the starting cards: five melts to cool, four rocks to melt, and **Mantle to magma:
Bowen's series**. A card fills in a composition, the process, the mode, the temperatures, the pressure
and the oxygen buffer. It does not change the engine or the list of minerals allowed. You can change
anything a card filled in before you press **Run**. Section 8 lists every card.

### Starting melt or rock

- **Composition**: a list of starting compositions (section 7). Picking one fills in the oxide boxes
  and shows its source below the list.
- **Oxide boxes**: SiO₂, TiO₂, Al₂O₃, Fe₂O₃, Cr₂O₃, FeO, MnO, MgO, NiO, CoO, CaO, Na₂O, K₂O, P₂O₅
  and H₂O, in wt%. Type your own analysis here. The total is shown; Magma normalizes to 100 g. Give
  iron as FeO and Fe₂O₃, or all of it as FeO. A composition must contain SiO₂, and no value may be
  negative.
- Water (H₂O) is part of the composition. There is no separate water control except in Mantle to
  magma.

### Process

- **Engine**: *Holland & Powell* (the default) or *MELTS*. Section 6 explains the difference.
- **Process**: *Cool a melt*, *Heat a rock* or *Mantle to magma*.
- **Mode** (not shown for Mantle to magma):

  | Process | Equilibrium | Fractional |
  | --- | --- | --- |
  | Cool a melt | The crystals stay with the melt and keep reacting with it (a closed system). | The crystals are taken out at each step; the next step starts from the remaining melt. |
  | Heat a rock | Batch melting: the melt stays with the rock and keeps reacting with it. | The melt above 1% of the rock is taken out at each step and pooled; the residue keeps melting. |

- **Minerals allowed**: the minerals the solver may form. All are allowed at the start. Leave one out
  to see what the magma does without it, for example a basalt with no augite. The list is: olivine,
  plagioclase, augite, pigeonite, orthopyroxene, spinel, titanomagnetite, ilmenite, garnet,
  hornblende, biotite, K-feldspar, muscovite, quartz, tridymite and water (fluid). With the Holland &
  Powell engine some minerals are always allowed (section 6).

### Conditions (Cool a melt, Heat a rock)

| Control | Meaning | Allowed values |
| --- | --- | --- |
| Start, °C | First temperature of the path | 600 to 2000 °C |
| End, °C | Last temperature | 600 to 2000 °C; below Start when cooling, above Start when heating |
| Step, °C | Temperature step | positive, at most 400 steps in the run |
| Pressure, bar | Constant along the path (1 kbar = 1000 bar) | 1 to 30,000 bar |
| Oxygen | QFM, NNO, IW or MH buffer, or *Closed to oxygen* | |
| Offset, log units | Shift from the buffer, for example +1 for one log unit above NNO | any number |

With a buffer, the system is open to oxygen: it takes up or gives off oxygen so that the melt's
Fe³⁺/Fe²⁺ follows the buffer as temperature changes. *Closed to oxygen* keeps the Fe₂O₃ the
composition started with. The buffers are those of Frost (1991): QFM (quartz-fayalite-magnetite),
NNO (nickel-nickel oxide), IW (iron-wüstite) and MH (magnetite-hematite).

### Mantle to magma

Two stages. First the mantle rock partly melts at the source (batch melting); that melt is the primary
magma. Then the melt rises to the crust and cools, and the crystals are taken out as they form
(fractional crystallization). The card loads KLB-1 peridotite with these defaults:

| Control | Default | Allowed values |
| --- | --- | --- |
| 1. Melting temperature, °C | 1380 | 1100 to 1700 °C |
| 1. Pressure, bar | 10,000 | 1 to 30,000 bar |
| 2. Pressure, bar | 2000 | 1 to 30,000 bar |
| 2. Water added, wt% | 2 | 0 to 10 wt% |
| 2. Cool to, °C | 800 | 600 °C up to the melting temperature |
| 2. Step, °C | 10 | positive, at most 400 steps |
| 2. Oxygen | NNO | NNO, QFM or NNO+1 |

Stage 1 is always run at the QFM buffer. The mantle rock is dry; the water added in stage 2 stands for
a wet source, as above a subduction zone. Stage 2 masses are grams per 100 g of primary melt. You can
pick another composition for stage 1, but the conditions above were chosen for a mantle rock.

## 4. Running and stopping

Press **Run**. A progress bar shows the step being computed; both windows follow the run live. Press
**Stop** to end the run after the current step; the steps already done are kept. When the run ends the
bar reports *done*, *stopped* or the error.

How long a run takes depends on the engine, the composition and the number of steps. The Holland &
Powell engine is usually much faster than the MELTS engine (section 8 gives the times of the starting
cards).
In Mantle to magma no step appears until stage 1 is done; meanwhile the window shows "Melting the
mantle ... (this takes a minute or two)".

Messages you may see at the end of a run:

- *Ended at ... °C: the melt is used up (the solidus).* (Holland & Powell) or *the melt is nearly gone
  (near the solidus), below which the MELTS engine cannot go* (MELTS). The path stopped before the End
  temperature because the melt ran out.
- *Melting starts near ... °C (the solidus).* Heating started below the solidus; the results begin at
  the first step with melt.
- *No melt at these conditions* or *The rock does not melt at these temperatures*: the whole path is
  below the solidus.
- An error in red: an input outside the limits above (the message says which), or a step the solver
  could not complete. Change the conditions and run again, or try the other engine.

## 5. The Results window

**Title and conditions.** The title names the composition and the process; the line below gives the
temperatures, step, pressure, oxygen buffer, engine and run number.

**Magma types.** For cooling runs and Mantle to magma, a chain such as
`basalt · 1180 °C → basaltic andesite · 1100 °C` marks the first step at which the melt reaches each
new magma type.

**What we found.** A short written account computed from the run: the first crystals and their
formula, the order in which minerals appeared, how much melt is left and its SiO₂ and Mg#. For melting
runs it gives the melt fraction at the first and last step, which minerals were used up and at what
temperature, and for fractional melting the amount and composition of the pooled melt.

**Views.** Pick one with the buttons above the plot.

| View | What it shows |
| --- | --- |
| Masses | Grams of melt and of each mineral per 100 g of starting material, against temperature. In fractional crystallization the crystal curves are the totals taken out so far. In fractional melting a dashed line shows the total melt removed. |
| Melt composition | Oxide wt% of the melt against temperature (H₂O not plotted). |
| Mineral compositions | Olivine Fo, plagioclase An, K-feldspar Or, augite Mg# and Wo, pigeonite and orthopyroxene Mg#, spinel Cr# and ulvöspinel, ilmenite FeTiO₃, hornblende and biotite Mg#, and the melt's Mg#, all in mol%, against temperature. |
| TAS diagram | The melt on the total alkali-silica diagram (Le Bas et al. 1986), SiO₂ against Na₂O + K₂O, volatile-free with all iron counted. Larger dots mark each new magma type; in Mantle to magma the primary melt is marked. |
| Bowen's series | The minerals of Bowen's reaction series that formed in this run and the temperature each first appeared. Gray boxes are minerals that did not form, or that the model does not include. Plagioclase is shown as Ca-rich (its first An) and Na-rich (its last An, if below An50). Spinel, titanomagnetite, ilmenite, garnet and water are listed below the boxes. |

Plots can be zoomed and hovered; the camera icon above a plot saves it as an image.

In fractional melting the mineral masses are what is left in the residue, which keeps 1% melt. With an
oxygen buffer the total mass can be a little over 100 g: the melt takes up oxygen as some of its iron
oxidizes. In fractional crystallization any water fluid is taken out with the crystals at each step.

**Every step.** A table with one row per temperature step: melt %, melt grams, magma type, TAS name,
grams of each mineral and the melt's oxides and Mg#. **Download the table (CSV)** saves the full table,
which also has the pressure and each mineral's structural formula.

**Notes on the models.** The data, models and limits of the engine used in this run.

**Sources.** The **Sources and how to cite** link at the foot of each window lists every paper, data set
and program behind Magma, with DOIs.

## 6. The two engines

Both engines minimize Gibbs energy; they differ in the thermodynamic models.

| | Holland & Powell (default) | MELTS |
| --- | --- | --- |
| Data | End members of Holland and Powell (2011), data set ds636 | Berman (1988) end members with the MELTS adjustments |
| Liquid | Holland, Green and Powell (2018) | Ghiorso and Sack (1995), as in MELTS v1.0.2 |
| Minerals | Olivine, feldspar, spinel, ilmenite, garnet, orthopyroxene, clinopyroxene (Holland et al. 2018); amphibole (Green et al. 2016); biotite, white mica, cordierite, chlorite (White et al. 2014); epidote; quartz and the other silica minerals, the Al₂SiO₅ minerals, rutile, sphene, corundum, nepheline | Olivine (Sack and Ghiorso 1989; Hirschmann 1991); plagioclase and K-feldspar (Elkins and Grove 1990); augite, pigeonite, orthopyroxene (Sack and Ghiorso 1994); spinel and titanomagnetite (Sack and Ghiorso 1991); ilmenite (Ghiorso and Evans 2008); garnet (Berman 1990); quartz, tridymite, muscovite, biotite, hornblende (MELTS v1.0.2) |
| Water | Aqueous fluid of Holland et al. (2018); H₂O from Pitzer and Sterner (1994) | Dissolved H₂O as in MELTS v1.0.2; pure-water fluid from IAPWS-95 (Wagner and Pruss 2002) |
| Oxides left out | MnO, P₂O₅, NiO, CoO | MnO, NiO, CoO |
| Below the solidus | Melting paths can start below the solidus and report from the first melt; cooling paths stop when the melt is gone | Cooling paths stop when the melt is nearly gone; melting starts at the solidus |
| Checked against | MAGEMin 2.1.0 (Riel et al. 2022) | ThermoEngine, MELTS v1.0.2 |

**How the Holland & Powell engine names minerals.** One model can stand for two minerals on screen: a
feldspar with more than 30% Or is shown as K-feldspar, a clinopyroxene with less than 25% Wo as
pigeonite, and a spinel with more than half magnetite and ulvöspinel as titanomagnetite. Cordierite,
epidote, chlorite, the silica and Al₂SiO₅ minerals, rutile, sphene, corundum and nepheline are always
allowed with this engine; the Minerals allowed list does not switch them off.

**How the MELTS engine names pyroxenes.** Augite and pigeonite are two copies of one clinopyroxene
model, as in MELTS. The name on screen follows the composition: 0.30 Ca or more per six oxygens is
augite. Pigeonite is left out above 15 kbar. K-feldspar is a second copy of the feldspar model, and
titanomagnetite a second copy of the spinel model.

**Iron and oxygen.** With the MELTS engine the melt's Fe³⁺/Fe²⁺ comes from Kress and Carmichael (1991)
at the chosen buffer. With either engine, *Closed to oxygen* keeps the starting Fe₂O₃.

**Which engine to use.** The Holland & Powell engine finds all phases at once, has the metamorphic
minerals (cordierite, garnet, micas, amphibole) and works from below the solidus, so it suits melting
of crustal rocks and wet magmas. The MELTS engine is the model behind the MELTS program, so its results
can be compared with MELTS. The engines differ, sometimes a great deal (section 8): comparing them on
the same card is itself a useful exercise.

## 7. Starting compositions

All values in wt%. "FeO" is total iron as FeO where the source gives total iron.

| Name in the list | Source | Notes |
| --- | --- | --- |
| MORB | Gale et al. (2013), Table 1, ALL MORB mean | |
| MELTS example basalt | ThermoEngine MELTS v1.0.2 fractionation example (high-alumina basalt) | 0.2 wt% H₂O, as there |
| KLB-1 peridotite | Davis et al. (2009), Table 1, recommended composition of KLB-1 (Takahashi 1986) | dry |
| Bishop Tuff rhyolite | Early Bishop Tuff average melt inclusion, ThermoEngine MELTS v1.0.2 equilibrium example (rhyolite-MELTS, Gualda et al. 2012) | 5.5 wt% H₂O, as there |
| Andesite (bulk crust) | Rudnick and Gao (2003), bulk continental crust | 3 wt% H₂O added by Magma |
| Dacite (upper crust) | Rudnick and Gao (2003), upper continental crust | 4 wt% H₂O added by Magma |
| Granite source (upper crust) | Rudnick and Gao (2003), upper continental crust | 1.5 wt% H₂O added by Magma |
| Lower crust | Rudnick and Gao (2003), lower continental crust | 0.5 wt% H₂O added by Magma |
| Amphibolite (wet MORB) | Gale et al. (2013) ALL MORB | 2 wt% H₂O added by Magma |

The water added to the crustal compositions is Magma's choice, not part of the published analyses.

## 8. Starting cards

Each card's settings, from `magma.composition.CARDS`. "Fractional" for a rock means fractional melting;
"equilibrium" means batch melting.

| Card | Composition | Process, mode | Temperatures, step | Pressure | Oxygen |
| --- | --- | --- | --- | --- | --- |
| Mid-ocean ridge basalt | MORB | Cool a melt, fractional | 1250 to 1050 °C, 10 °C | 1000 bar | QFM |
| High-alumina basalt | MELTS example basalt | Cool a melt, fractional | 1300 to 1100 °C, 10 °C | 1000 bar | QFM |
| Arc andesite | Andesite (bulk crust) | Cool a melt, fractional | 1150 to 700 °C, 10 °C | 2000 bar | NNO |
| Dacite | Dacite (upper crust) | Cool a melt, fractional | 1050 to 680 °C, 10 °C | 2000 bar | NNO |
| Rhyolite (Bishop Tuff) | Bishop Tuff rhyolite | Cool a melt, equilibrium | 800 to 700 °C, 5 °C | 1750 bar | closed |
| Mantle peridotite | KLB-1 peridotite | Heat a rock, batch | 1250 to 1500 °C, 10 °C | 10,000 bar | QFM |
| Amphibolite | Amphibolite (wet MORB) | Heat a rock, batch | 850 to 1200 °C, 10 °C | 10,000 bar | NNO |
| Lower crust | Lower crust | Heat a rock, batch | 950 to 1300 °C, 10 °C | 8000 bar | QFM |
| Upper crust | Granite source (upper crust) | Heat a rock, batch | 700 to 1000 °C, 10 °C | 5000 bar | NNO |
| Mantle to magma | KLB-1 peridotite | stage 1 batch melting, stage 2 fractional crystallization | stage 1 at 1380 °C; stage 2 from 1380 to 800 °C, 10 °C, 2 wt% water added | 10,000 bar, then 2000 bar | QFM, then NNO |

What the cards give in Magma 0.1.0, run as the app runs them, with all minerals allowed:

Melts to cool. "First crystals" is the first step with a mineral; the order lists each mineral at the
first step it appears; the last column gives the melt left at the End temperature, its SiO₂ (wt%, as
shown in *What we found*) and its magma type.

| Card | Engine | First crystals | Then | At the end |
| --- | --- | --- | --- | --- |
| Mid-ocean ridge basalt | H&P | plagioclase, 1180 °C | olivine 1170, augite 1150, ilmenite 1080, titanomagnetite 1070 °C | 21.1 g melt, 58.3% SiO₂, andesite |
| | MELTS | plagioclase, 1200 °C | augite 1190, pigeonite 1180, ilmenite 1110, spinel 1100 °C | 12.5 g melt, 59.3% SiO₂, andesite |
| High-alumina basalt | H&P | plagioclase, 1210 °C | olivine 1200, augite 1160 °C | 41.5 g melt, 49.4% SiO₂, basalt |
| | MELTS | Cr-spinel, 1240 °C | plagioclase 1220, olivine 1210, augite 1200, titanomagnetite and pigeonite 1100 °C | 25.8 g melt, 49.9% SiO₂, basalt |
| Arc andesite | H&P | plagioclase, 1070 °C | orthopyroxene 1010, augite 960, titanomagnetite 910, ilmenite 900, hornblende 800, biotite and water 760 °C | 36.4 g melt, 74.4% SiO₂, rhyolite |
| | MELTS | orthopyroxene, 1130 °C | augite 1050, spinel 1020, plagioclase 990, ilmenite 860, water 810, quartz and biotite 730, K-feldspar 720 °C | 6.1 g melt, 65.0% SiO₂, dacite |
| Dacite | H&P | plagioclase, 980 °C | orthopyroxene 940, titanomagnetite 880, ilmenite 860, biotite 770, water 740, augite 690 °C | 54.2 g melt, 75.9% SiO₂, rhyolite |
| | MELTS | orthopyroxene, 1050 °C | spinel 990, augite 920, plagioclase 910, ilmenite 850, water 790, biotite, quartz and K-feldspar 740, garnet 680 °C | 7.0 g melt, 63.7% SiO₂, dacite |
| Rhyolite (Bishop Tuff), equilibrium | H&P | fayalite-rich olivine and titanomagnetite, 735 °C | plagioclase 725, biotite and water 715, K-feldspar 710 °C | 85.4% melt, 75.5% SiO₂ |
| | MELTS | water, 800 °C (saturated from the start) | K-feldspar, quartz and plagioclase 755, spinel, orthopyroxene, olivine and ilmenite 750, augite 715 °C | 0.7% melt |

Rocks to melt (batch melting). The melt fraction at the first and last step, the minerals left in the
residue at the first step, and the minerals used up on heating.

| Card | Engine | First step | Used up | Last step |
| --- | --- | --- | --- | --- |
| Mantle peridotite | H&P | 1250 °C, 3.1% melt; olivine, orthopyroxene, augite, spinel | augite 1340, orthopyroxene 1490 °C | 1500 °C, 35.8% melt, 51.7% SiO₂, 19.1% MgO |
| | MELTS | 1250 °C, 1.1% melt; olivine, spinel, augite, orthopyroxene | orthopyroxene 1400 °C | 1500 °C, 27.0% melt, 47.0% SiO₂, 19.6% MgO |
| Amphibolite | H&P | 850 °C, 17.0% melt; hornblende, augite, plagioclase, ilmenite, quartz | quartz 860, ilmenite 980, garnet 1020, hornblende 1050, orthopyroxene 1070, plagioclase 1120 °C | 1200 °C, 99.2% melt |
| | MELTS | 850 °C, 23.9% melt; augite, garnet, plagioclase | plagioclase 900, orthopyroxene 1100, garnet 1130, pigeonite 1160 °C | 1200 °C, 69.2% melt |
| Lower crust | H&P | 950 °C, 14.0% melt; plagioclase, orthopyroxene, augite, ilmenite | ilmenite 1020, augite 1190, orthopyroxene 1210, plagioclase 1230 °C | fully molten at 1230 °C |
| | MELTS | 950 °C, 14.1% melt; orthopyroxene, plagioclase, augite | augite 1110, plagioclase 1230 °C | 1300 °C, 98.1% melt |
| Upper crust | H&P | 700 °C, 11.0% melt; plagioclase, quartz, biotite, ilmenite | augite 760, biotite 810, K-feldspar 830, quartz 880, ilmenite 900 °C | 1000 °C, 85.3% melt |
| | MELTS | 700 °C, 9.5% melt; orthopyroxene, plagioclase, ilmenite, quartz, K-feldspar, biotite | biotite 760, K-feldspar 810, augite 870, ilmenite and quartz 920, spinel 970 °C | 1000 °C, 79.6% melt |

Some minerals are used up that were not in the residue at the first step: they formed on the way.

Mantle to magma, with the defaults:

| | H&P | MELTS |
| --- | --- | --- |
| Stage 1, 1380 °C and 10 kbar | 20.6% melt; residue olivine, orthopyroxene, spinel | 10.5% melt; residue olivine, spinel, augite, orthopyroxene |
| Primary melt | 48.6% SiO₂, 14.3% MgO, Mg# 76, basalt | 46.0% SiO₂, 13.9% MgO, Mg# 77, basalt |
| Minerals in stage 2 | olivine 1280, spinel 1270, augite 1150, plagioclase 1100, water 1020, titanomagnetite 1000, orthopyroxene 970, ilmenite 940 °C | olivine 1290, augite and spinel 1170, titanomagnetite 1130, plagioclase 1060, water 950, orthopyroxene 910, biotite 830 °C |
| Magma types | basaltic andesite 1000, andesite 970, dacite 930, rhyolite 880 °C | basaltic andesite 1030, andesite 940 °C |
| At 800 °C | 10.4 g melt, 74.6% SiO₂ and 7.8% Na₂O + K₂O (volatile-free), rhyolite | 10.6 g melt, 61.3% SiO₂ and 14.0% Na₂O + K₂O (volatile-free), andesite by silica, off the TAS diagram in alkalis |

The two engines often disagree, most of all for the wet, silica-rich melts: the MELTS engine grows
quartz in the Arc andesite, Dacite and Rhyolite cards, and the Holland & Powell engine grows none within
the cards' temperatures. The MELTS mantle-to-magma melt ends with 14.0 wt% Na₂O + K₂O at 61.3 wt%
SiO₂, off the TAS diagram: the sodium problem of section 9.

The runs took from about 20 seconds to 4 minutes each with the Holland & Powell engine and from about
3 to 26 minutes with the MELTS engine, four runs at a time on one server.

## 9. Limits

- **Sodium in the MELTS liquid.** The Ghiorso and Sack (1995) liquid holds sodium too well, as
  ThermoEngine's MELTS does too. Evolved melts build up Na₂O and can plot in the trachy- fields of the
  TAS diagram where real melts would be andesite or dacite. When the melt's silica falls while
  Na₂O + K₂O passes 10 wt%, *What we found* says so and advises reading the path no further.
- **Simple hydrous minerals in MELTS.** MELTS's hornblende and biotite are simple models and seldom
  saturate; treat them as rough.
- **Missing minerals.** Neither engine has apatite or zircon. The MELTS engine has no sphene,
  cordierite or Na-rich amphibole, and its garnet is the mantle (pyrope-almandine-grossular) garnet.
- **Missing components.** CO₂ is not modelled. MnO, NiO and CoO are left out of both engines, and P₂O₅
  of the Holland & Powell engine; they are dropped before the calculation.
- **Range of the models.** Holland, Green and Powell (2018) give their models for 0.001 to 70 kbar and
  from 650 °C to peridotite liquidus temperatures. The Holland & Powell engine considers no melt below
  500 °C. Magma's own input limits are 600 to 2000 °C and 1 bar to 30 kbar.
- **Constant pressure.** Each path is isobaric; decompression melting and ascent are not modelled
  except as the two separate pressures of Mantle to magma.
- **Water at 1 bar.** At 1 bar any water leaves the melt at once.
- **Equilibrium at every step.** Magma computes equilibrium; it has no kinetics, undercooling or
  crystal settling rates. Fractional crystallization removes all crystals at every step, and fractional
  melting keeps 1% melt in the residue.
- **Holland et al. (2018) corrigendum.** Green et al. (2025) correct an algebraic error in the Holland
  et al. (2018) melt model. Magma's liquid reproduces MAGEMin 2.1.0; whether that release includes the
  correction has not been checked.

## 10. How Magma is checked

The tests in `tests/` (run with `python -m pytest`) compare Magma with reference codes and data:

- MELTS engine: affinities and assemblages of the ThermoEngine MELTS v1.0.2 MORB and Early Bishop Tuff
  examples; the README tabulates the comparison.
- Holland & Powell engine: Gibbs energies of 59 end members at three temperature-pressure points, and
  seventeen equilibria (basalt, wet basalt, KLB-1, granite, andesite and pelite, 2 to 20 kbar,
  550 to 1350 °C), against MAGEMin 2.1.0; and three end-member reactions against experiments: albite =
  jadeite + quartz (Holland 1980), quartz = coesite (Bose and Ganguly 1995) and the Al₂SiO₅ triple point
  (Pattison 1992).
- Water: IAPWS-95 against the verification values of Wagner and Pruss (2002).

## 11. How to cite

If you use Magma in a study or a class, cite Magma itself and the models and data behind your result.

- **Magma**: the repository's `CITATION.cff` gives the reference (GitHub's "Cite this repository"
  button formats it). Give the version you used (shown at the foot of each window).
- **The models**: for the Holland & Powell engine, Holland and Powell (2011) and Holland, Green and
  Powell (2018), plus Green et al. (2016) or White et al. (2014) if amphibole, micas, cordierite or
  chlorite matter to your result. For the MELTS engine, Ghiorso and Sack (1995) and the mineral models
  that formed in your run (section 6).
- **The starting composition**: its source (section 7).

Full references with DOIs are on the app's Sources page and in `docs/REFERENCES.md`.

Magma is free software under the GNU General Public License, version 3 or later. It comes with no
warranty.

## 12. Terms used on screen

| Term | Meaning |
| --- | --- |
| Liquidus, solidus | The temperatures at which the first crystal appears on cooling, and the last melt freezes. |
| Equilibrium, fractional | Crystals (or melt) stay and react, or are removed at each step. |
| Batch melting | Equilibrium melting: the melt stays with the rock. |
| Pooled melt | All the melt removed so far in fractional melting, mixed. |
| Primary melt | The melt that leaves the mantle source in Mantle to magma. |
| Oxygen buffer | A mineral assemblage that fixes oxygen fugacity at a given temperature and pressure: QFM, NNO, IW, MH (Frost 1991). |
| Mg# | 100 Mg / (Mg + Fe²⁺), molar. For the melt, Fe²⁺ is the FeO after the iron split at the buffer. |
| Fo, An, Or, Wo | Mol% forsterite in olivine, anorthite and orthoclase in feldspar, wollastonite (Ca) in pyroxene. |
| Cr#, ulvöspinel | 100 Cr / (Cr + Al) in spinel; mol% Fe₂TiO₄ in spinel. |
| Magma type | Basalt (45 to 52 wt% SiO₂), basaltic andesite (52 to 57), andesite (57 to 63), dacite (63 to 69), rhyolite (above 69), volatile-free. Below 45 the melt is called ultramafic. |
| TAS name | The field of the total alkali-silica diagram (Le Bas et al. 1986) the melt plots in. |
