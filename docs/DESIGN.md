# Magma: design notes

Agreed 2026-10-31: MELTS-style Gibbs energy minimization (not empirical mineral-melt models).
The full plan is the doc "Magma: architecture and scope"
(https://claude.ai/artifact/CT3EFvpLiowxUxAfx9JXr5).  Summary:

* Own code in Python (numpy, scipy), laid out like Reaction.  Never copied from MELTS.
* Liquid: Ghiorso and Sack (1995); pMELTS (Ghiorso et al. 2002) later for mantle melting.
* Solids: Berman (1988) end-members with the MELTS adjustments + solid-solution models.
* Iron redox: Kress and Carmichael (1991).  Buffers: Frost (1991).
* Checks: ThermoEngine (ENKI) and experimental liquid lines of descent.

## Milestones

1. Thermodynamic layer (liquid, Fe redox, olivine, plagioclase, clinopyroxene, spinel)
2. Equilibrium at one T, P (minimizer, liquidus finder)
3. Paths: equilibrium/fractional cooling, batch/fractional melting
4. Two-window app and Windows exe, as in Reaction
5. H2O-CO2, more phases, pMELTS; hand exsolved fluid to Reaction

## Layout

```
src/magma/
  composition.py   oxides, moles, fractions, presets
  redox.py         fO2 buffers, melt Fe3+/Fe2+
  data/            parameter tables (.tsv) with citations
  thermo/          liquid, standard states, minerals, solid solutions
  equilibrium.py   phase equilibrium (affinities + equilibrium conditions)
  path.py          cooling paths (melting coming)
  app/             local two-window app                       (coming)
```
