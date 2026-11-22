# Magma: design

Magma computes phase equilibria of silicate melts and rocks by Gibbs energy minimization, with
thermodynamic models taken from the literature and implemented from their published equations.
It does not use empirical mineral-melt partitioning models.

## Engines

* MELTS (`magma.thermo`, `magma.equilibrium`): the liquid of Ghiorso and Sack (1995), Berman (1988)
  end members with the MELTS adjustments, and the MELTS solid-solution models. Phases are added and
  dropped by their affinities (Ghiorso 1994). Iron redox follows Kress and Carmichael (1991); oxygen
  buffers follow Frost (1991).
* Holland and Powell (`magma.hp`): the ds636 end-member dataset (Holland and Powell 2011) and the
  igneous solution models (Holland, Green and Powell 2018; Green et al. 2016), with a Gibbs energy
  minimizer of Magma's own.

Both are checked against reference codes in the tests: ThermoEngine (MELTS v1.0.2) for the MELTS
engine and MAGEMin for the Holland and Powell engine.

## Layout

```
src/magma/
  composition.py     oxides, moles, fractions, preset compositions and starting cards
  redox.py           fO2 buffers, melt Fe3+/Fe2+
  data/              parameter tables (.tsv), each with its source
  thermo/            MELTS liquid, standard states, minerals and solid solutions, water
  equilibrium.py     MELTS phase equilibrium at one T, P and fO2
  path.py            cooling paths (equilibrium or fractional crystallization)
  melting.py         batch and fractional melting
  differentiation.py, classify.py   mantle to magma, magma-type and TAS names
  hp/                Holland and Powell end members, solutions, minimizer and engine adapter
  app/               local two-window app (standard-library server and static pages)
  web/               Streamlit front end with the same two windows
streamlit_app.py     Streamlit entry point
scripts/make_sources.py   writes docs/REFERENCES.md and the app's Sources page
```
