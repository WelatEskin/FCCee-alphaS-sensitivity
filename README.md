# FCC-ee Generator-Level Study of Event-Shape and Jet Observables

This repository contains the analysis codes used for a generator-level study of event-shape and jet observables for strong-coupling sensitivity at an FCC-ee-like electron-positron collider environment.

The study uses PYTHIA 8 for event generation, FastJet for jet reconstruction, ROOT for histogramming and analysis, and Python for plotting and event-display visualization.

## Physics process

The simulated process is:

e+ e- -> Z/gamma* -> q qbar -> hadrons

at a center-of-mass energy of:

sqrt(s) = 91.2 GeV

## Main observables

The main observables studied in this analysis are:

- thrust and 1 - T
- Durham jet multiplicity
- Durham leading jet energy
- Durham y23
- algorithm-dependent jet multiplicities

The jet algorithms considered are:

- Durham
- anti-kT-like
- kT-like
- Cambridge/Aachen-like

## Strong coupling scan

Template samples were generated for the following input values:

alphaS = 0.110, 0.114, 0.118, 0.122, 0.126, 0.130

An independent pseudo-data sample was generated at:

alphaS,true = 0.118

## Closure test

The closure test uses the 1 - T distribution in the fit region:

0.02 < 1 - T < 0.25

The number of bins used in the fit region is 23.

The parabolic interpolation around the chi-square minimum gives:

alphaS,fit = 0.11794 +/- 0.00014

This result should be interpreted as a generator-level closure-test result, not as a standalone experimental or theoretical precision determination of alphaS.

## Repository structure

```text
src/      C++ source codes for event generation and event-display preparation
scripts/  Shell scripts for building and running the analysis
macros/   ROOT macros for plotting and closure testing
python/   Python scripts for checks and event-display visualization
figures/  Example output figures
logs/     Example terminal outputs and closure-test logs
docs/     Additional notes
