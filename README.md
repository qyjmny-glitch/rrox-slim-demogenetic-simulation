# Rhinopithecus roxellana SLiM extinction simulation

This repository contains scripts and input files for nonWF SLiM simulations of demographic and genetic extinction processes in Rhinopithecus roxellana.

## Project structure

- `slim/`: SLiM simulation scripts.
- `scripts/`: Python and shell scripts for preparing inputs, running simulations, and summarizing outputs.
- `config/`: parameter tables and scenario settings.
- `input/`: small input tables used by SLiM.
- `docs/`: workflow notes and parameter documentation.
- `tests/`: smoke tests and small test inputs.
- `results/`: simulation outputs, ignored by Git.
- `logs/`: runtime logs, ignored by Git.

## Main workflow

1. Prepare yearly demographic and environmental inputs.
2. Run burn-in simulations.
3. Run forward extinction scenarios from 5200 BP to present.
4. Summarize demographic, genetic, and extinction outputs.

## Notes

Large files such as BAM, VCF, tree-sequence files, full simulation outputs, and logs are not tracked in GitHub.
