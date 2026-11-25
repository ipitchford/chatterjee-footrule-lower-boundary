#!/bin/sh
# Regenerate every result file and the figure. Run from the package root.
# Uses python -O to confirm that no check depends on assert statements.
set -e
PY="${PYTHON:-python3}"
$PY -I -B -O code/verify_boundary.py > results/verification.csv
$PY -I -B -O code/finite_modes.py > results/finite-modes-check.csv
$PY -I -B -O code/two_predictors.py > results/two-predictors.txt
$PY -I -B -O code/boundary_table.py --output results/table.json
$PY -I -B -O code/sd_probe.py --output results/sd_probe.json
$PY -E -s -B code/plot_boundary.py
echo "run_all: done"
