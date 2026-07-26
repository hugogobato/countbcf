#!/usr/bin/env bash
# Reproduce the whole applied study.  Run from the repository root:
#     bash Applied_Study_v2/code/run_all.sh
# Optional: CBCF_CORES=8 bash Applied_Study_v2/code/run_all.sh
set -euo pipefail
export CBCF_ROOT="${CBCF_ROOT:-$(pwd)}"
export CBCF_CORES="${CBCF_CORES:-6}"
LOG="$CBCF_ROOT/Applied_Study_v2/analysis/logs"
mkdir -p "$LOG"

for s in 01_eda_parametric 02_countbcf_main 03_sensitivity 04_overlap_cluster 06_attributable 07_lowdose_city 05_figures; do
  echo ">>> $s"
  Rscript "$CBCF_ROOT/Applied_Study_v2/code/$s.R" 2>&1 | tee "$LOG/$s.log"
done
echo ">>> done.  Tables: Applied_Study_v2/analysis/tables  Plots: Applied_Study_v2/analysis/plots"
