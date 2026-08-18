#!/usr/bin/env bash
# Reproduce the whole applied study.  Run from the repository root:
#     bash Applied_Study/code/run_all.sh
# Optional: CBCF_CORES=8 bash Applied_Study/code/run_all.sh
set -euo pipefail
export CBCF_ROOT="${CBCF_ROOT:-$(pwd)}"
export CBCF_CORES="${CBCF_CORES:-6}"
LOG="$CBCF_ROOT/Applied_Study/analysis/logs"
mkdir -p "$LOG"

# Record the exact environment first: the seeds are fixed (common.R:44, and
# they flow through common.R:202, 02_countbcf_main.R:42 and
# 04_overlap_cluster.R:166), so a clean clone reproduces the pipeline, but
# reproducing the NUMBERS also needs the package versions.
Rscript -e 'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/analysis/logs/sessionInfo.txt"))
ip <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
write.csv(ip[order(ip$Package), ], file.path(Sys.getenv("CBCF_ROOT", "."), "Applied_Study/analysis/logs/package_versions.csv"), row.names = FALSE)'

for s in 01_eda_parametric 02_countbcf_main 03_sensitivity 04_overlap_cluster 06_attributable 07_lowdose_city 08_city_footing 09_stomach_sweep 05_figures; do
  echo ">>> $s"
  Rscript "$CBCF_ROOT/Applied_Study/code/$s.R" 2>&1 | tee "$LOG/$s.log"
done
echo ">>> done.  Tables: Applied_Study/analysis/tables  Plots: Applied_Study/analysis/plots"
echo ">>> environment: Applied_Study/analysis/logs/sessionInfo.txt, package_versions.csv"
