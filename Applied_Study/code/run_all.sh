#!/usr/bin/env bash
# Reproduce the whole applied study.  Run from the repository root:
#     bash Applied_Study/code/run_all.sh
# Keep one worker unless a measured memory budget supports more.
set -euo pipefail
export CBCF_ROOT="${CBCF_ROOT:-$(pwd)}"
export CBCF_CORES="${CBCF_CORES:-1}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CBCF_OUTDIR="${CBCF_OUTDIR:-$CBCF_ROOT/Applied_Study/analysis}"
export CBCF_MAX_MEMORY_GB="${CBCF_MAX_MEMORY_GB:-6}"
LOG="$CBCF_OUTDIR/logs"
mkdir -p "$LOG"

# Refuse a new job if its cap would consume the desktop/other-work reserve.
check_memory_budget() {
  python3 -c 'import os, pathlib, sys
p = pathlib.Path("/proc/meminfo")
if p.exists():
    available = next(int(line.split()[1]) for line in p.read_text().splitlines() if line.startswith("MemAvailable:")) / 1048576
    needed = float(os.environ["CBCF_MAX_MEMORY_GB"]) + 8
    if available < needed:
        sys.exit(f"Only {available:.1f} GiB available; need {needed:.1f} GiB including the 8 GiB desktop reserve.")'
}
check_memory_budget

# Preserve the executable source and base commit at execution time, rather
# than trying to reconstruct provenance from a later working tree.
python3 "$CBCF_ROOT/Applied_Study/code/export_reproducibility.py" --root "$CBCF_ROOT" --source-only --destination "$LOG/execution_source"

# Build the captured local source in an isolated library. Preserve R's normal
# library search paths for dependencies; common.R prepends CBCF_R_LIBRARY.
export CBCF_R_LIBRARY="${CBCF_R_LIBRARY:-$(mktemp -d /tmp/countbcf-library-XXXXXX)}"
mkdir -p "$CBCF_R_LIBRARY"
if [ "${CBCF_SKIP_PACKAGE_INSTALL:-0}" != "1" ]; then
  if command -v prlimit >/dev/null 2>&1; then
    MAKEFLAGS=-j1 prlimit --as="$((CBCF_MAX_MEMORY_GB * 1024 * 1024 * 1024))" -- R CMD INSTALL --no-multiarch --library="$CBCF_R_LIBRARY" "$CBCF_ROOT" > "$LOG/package_install.log" 2>&1
  else
    MAKEFLAGS=-j1 R CMD INSTALL --no-multiarch --library="$CBCF_R_LIBRARY" "$CBCF_ROOT" > "$LOG/package_install.log" 2>&1
  fi
fi

# Record the environment before any fitting. Fixed seeds alone do not
# identify historical package versions or execution environments.
Rscript -e '.libPaths(c(Sys.getenv("CBCF_R_LIBRARY"), .libPaths()))
library(countbcf)
writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("CBCF_OUTDIR"), "logs/sessionInfo.txt"))
ip <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
write.csv(ip[order(ip$Package), ], file.path(Sys.getenv("CBCF_OUTDIR"), "logs/package_versions.csv"), row.names = FALSE)
writeLines(c(find.package("countbcf"), as.character(packageVersion("countbcf"))), file.path(Sys.getenv("CBCF_OUTDIR"), "logs/loaded_package.txt"))'

for s in 01_eda_parametric 02_countbcf_main 03_sensitivity 04_overlap_cluster 06_attributable 07_lowdose_city 08_city_footing 09_stomach_sweep 10_profile_comparison 05_figures; do
  check_memory_budget
  echo ">>> $s"
  if command -v prlimit >/dev/null 2>&1; then
    prlimit --as="$((CBCF_MAX_MEMORY_GB * 1024 * 1024 * 1024))" -- Rscript "$CBCF_ROOT/Applied_Study/code/$s.R" 2>&1 | tee "$LOG/$s.log"
  else
    Rscript "$CBCF_ROOT/Applied_Study/code/$s.R" 2>&1 | tee "$LOG/$s.log"
  fi
done
python3 "$CBCF_ROOT/Applied_Study/code/export_reproducibility.py" --root "$CBCF_ROOT" --analysis "$CBCF_OUTDIR"
echo ">>> done. Outputs: $CBCF_OUTDIR"
echo ">>> environment: $LOG/sessionInfo.txt, package_versions.csv"
