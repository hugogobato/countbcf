#!/usr/bin/env python3
"""Export a whitelist of final replication artifacts, without private papers/notes."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile
from datetime import datetime, timezone


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--analysis", type=Path, help="Applied output directory")
    parser.add_argument("--destination", type=Path, help="Public replication directory")
    parser.add_argument("--source-only", action="store_true", help="Capture execution source before a run")
    args = parser.parse_args()
    root = args.root.resolve()
    analysis = (args.analysis or root / "Applied_Study/analysis").resolve()
    dest = (args.destination or root / "reproducibility").resolve()
    dest.mkdir(parents=True, exist_ok=True)

    def copy(source: Path, target: Path) -> None:
        if source.exists() and not args.source_only:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)

    # These are the final data/model objects and compact posterior diagnostics,
    # not internal review notes, correspondence, drafts, or preliminary runs.
    for name in ("main_cells.rds", "preston_fit.rds", "city_statistic_draws.rds",
                 "profile_age30_70_draws.rds"):
        copy(analysis / name, dest / "applied/models" / name)
    for pattern in ("summary_*.csv",):
        for source in (root / "tests/analysis").glob(pattern):
            copy(source, dest / "simulations" / source.name)
    copy(root / "tests/analysis/replicate_completeness.csv", dest / "simulations/replicate_completeness.csv")
    for source in (root / "tests/analysis/tables").glob("table_*_ref.*"):
        if source.suffix in (".csv", ".md"):
            copy(source, dest / "simulations/tables" / source.name)
    for source in (analysis / "tables").glob("T*.csv*"):
        copy(source, dest / "applied/tables" / source.name)
    # Explicitly select completed revision environments and execution logs.
    for name in ("sessionInfo_08.txt", "sessionInfo_10.txt", "sessionInfo_11.txt", "sessionInfo.txt",
                 "package_versions.csv", "loaded_package.txt", "package_install.log",
                 "08_city_diagnostics_revision.log", "10_profile_comparison.log",
                 "city_run_manifest.json", "10_profile_summary_revision.log",
                 "10_profile_chain1_revision.log", "10_profile_chain2_revision.log"):
        copy(analysis / "logs" / name, dest / "environment" / name)

    execution_source = analysis / "logs/profile_execution_source"
    if execution_source.exists() and not args.source_only:
        for source in execution_source.iterdir():
            if source.is_file():
                copy(source, dest / "environment/profile_execution_source" / source.name)

    # Freeze executable package/analysis source only. Exclude manuscript and
    # all private review materials, independently of repository ignore rules.
    sources = [root / "DESCRIPTION", root / "NAMESPACE"]
    for directory in ("R", "src", "Applied_Study/code"):
        sources.extend(p for p in (root / directory).rglob("*") if p.is_file()
                       and (p.suffix in (".R", ".cpp", ".h", ".hpp", ".c", ".py", ".sh")
                            or p.name in ("Makevars", "Makevars.win")))
    sources.extend(root / f for f in ("tests/_generate_runners.R", "tests/analyze_results.R",
        "tests/sim_count.R", "tests/sim_zi.R", "tests/check_prediction_rows.R",
        "Applied_Study/code/README.md", "README.md"))
    sources.extend((root / "man").glob("*.Rd"))
    sources = sorted(set(p for p in sources if p.exists()))
    with zipfile.ZipFile(dest / "analysis_source_snapshot.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for source in sources:
            archive.write(source, source.relative_to(root))

    def git(*arguments: str) -> str:
        return subprocess.check_output(["git", "-C", str(root), *arguments], text=True).strip()

    manifest = {
        "exported_utc": datetime.now(timezone.utc).isoformat(),
        "base_commit": git("rev-parse", "HEAD"),
        "tags_at_base_commit": git("tag", "--points-at", "HEAD").splitlines(),
        "working_tree_modified": bool(git("status", "--porcelain", "--untracked-files=no")),
        "source_sha256": {str(p.relative_to(root)): sha256(p) for p in sources},
        "data_sha256": {"Applied_Study/lssinc07.csv": sha256(root / "Applied_Study/lssinc07.csv")},
        "capture_phase": "execution_start" if args.source_only else "export",
        "provenance_note": (
            "Source snapshot at the indicated capture phase. Historical results predate this revision and their "
            "execution environments cannot be reconstructed retroactively. Fresh targeted "
            "runs are distinguished in README.md; this export is not a clean tagged full rerun."
        ),
    }
    (dest / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (dest / "artifact_sha256.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["path", "sha256", "bytes"])
        for path in sorted(dest.rglob("*")):
            if path.is_file() and path.name != "artifact_sha256.csv":
                writer.writerow([path.relative_to(dest), sha256(path), path.stat().st_size])
    print(f"Exported final replication artifacts to {dest}")


if __name__ == "__main__":
    main()
