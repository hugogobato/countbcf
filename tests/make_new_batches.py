"""Group the NEW benchmark-method runners into Colab batch folders of three.

The original two methods (countbcf, bcf_gauss) already live in
Experiment_Code_in_Notebooks/batch_01 .. batch_38.  This script converts only
the runners for the added parametric / causal-forest baselines

    glm_poisson, glm_nb   (count study)
    glm_zip, glm_zinb     (zi study)
    causal_forest         (both studies)

into notebooks and writes them, three per folder, into fresh batch_NN
directories that continue the existing numbering.  It is safe to re-run: any
batch folder that contains *only* new-method notebooks is rebuilt from scratch,
while the original batches are left untouched.

    python make_new_batches.py
"""

import json
import re
from pathlib import Path

from convert_to_notebooks import r_to_notebook

NEW_METHODS = ["glm_poisson", "glm_nb", "glm_zip", "glm_zinb", "causal_forest"]
NEW_RE = re.compile(r"^run_(?:count|zi)_(?:" + "|".join(map(re.escape, NEW_METHODS)) + r")_")
BATCH_RE = re.compile(r"^batch_(\d+)$")
CHUNK = 3

BASE = Path(__file__).parent
RUNNER_DIRS = [BASE / "runners_count", BASE / "runners_zi"]
OUT_ROOT = BASE / "Experiment_Code_in_Notebooks"


def is_new(name: str) -> bool:
    return bool(NEW_RE.match(name))


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    # 1. Tear down any previously generated all-new-method batch folders so the
    #    script is idempotent without ever touching the original batches.
    for folder in sorted(OUT_ROOT.glob("batch_*")):
        if not folder.is_dir():
            continue
        nbs = list(folder.glob("*.ipynb"))
        if nbs and all(is_new(p.name) for p in nbs):
            for p in nbs:
                p.unlink()
            folder.rmdir()

    # 2. Next free batch number after the surviving (original) batches.
    nums = []
    for folder in OUT_ROOT.glob("batch_*"):
        m = BATCH_RE.match(folder.name)
        if m:
            nums.append(int(m.group(1)))
    batch = max(nums) + 1 if nums else 1

    # 3. Collect the new-method runners: count study first, then zi; sorted.
    runners = []
    for d in RUNNER_DIRS:
        runners.extend(sorted(p for p in d.glob("*.R") if is_new(p.name)))

    # 4. Emit notebooks, CHUNK per folder.
    first = batch
    for i in range(0, len(runners), CHUNK):
        folder = OUT_ROOT / f"batch_{batch:02d}"
        folder.mkdir(parents=True, exist_ok=True)
        for r_file in runners[i:i + CHUNK]:
            nb = r_to_notebook(r_file)
            out = folder / r_file.with_suffix(".ipynb").name
            out.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
        batch += 1

    n_batches = batch - first
    print(f"Wrote {len(runners)} notebooks into "
          f"batch_{first:02d}..batch_{batch - 1:02d} ({n_batches} folders).")


if __name__ == "__main__":
    main()
