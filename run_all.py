"""Run the whole GPC3 analysis.

    python run_all.py            # everything (needs internet)
    python run_all.py --offline  # steps 2-7 from the data already downloaded

The GitHub workflow in .github/workflows/pipeline.yml runs the full version.
"""
import argparse
import importlib
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "pipeline"))

STEPS = [
    ("s01_fetch", "download TCGA, GTEx, survival, copy number, methylation, mutations, UniProt, AlphaFold", True),
    ("s02_pancancer", "GPC3 in every cancer vs normal tissue", False),
    ("s03_survival", "survival: Kaplan-Meier, Cox, cutpoint permutation, pan-cancer", False),
    ("s04_drivers", "copy number, methylation, AFP, sex, grade and stage in LIHC", False),
    ("s05_coexpression", "co-expressed genes and Hallmark GSEA in LIHC", False),
    ("s06_structure", "protein features, AlphaFold confidence and mutations", False),
    ("s07_export", "write docs/results.js", False),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true", help="skip the steps that need internet")
    ap.add_argument("--from", dest="start", type=int, default=1, help="first step to run (1-7)")
    args = ap.parse_args()
    t0 = time.time()
    for i, (mod, what, online) in enumerate(STEPS, 1):
        if i < args.start or (online and args.offline):
            continue
        print(f"\n=== step {i}: {what} ===", flush=True)
        importlib.import_module(mod).main()
    print(f"\nfinished in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
