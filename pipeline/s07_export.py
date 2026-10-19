"""Step 7: collect results into docs/results.js (and the AlphaFold model for the 3D viewer)."""
import json
import shutil

import numpy as np
import pandas as pd

from util import DERIVED, RESULTS, ROOT, config, load_json, log


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, np.integer):
        return int(x)
    return x


def table(path):
    return json.loads(pd.read_csv(path).replace({np.nan: None}).to_json(orient="records")) if path.exists() else []


def main():
    j = lambda name: load_json(RESULTS / name) if (RESULTS / name).exists() else None
    prov = load_json(DERIVED / "provenance.json")
    out = {"meta": {"gene": config.GENE, "focus": config.FOCUS_CANCER, "downloaded": prov.get("downloaded"),
                    "sources": {k: v.get("source") for k, v in prov["files"].items()}},
           "pancancer": table(RESULTS / "pancancer.csv"), "pancancer_extra": j("pancancer.json"),
           "survival": j("survival.json"), "survival_pancancer": table(RESULTS / "survival_pancancer.csv"),
           "drivers": j("drivers.json"), "coexpression": j("coexpression.json"), "structure": j("structure.json")}
    js = "window.GPC3 = " + json.dumps(clean(out), separators=(",", ":"), allow_nan=False) + ";\n"
    (ROOT / "docs" / "results.js").write_text(js)
    shutil.copy(DERIVED / "alphafold.pdb", ROOT / "docs" / "alphafold.pdb")
    log(f"docs/results.js written ({len(js) / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
