"""Checks on the pipeline output against what is already known. Skipped until the pipeline has run."""
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RES, DER = ROOT / "results", ROOT / "data" / "derived"
needs = pytest.mark.skipif(not (RES / "pancancer.csv").exists(), reason="pipeline not run yet")


@needs
def test_lihc_is_strongly_up_with_high_auc():
    df = pd.read_csv(RES / "pancancer.csv").set_index("cancer")
    lihc = df.loc["LIHC"]
    assert lihc.q < 0.05 and lihc.log2fc > 2 and lihc.auc > 0.8


@needs
def test_sample_counts_are_plausible():
    df = pd.read_csv(RES / "pancancer.csv").set_index("cancer")
    assert 330 <= df.loc["LIHC", "n_tumour"] <= 380        # TCGA-LIHC has ~370 primary tumours
    assert len(df) >= 30                                    # 33 TCGA cancer types


@needs
def test_survival_endpoints_present():
    s = json.loads((RES / "survival.json").read_text())
    assert s["n_patients"] > 300
    for ep in ("OS", "DSS", "PFI"):
        assert s["endpoints"][ep]["events"] > 50


@needs
def test_protein_matches_uniprot_length():
    p = json.loads((DER / "protein.json").read_text())
    assert p["length"] == 580 and p["alphafold"]["same_sequence"]
    st = json.loads((RES / "structure.json").read_text())
    assert len(st["plddt"]) == 580


@needs
def test_results_js_is_valid():
    js = (ROOT / "docs" / "results.js").read_text()
    data = json.loads(js[len("window.GPC3 = "):].rstrip().rstrip(";"))
    assert data["pancancer"] and data["survival"]
