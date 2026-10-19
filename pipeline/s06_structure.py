"""Step 6: the GPC3 protein, its AlphaFold model and where cancer mutations fall on it.

- UniProt features: signal peptide, chain, cleavage site, heparan sulfate attachment sites, GPI anchor
- AlphaFold per-residue confidence (pLDDT); below 50 usually means a flexible, disordered region
- Somatic mutations in ~10,000 TCGA tumours (PanCancer Atlas, cBioPortal): how common, in which
  cancers, and whether missense changes fall in the folded core more or less than expected by chance
"""
import numpy as np
import pandas as pd
from scipy import stats as st

from util import DERIVED, RESULTS, config, load_json, log, save_json


def plddt(pdb_text):
    vals = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            vals[int(line[22:26])] = float(line[60:66])
    return [vals[i] for i in sorted(vals)]


def main():
    prot = load_json(DERIVED / "protein.json")
    pl = plddt((DERIVED / "alphafold.pdb").read_text())
    L = prot["length"]
    feats = [f for f in prot["features"] if f["type"] in ("Signal", "Chain", "Propeptide", "Lipidation", "Glycosylation",
                                                          "Site", "Disulfide bond", "Domain", "Region", "Mutagenesis")]
    out = {"length": L, "plddt": pl, "features": feats, "alphafold": prot["alphafold"],
           "ordered_fraction": float(np.mean(np.array(pl) >= 70)), "disordered_fraction": float(np.mean(np.array(pl) < 50))}

    mf = DERIVED / "mutations.csv"
    if mf.exists() and mf.stat().st_size > 10:
        m = pd.read_csv(mf)
        den = load_json(DERIVED / "mutation_denominators.json") if (DERIVED / "mutation_denominators.json").exists() else {}
        m["cancer"] = m.study.str.replace("_tcga_pan_can_atlas_2018", "").str.upper()
        per = m.groupby("cancer")["sample"].nunique().rename("mutated").reset_index()
        per["sequenced"] = per.cancer.map({k.replace("_tcga_pan_can_atlas_2018", "").upper(): v for k, v in den.items()})
        per["percent"] = 100 * per.mutated / per.sequenced
        total_seq = sum(v for v in den.values() if v)
        mis = m[m.type == "Missense_Mutation"].dropna(subset=["start"])
        pos = mis.start.astype(int).clip(1, L)
        ordered = np.array(pl) >= 70
        obs = int(ordered[pos - 1].sum())
        expected_frac = float(ordered.mean())
        binom = st.binomtest(obs, len(pos), expected_frac) if len(pos) else None
        hot = pos.value_counts()
        out["mutations"] = {
            "total": int(len(m)), "samples": int(m["sample"].nunique()), "sequenced": int(total_seq),
            "percent_samples": 100 * m["sample"].nunique() / total_seq if total_seq else None,
            "by_type": m.type.value_counts().to_dict(), "by_cancer": per.sort_values("mutated", ascending=False).round(2).to_dict(orient="records"),
            "missense_in_ordered": obs, "missense_total": int(len(pos)), "ordered_fraction_expected": expected_frac,
            "ordered_binom_p": float(binom.pvalue) if binom else None,
            "recurrent": [{"position": int(k), "n": int(v)} for k, v in hot[hot >= 3].items()],
            "list": m[["cancer", "protein_change", "type", "start"]].fillna("").to_dict(orient="records")}
        log(f"{len(m)} mutations in {m["sample"].nunique()} of {total_seq} tumours; missense in ordered core {obs}/{len(pos)} "
            f"(expected {expected_frac:.0%})")
    save_json(RESULTS / "structure.json", out)
    log(f"AlphaFold: {out['ordered_fraction']:.0%} residues pLDDT>=70, {out['disordered_fraction']:.0%} below 50")


if __name__ == "__main__":
    main()
