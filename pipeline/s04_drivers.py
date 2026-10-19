"""Step 4: what goes with high GPC3 in liver tumours?

- Copy number (GISTIC2 thresholded calls): is GPC3 high because its region is gained?
- Promoter methylation (450k probes in the gene): is it switched on by losing methylation?
- AFP, the classic liver-cancer marker: are the two markers high in the same tumours?
- Sex (GPC3 is on the X chromosome), tumour grade and stage.
"""
import numpy as np
import pandas as pd
from scipy import stats as st

from s03_survival import stage_group
from stats import bh
from util import DERIVED, RESULTS, config, log, save_json


def main():
    g, m = config.GENE, config.PARTNER_MARKER
    pan = pd.read_csv(DERIVED / "pancancer_expression.csv.gz")
    surv = pd.read_csv(DERIVED / "survival.csv.gz").drop_duplicates("_PATIENT")
    extra = pd.read_csv(DERIVED / "lihc_cnv_methylation.csv")
    t = pan[(pan.study == "TCGA") & (pan.sample_type == "Primary Tumor")].copy()
    t["patient"] = t["sample"].str[:12]
    t = t.merge(surv[["_PATIENT", "cancer type abbreviation", "gender", "histological_grade", "ajcc_pathologic_tumor_stage"]],
                left_on="patient", right_on="_PATIENT")
    t = t[t["cancer type abbreviation"] == config.FOCUS_CANCER].drop_duplicates("patient")
    t = t.merge(extra, on="sample", how="left")
    out = {"n": len(t)}

    if "gpc3_cnv" in t and t.gpc3_cnv.notna().sum() > 30:
        d = t.dropna(subset=["gpc3_cnv"])
        groups = {int(k): v[g].tolist() for k, v in d.groupby(d.gpc3_cnv.round().astype(int))}
        out["cnv"] = {"n": len(d), "spearman": float(st.spearmanr(d.gpc3_cnv, d[g]).statistic),
                      "p": float(st.spearmanr(d.gpc3_cnv, d[g]).pvalue),
                      "groups": {k: {"n": len(v), "median": float(np.median(v)), "q": np.percentile(v, [5, 25, 50, 75, 95]).round(3).tolist()}
                                 for k, v in sorted(groups.items())}}
        log(f"copy number: rho={out['cnv']['spearman']:.2f}, calls {{k: n}} = {{{', '.join(f'{k}: {v['n']}' for k, v in out['cnv']['groups'].items())}}}")

    meth = [c for c in t.columns if c.startswith("meth_")]
    if meth:
        rows = []
        for c in meth:
            d = t[[c, g]].dropna()
            if len(d) > 30:
                r = st.spearmanr(d[c], d[g])
                rows.append({"probe": c[5:], "n": len(d), "rho": float(r.statistic), "p": float(r.pvalue),
                             "median_beta": float(d[c].median())})
        md = pd.DataFrame(rows)
        if len(md):
            md["q"] = bh(md.p.values)
            md.sort_values("rho").to_csv(RESULTS / "methylation_probes.csv", index=False)
            best = md.sort_values("rho").iloc[0]
            d = t[["meth_" + best.probe, g]].dropna()
            out["methylation"] = {"probes": len(md), "negative_sig": int(((md.q < config.FDR) & (md.rho < 0)).sum()),
                                  "positive_sig": int(((md.q < config.FDR) & (md.rho > 0)).sum()),
                                  "best": best.to_dict(), "scatter": d.round(3).values.tolist(),
                                  "table": md.round(4).to_dict(orient="records")}
            log(f"methylation: {len(md)} probes, strongest rho={best.rho:.2f} ({best.probe})")

    d = t[[g, m]].dropna()
    r = st.spearmanr(d[g], d[m])
    hi_g, hi_m = d[g] > d[g].median(), d[m] > d[m].median()
    out["afp"] = {"n": len(d), "rho": float(r.statistic), "p": float(r.pvalue), "scatter": d.round(3).values.tolist(),
                  "both_high": float((hi_g & hi_m).mean()), "gpc3_only": float((hi_g & ~hi_m).mean()),
                  "afp_only": float((~hi_g & hi_m).mean())}
    log(f"AFP: rho={r.statistic:.2f}")

    sex = t.dropna(subset=["gender"])
    a, b = sex[sex.gender.str.upper() == "MALE"][g], sex[sex.gender.str.upper() == "FEMALE"][g]
    out["sex"] = {"n_male": len(a), "n_female": len(b), "median_male": float(a.median()), "median_female": float(b.median()),
                  "p": float(st.mannwhitneyu(a, b).pvalue)}

    gr = t.assign(grade=t.histological_grade.astype(str).str.extract(r"G(\d)")[0].astype(float)).dropna(subset=["grade"])
    out["grade"] = {"rho": float(st.spearmanr(gr.grade, gr[g]).statistic), "p": float(st.spearmanr(gr.grade, gr[g]).pvalue),
                    "groups": {f"G{int(k)}": {"n": len(v), "median": float(v[g].median())} for k, v in gr.groupby("grade")}}
    stg = t.assign(late=t.ajcc_pathologic_tumor_stage.map(stage_group)).dropna(subset=["late"])
    out["stage"] = {"n_early": int((stg.late == 0).sum()), "n_late": int((stg.late == 1).sum()),
                    "median_early": float(stg[stg.late == 0][g].median()), "median_late": float(stg[stg.late == 1][g].median()),
                    "p": float(st.mannwhitneyu(stg[stg.late == 0][g], stg[stg.late == 1][g]).pvalue)}
    log(f"grade rho={out['grade']['rho']:.2f} p={out['grade']['p']:.3g}; sex p={out['sex']['p']:.3g}; stage p={out['stage']['p']:.3g}")
    save_json(RESULTS / "drivers.json", out)


if __name__ == "__main__":
    main()
