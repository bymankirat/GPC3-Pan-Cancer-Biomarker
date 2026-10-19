"""Step 2: GPC3 in every TCGA cancer vs normal tissue.

Tumour = TCGA primary tumours of one cancer type. Normal = TCGA adjacent normal
samples of the same cancer plus GTEx samples of the same organ (the same idea
as GEPIA2, but every test is done here and corrected for testing 30+ cancers).

For each cancer: median difference (log2 fold change), Mann-Whitney test,
Benjamini-Hochberg q, and AUC with a bootstrap 95% CI: how well GPC3 alone tells
a tumour from normal tissue. A sensitivity check repeats the test with TCGA
normals only (no GTEx), and LIHC also gets a paired test on patients who gave
both a tumour and a normal sample.
"""
import numpy as np
import pandas as pd
from scipy import stats as st

from stats import auc, bh, bootstrap_auc
from util import DERIVED, RESULTS, config, log, save_json

TUMOUR_TYPES = {"Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"}
NORMAL_TYPES = {"Solid Tissue Normal", "Normal Tissue"}


def load():
    pan = pd.read_csv(DERIVED / "pancancer_expression.csv.gz")
    surv = pd.read_csv(DERIVED / "survival.csv.gz")
    abbr = dict(zip(surv["_PATIENT"], surv["cancer type abbreviation"]))
    pan["patient"] = pan["sample"].str[:12]
    pan["cancer"] = np.where(pan.study == "TCGA", pan.patient.map(abbr), None)
    return pan


def gtex_match(pan, site):
    """GTEx samples of the organ a TCGA cancer comes from (names differ for a few organs)."""
    gt = pan.study == "GTEX"
    if site in config.GTEX_SITE:
        m = config.GTEX_SITE[site]
        if m is None:
            return gt & False
        tissue, detail = m
        return gt & (pan.site == tissue) & ((pan.detail == detail) if detail else True)
    return gt & (pan.site.str.lower() == str(site).lower())


def quantiles(v):
    return [round(float(x), 3) for x in np.percentile(v, [5, 25, 50, 75, 95])]


def main():
    rng = np.random.default_rng(config.SEED)
    pan = load()
    g = config.GENE
    tum = pan[(pan.study == "TCGA") & pan.sample_type.isin(TUMOUR_TYPES) & pan.cancer.notna() & pan[g].notna()]
    rows, dist = [], {}
    for c, t in tum.groupby("cancer"):
        site = t.site.mode().iat[0]
        normal = pan[pan.sample_type.isin(NORMAL_TYPES) & pan[g].notna() & (((pan.study == "TCGA") & (pan.cancer == c)) | gtex_match(pan, site))]
        tcga_n = normal[normal.study == "TCGA"]
        row = {"cancer": c, "disease": t.disease.mode().iat[0], "site": site, "n_tumour": len(t), "n_normal": len(normal),
               "n_normal_tcga": len(tcga_n), "n_normal_gtex": int((normal.study == "GTEX").sum()),
               "median_tumour": float(t[g].median())}
        if len(normal) >= config.MIN_GROUP and len(t) >= config.MIN_GROUP:
            a = auc(t[g], normal[g])
            lo, hi = bootstrap_auc(t[g].values, normal[g].values, config.N_BOOT, rng)
            row.update({"median_normal": float(normal[g].median()), "log2fc": float(t[g].median() - normal[g].median()),
                        "p": float(st.mannwhitneyu(t[g], normal[g], alternative="two-sided").pvalue),
                        "auc": a, "auc_lo": lo, "auc_hi": hi})
            dist[c] = {"tumour": quantiles(t[g]), "normal": quantiles(normal[g])}
        else:
            dist[c] = {"tumour": quantiles(t[g])}
        if len(tcga_n) >= config.MIN_GROUP:
            row.update({"log2fc_tcga_only": float(t[g].median() - tcga_n[g].median()),
                        "p_tcga_only": float(st.mannwhitneyu(t[g], tcga_n[g], alternative="two-sided").pvalue)})
        rows.append(row)
    df = pd.DataFrame(rows)
    df["q"] = bh(df["p"].values) if "p" in df else np.nan
    if "p_tcga_only" in df:
        df["q_tcga_only"] = bh(df["p_tcga_only"].values)
    df["direction"] = np.where(df.q < config.FDR, np.where(df.log2fc > 0, "up", "down"), "ns")
    df = df.sort_values("log2fc", ascending=False, na_position="last")
    df.to_csv(RESULTS / "pancancer.csv", index=False)

    # paired tumour vs normal in the focus cancer
    f = config.FOCUS_CANCER
    ft = tum[tum.cancer == f].drop_duplicates("patient").set_index("patient")[g]
    fn = pan[(pan.study == "TCGA") & (pan.cancer == f) & (pan.sample_type == "Solid Tissue Normal")].drop_duplicates("patient").set_index("patient")[g]
    both = ft.index.intersection(fn.index)
    w = st.wilcoxon(ft[both], fn[both])
    paired = {"n_pairs": len(both), "median_difference": float((ft[both] - fn[both]).median()),
              "fraction_higher_in_tumour": float((ft[both] > fn[both]).mean()), "p": float(w.pvalue),
              "pairs": [[round(float(fn[p]), 3), round(float(ft[p]), 3)] for p in both]}
    # AFP, the marker GPC3 is usually compared with
    afp = None
    if config.PARTNER_MARKER in pan:
        m = config.PARTNER_MARKER
        t = tum[tum.cancer == f]
        site = t.site.mode().iat[0]
        normal = pan[pan.sample_type.isin(NORMAL_TYPES) & (((pan.study == "TCGA") & (pan.cancer == f)) | gtex_match(pan, site))]
        both_t = t.dropna(subset=[g, m])
        both_n = normal.dropna(subset=[g, m])
        afp = {"auc_gpc3": auc(both_t[g], both_n[g]), "auc_afp": auc(both_t[m], both_n[m]),
               "auc_gpc3_ci": bootstrap_auc(both_t[g].values, both_n[g].values, config.N_BOOT, rng),
               "auc_afp_ci": bootstrap_auc(both_t[m].values, both_n[m].values, config.N_BOOT, rng),
               "spearman_in_tumours": float(st.spearmanr(both_t[g], both_t[m]).statistic)}
    save_json(RESULTS / "pancancer.json", {"distributions": dist, "paired": paired, "afp": afp,
                                           "n_cancers_tested": int(df.p.notna().sum()),
                                           "n_up": int((df.direction == "up").sum()), "n_down": int((df.direction == "down").sum())})
    top = df.dropna(subset=["auc"]).iloc[0]
    log(f"{int(df.p.notna().sum())} cancers tested; up in {int((df.direction == 'up').sum())}, down in {int((df.direction == 'down').sum())}; "
        f"largest: {top.cancer} log2FC {top.log2fc:.2f}, AUC {top.auc:.3f}")
    log(f"{f} paired: {paired['n_pairs']} pairs, tumour higher in {paired['fraction_higher_in_tumour']:.0%}")


if __name__ == "__main__":
    main()
