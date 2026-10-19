"""Step 5: which genes rise and fall with GPC3 in liver tumours, and which pathways do they belong to?

Spearman correlation of GPC3 with every expressed gene in TCGA-LIHC primary tumours,
then a pre-ranked gene set enrichment analysis (GSEA) on the 50 MSigDB Hallmark sets.
Phase 1's STRING network pointed to Wnt signalling; this asks whether tumours with more
GPC3 actually show more Wnt/beta-catenin activity.
"""
import numpy as np
import pandas as pd
from scipy import stats as st

from stats import bh
from util import RAW, RESULTS, config, log, save_json


def main():
    src = RAW / "lihc_expression.csv.gz"
    out_json = RESULTS / "coexpression.json"
    if not src.exists():
        if out_json.exists():
            log("LIHC matrix not downloaded; keeping existing co-expression results")
            return
        raise SystemExit("run step 1 first")
    ex = pd.read_csv(src, index_col=0)
    tum = ex[[c for c in ex.columns if c.endswith("-01")]]
    g = config.GENE
    x = tum.loc[g]
    others = tum.drop(index=g)
    xr = st.rankdata(x)
    ranks = np.apply_along_axis(st.rankdata, 1, others.values)
    xr = (xr - xr.mean()) / xr.std()
    rr = (ranks - ranks.mean(axis=1, keepdims=True)) / ranks.std(axis=1, keepdims=True)
    rho = (rr * xr).mean(axis=1)
    n = len(x)
    tstat = rho * np.sqrt((n - 2) / np.clip(1 - rho ** 2, 1e-12, None))
    p = 2 * st.t.sf(np.abs(tstat), n - 2)
    co = pd.DataFrame({"gene": others.index, "rho": rho, "p": p}).dropna()
    co["q"] = bh(co.p.values)
    co = co.sort_values("rho", ascending=False)
    co.round(5).to_csv(RESULTS / "coexpression.csv.gz", index=False)
    log(f"{len(co)} genes; {int(((co.q < config.FDR) & (co.rho > 0)).sum())} positive and "
        f"{int(((co.q < config.FDR) & (co.rho < 0)).sum())} negative at q<0.05")

    gsea = []
    try:
        import gseapy as gp
        rnk = co.drop_duplicates("gene").set_index("gene").rho.sort_values(ascending=False)
        res = gp.prerank(rnk=rnk, gene_sets=config.HALLMARK_LIBRARY, permutation_num=config.GSEA_PERM, seed=config.SEED,
                         min_size=15, max_size=500, threads=4, outdir=None, verbose=False)
        r = res.res2d.copy()
        r.columns = [c.lower().replace(" ", "_").replace("-", "_") for c in r.columns]
        for _, row in r.iterrows():
            gsea.append({"term": row.get("term"), "nes": float(row.get("nes")), "es": float(row.get("es")),
                         "p": float(row.get("nom_p_val")), "fdr": float(row.get("fdr_q_val")),
                         "lead": str(row.get("lead_genes", ""))[:400]})
        gsea.sort(key=lambda d: -d["nes"])
        log(f"GSEA: {sum(d['fdr'] < 0.05 for d in gsea)} Hallmark sets at FDR<0.05")
    except Exception as e:
        log("GSEA not run:", str(e)[:160])

    save_json(out_json, {"n_tumours": n, "n_genes": len(co),
                         "top_positive": co.head(25).round(4).to_dict(orient="records"),
                         "top_negative": co.tail(25).iloc[::-1].round(4).to_dict(orient="records"),
                         "histogram": np.histogram(co.rho, bins=np.linspace(-1, 1, 41))[0].tolist(),
                         "gsea": gsea})


if __name__ == "__main__":
    main()
