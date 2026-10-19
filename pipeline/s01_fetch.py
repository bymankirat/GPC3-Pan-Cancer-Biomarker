"""Step 1 (needs internet): get every public dataset the analysis uses.

- GPC3 and AFP expression for every TCGA tumour, TCGA normal and GTEx normal sample
  (UCSC Xena Toil hub, processed together so tumour and normal are comparable)
- sample annotations (cancer type, tissue, tumour or normal)
- TCGA Clinical Data Resource: overall, disease-specific and progression-free survival, age, sex, stage
- TCGA-LIHC whole transcriptome, GPC3 copy number and promoter methylation
- GPC3 somatic mutations in the TCGA PanCancer Atlas (cBioPortal)
- GPC3 protein features (UniProt) and AlphaFold model

Small, derived tables are saved to data/derived/ with MD5 fingerprints in provenance.json.
"""
import gzip
import io
import time

import numpy as np
import pandas as pd

from util import DERIVED, RAW, config, download, http, log, md5, save_json

prov = {"downloaded": time.strftime("%Y-%m-%d"), "files": {}}


def xena():
    import xenaPython as x
    return x


def note(name, path, source):
    prov["files"][name] = {"source": source, "md5": md5(path), "bytes": path.stat().st_size}


def read_table(url, raw_name, hub=None, dataset=None):
    """Download a Xena file; if the direct link fails, ask the hub for it."""
    dest = RAW / raw_name
    try:
        download(url, dest)
        data = dest.read_bytes()
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        return pd.read_csv(io.BytesIO(data), sep="\t", low_memory=False, encoding="latin-1")
    except Exception as e:
        log(f"  direct download failed ({str(e)[:80]}); using the Xena API")
    x = xena()
    samples = x.dataset_samples(hub, dataset, None)
    fields = [f for f in x.dataset_field(hub, dataset) if f not in ("sampleID", "sample")]
    out = {"sample": samples}
    for f in fields:
        codes = x.field_codes(hub, dataset, [f])[0].get("code")
        vals = x.dataset_probe_values(hub, dataset, samples, [f])[1][0]
        if codes:
            codes = codes.split("\t")
            vals = [codes[int(float(v))] if v is not None and str(v) != "NaN" and float(v) == float(v) else None for v in vals]
        out[f] = vals
    return pd.DataFrame(out)


def read_matrix(url, raw_name, hub, dataset):
    """Genomic matrix (genes x samples). Direct download first, chunked API queries as fallback."""
    dest = RAW / raw_name
    try:
        download(url, dest)
        df = pd.read_csv(dest, sep="\t", index_col=0, compression="gzip" if dest.read_bytes()[:2] == b"\x1f\x8b" else None)
        return df
    except Exception as e:
        log(f"  direct download failed ({str(e)[:80]}); using the Xena API")
    x = xena()
    samples = x.dataset_samples(hub, dataset, None)
    genes = x.dataset_field(hub, dataset)
    rows = []
    for i in range(0, len(genes), 1000):
        chunk = genes[i:i + 1000]
        vals = x.dataset_probe_values(hub, dataset, samples, chunk)[1]
        rows += [[np.nan if v is None else v for v in r] for r in vals]
    return pd.DataFrame(rows, index=genes, columns=samples)


def toil_expression(samples):
    x = xena()
    fields = x.dataset_field(config.TOIL_HUB, config.TOIL_EXPR)
    probes = {g: next(f for f in fields if f.split(".")[0] == ens) for g, ens in
              ((config.GENE, config.GENE_ENSEMBL), (config.PARTNER_MARKER, config.PARTNER_ENSEMBL))}
    out = {"sample": samples}
    for gene, probe in probes.items():
        vals = []
        for i in range(0, len(samples), 3000):
            vals += x.dataset_probe_values(config.TOIL_HUB, config.TOIL_EXPR, samples[i:i + 3000], [probe])[1][0]
        v = np.array([np.nan if (a is None or a != a) else a for a in vals], dtype=float)
        out[gene] = np.log2(np.clip(2 ** v - 0.001, 0, None) + 1)   # log2(TPM + 0.001) -> log2(TPM + 1)
    log(f"  toil probes: {probes}")
    return pd.DataFrame(out), probes


def main():
    log("sample annotations (Toil)")
    pheno = read_table(config.TOIL_PHENO_URL, "toil_phenotype.txt.gz", config.TOIL_HUB, config.TOIL_PHENO)
    pheno.columns = [c.strip() for c in pheno.columns]
    pheno = pheno.rename(columns={"primary disease or tissue": "disease", "_primary_site": "site", "_sample_type": "sample_type",
                                  "_study": "study", "_gender": "sex", "detailed_category": "detail"})
    log(f"  {len(pheno)} samples")

    log("GPC3 and AFP expression across TCGA and GTEx (Toil)")
    x = xena()
    expr_samples = x.dataset_samples(config.TOIL_HUB, config.TOIL_EXPR, None)
    expr, probes = toil_expression(expr_samples)
    pan = pheno.merge(expr, on="sample", how="inner")
    pan = pan[pan.study.isin(["TCGA", "GTEX"])]
    pan.to_csv(DERIVED / "pancancer_expression.csv.gz", index=False)
    note("pancancer_expression", DERIVED / "pancancer_expression.csv.gz", f"{config.TOIL_HUB} {config.TOIL_EXPR} probes {probes}")
    log(f"  {len(pan)} TCGA + GTEx samples")

    log("survival (TCGA Clinical Data Resource)")
    surv = read_table(config.SURVIVAL_URL, "tcga_cdr.tsv", config.PANCAN_HUB, config.SURVIVAL)
    keep = [c for c in surv.columns if c in ("sample", "_PATIENT", "cancer type abbreviation", "age_at_initial_pathologic_diagnosis",
                                             "gender", "ajcc_pathologic_tumor_stage", "histological_grade", "OS", "OS.time",
                                             "DSS", "DSS.time", "PFI", "PFI.time", "DFI", "DFI.time")]
    surv[keep].to_csv(DERIVED / "survival.csv.gz", index=False)
    note("survival", DERIVED / "survival.csv.gz", f"{config.PANCAN_HUB} {config.SURVIVAL}")

    log("TCGA-LIHC transcriptome")
    lihc = read_matrix(config.LIHC_EXPR_URL, "lihc_hiseqv2.gz", config.TCGA_HUB, config.LIHC_EXPR)
    lihc.index.name = "gene"
    lihc = lihc[(lihc > 1).mean(axis=1) >= 0.2]   # keep genes expressed in at least 20% of samples
    lihc.round(3).to_csv(RAW / "lihc_expression.csv.gz")   # too large to commit; rebuilt by this step
    note("lihc_expression", RAW / "lihc_expression.csv.gz", config.LIHC_EXPR_URL)
    log(f"  {lihc.shape[0]} expressed genes x {lihc.shape[1]} samples")

    log("GPC3 copy number and methylation in LIHC")
    samples = list(lihc.columns)
    extra = {"sample": samples}
    try:
        cnv_samples = x.dataset_samples(config.TCGA_HUB, config.LIHC_CNV, None)
        common = [s for s in samples if s in set(cnv_samples)]
        vals = x.dataset_probe_values(config.TCGA_HUB, config.LIHC_CNV, common, [config.GENE])[1][0]
        cn = dict(zip(common, vals))
        extra["gpc3_cnv"] = [cn.get(s) for s in samples]
        prov["files"]["cnv"] = {"source": f"{config.TCGA_HUB} {config.LIHC_CNV}"}
    except Exception as e:
        log("  copy number not available:", str(e)[:120])
    meth_probes = []
    try:
        meth_samples = x.dataset_samples(config.TCGA_HUB, config.LIHC_METH, None)
        common = [s for s in samples if s in set(meth_samples)]
        res = x.dataset_gene_probes_values(config.TCGA_HUB, config.LIHC_METH, common, [config.GENE])
        # result: [ {name, position...} list , values ] in the hub's format; handle the common shapes
        info, values = res[0], res[1]
        names = [p["name"] if isinstance(p, dict) else p for p in (info if isinstance(info, list) else info.get("name", []))]
        for name, v in zip(names, values):
            d = dict(zip(common, v))
            extra[f"meth_{name}"] = [d.get(s) for s in samples]
            meth_probes.append(name)
        prov["files"]["methylation"] = {"source": f"{config.TCGA_HUB} {config.LIHC_METH}", "probes": meth_probes}
        log(f"  {len(meth_probes)} methylation probes in {config.GENE}")
    except Exception as e:
        log("  methylation not available:", str(e)[:120])
    pd.DataFrame(extra).to_csv(DERIVED / "lihc_cnv_methylation.csv", index=False)

    log("GPC3 mutations (cBioPortal, TCGA PanCancer Atlas)")
    try:
        studies = http("GET", f"{config.CBIOPORTAL}/studies", params={"keyword": "tcga_pan_can_atlas_2018", "projection": "SUMMARY"}).json()
        ids = [s["studyId"] for s in studies if s["studyId"].endswith("_tcga_pan_can_atlas_2018")]
        muts = http("POST", f"{config.CBIOPORTAL}/mutations/fetch", params={"projection": "DETAILED"},
                    json={"molecularProfileIds": [f"{s}_mutations" for s in ids], "entrezGeneIds": [config.GENE_ENTREZ]}).json()
        rows = [{"study": m.get("studyId"), "sample": m.get("sampleId"), "protein_change": m.get("proteinChange"),
                 "type": m.get("mutationType"), "start": m.get("proteinPosStart"), "end": m.get("proteinPosEnd"),
                 "gene": (m.get("gene") or {}).get("hugoGeneSymbol")} for m in muts]
        mdf = pd.DataFrame(rows)
        if len(mdf) and set(mdf.gene.dropna()) != {config.GENE}:
            raise RuntimeError(f"unexpected genes in cBioPortal response: {set(mdf.gene)}")
        mdf.to_csv(DERIVED / "mutations.csv", index=False)
        # sequenced samples per study, for mutation frequency
        counts = {}
        for s in ids:
            lst = http("GET", f"{config.CBIOPORTAL}/sample-lists/{s}_sequenced", ok404=True)
            counts[s] = len(lst.json().get("sampleIds", [])) if lst is not None else None
        save_json(DERIVED / "mutation_denominators.json", counts)
        prov["files"]["mutations"] = {"source": f"{config.CBIOPORTAL} {len(ids)} PanCancer Atlas studies", "n": len(mdf)}
        log(f"  {len(mdf)} mutations in {len(ids)} studies")
    except Exception as e:
        log("  mutations not available:", str(e)[:160])

    log("UniProt features and AlphaFold model")
    up = http("GET", config.UNIPROT_JSON.format(acc=config.UNIPROT)).json()
    feats = [{"type": f["type"], "start": f["location"]["start"].get("value"), "end": f["location"]["end"].get("value"),
              "description": f.get("description", "")} for f in up.get("features", [])]
    seq = up["sequence"]["value"]
    af = http("GET", config.ALPHAFOLD_API.format(acc=config.UNIPROT)).json()[0]
    pdb = http("GET", af["pdbUrl"]).text
    (DERIVED / "alphafold.pdb").write_text(pdb)
    save_json(DERIVED / "protein.json", {"accession": config.UNIPROT, "sequence": seq, "length": len(seq), "features": feats,
                                         "alphafold": {"model": af.get("entryId"), "version": af.get("latestVersion"),
                                                       "same_sequence": af.get("uniprotSequence", seq) == seq}})
    prov["files"]["uniprot"] = {"source": config.UNIPROT_JSON.format(acc=config.UNIPROT), "length": len(seq)}
    prov["files"]["alphafold"] = {"source": af["pdbUrl"], "md5": md5(DERIVED / "alphafold.pdb")}

    save_json(DERIVED / "provenance.json", prov)
    log("done")


if __name__ == "__main__":
    main()
