# Glypican-3 (GPC3) across 33 cancers

**Live report:** https://bymankirat.github.io/GPC3-Pan-Cancer-Biomarker/ · **Phase 1 summary:** https://bymankirat.github.io/GPC3-Pan-Cancer-Biomarker/phase1.html

GPC3 is a cell-surface protein that is almost silent in adult tissue but switched back on in liver cancer, which makes it a target for
diagnosis, antibodies and CAR-T cells. In my PG Diploma project (St. Xavier's College, Mumbai, 2025; with Amarthya Vishnu and Aditya) we
tested five cancers with web tools (phase 1). This repository is phase 2: the same questions as one reproducible Python pipeline on
public data, extended to every TCGA cancer, adjusted survival models, regulation, co-expression and structure.

## What I found

| Question | Result |
|---|---|
| Where is GPC3 changed? | 26 cancers could be tested against TCGA + GTEx normal tissue: **up in 5, down in 18** (FDR 5%). Liver cancer stands out (log2 fold change **6.3**); outside the liver GPC3 is more often lower in tumours |
| Can it spot a liver tumour? | **AUC 0.94** (95% CI 0.92–0.96) against normal liver, vs 0.75 for AFP; higher in the tumour in **92%** of the 50 patients with both samples (p < 0.001); same result with TCGA normals only |
| Does it predict survival? | A median split finds nothing (p = 0.25), as in phase 1. Adjusted for age, sex and stage, more GPC3 goes with shorter **overall** (HR 1.22 per SD, p = 0.037) and **disease-specific** survival (HR 1.38, p = 0.013). Late-stage tumours have less GPC3, which hides the effect in crude analyses |
| The "best cutpoint" trap | Searching for the cutpoint with the smallest p gives p = 0.068 (OS) and 0.034 (DSS); corrected by 1,000 permutations: **0.37** and **0.24** |
| What switches it on? | Not copy number (ρ = 0.01). **Methylation**: 11 of 12 probes in GPC3 are less methylated in high-GPC3 tumours (strongest ρ = −0.61) |
| What moves with it? | Foetal-liver genes (AFP, IGF2BP1, IGF2BP2, MSI1, DUSP9) and proliferation (E2F, G2-M, MYC Hallmarks up); normal liver metabolism down. Wnt/β-catenin is not significant (NES 1.25, FDR 0.13) |
| Is it mutated? | Rarely: **1.0%** of 10,443 tumours, scattered, no hotspot; missense changes fall in the folded core at the rate expected by chance |

## Pipeline

| Step | What | Data |
|---|---|---|
| 1 | Download everything, MD5 fingerprints | UCSC Xena (Toil TCGA + GTEx, TCGA Clinical Data Resource, TCGA-LIHC RNA-seq, GISTIC2, 450k methylation), cBioPortal, UniProt, AlphaFold DB |
| 2 | Tumour vs normal in every cancer: Mann–Whitney, BH FDR, AUC with bootstrap CI, TCGA-only check, paired test, AFP | |
| 3 | Survival: KM, Cox per SD (alone and adjusted), proportional-hazards test, cutpoint scan + permutation, all cancers | lifelines |
| 4 | Copy number, methylation, AFP, grade, stage, sex in LIHC | |
| 5 | Genome-wide Spearman co-expression and pre-ranked Hallmark GSEA | gseapy |
| 6 | UniProt features, AlphaFold pLDDT, TCGA mutations on the structure | |
| 7 | Interactive report with a 3D structure viewer | `docs/` |

```bash
pip install -r requirements.txt
python run_all.py          # ~10 minutes, needs internet
python -m pytest -q        # 11 tests
```

The GitHub workflow runs the whole pipeline and the tests, and saves the results here.

## Limits

TCGA adjacent normal tissue is not healthy tissue and GTEx donors differ from patients. RNA is not surface protein. The survival result comes
from one cohort with three endpoints tested and needs an independent cohort. Correlation is not causation.

Mankirat Singh · [Portfolio](https://bymankirat.github.io/myportfolio/) · [GitHub](https://github.com/bymankirat)
