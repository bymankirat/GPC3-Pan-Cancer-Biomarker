"""Settings for the GPC3 pan-cancer pipeline. Everything that decides what the analysis does lives here."""

GENE = "GPC3"
GENE_ENSEMBL = "ENSG00000147257"
GENE_ENTREZ = 2719
UNIPROT = "P51654"
PARTNER_MARKER = "AFP"          # alpha-fetoprotein, the classic serum marker of liver cancer
PARTNER_ENSEMBL = "ENSG00000081051"
FOCUS_CANCER = "LIHC"

# UCSC Xena hubs
TOIL_HUB = "https://toil.xenahubs.net"
TOIL_EXPR = "TcgaTargetGtex_rsem_gene_tpm"            # log2(TPM + 0.001), TCGA + GTEx processed together
TOIL_PHENO = "TcgaTargetGTEX_phenotype.txt"
TOIL_PHENO_URL = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz"
PANCAN_HUB = "https://pancanatlas.xenahubs.net"
SURVIVAL = "Survival_SupplementalTable_S1_20171025_xena_sp"   # TCGA Clinical Data Resource (Liu et al. 2018)
SURVIVAL_URL = "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp"
TCGA_HUB = "https://tcga.xenahubs.net"
LIHC_EXPR = "TCGA.LIHC.sampleMap/HiSeqV2"                   # RSEM, log2(x+1), gene symbols
LIHC_EXPR_URL = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FHiSeqV2.gz"
LIHC_CNV = "TCGA.LIHC.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
LIHC_METH = "TCGA.LIHC.sampleMap/HumanMethylation450"

CBIOPORTAL = "https://www.cbioportal.org/api"
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{acc}.json"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
HALLMARK_LIBRARY = "MSigDB_Hallmark_2020"   # from Enrichr via gseapy

# TCGA primary site -> matching GTEx tissue where the two use different names (None = no fair GTEx match)
GTEX_SITE = {"Thyroid Gland": ("Thyroid", None), "Endometrium": ("Uterus", None), "Cervix": ("Cervix Uteri", None),
             "Adrenal gland": ("Adrenal Gland", None), "White blood cell": ("Blood", "Whole Blood"),
             "Lymphatic tissue": None, "Head and Neck region": None, "Soft tissue,Bone": None}

# Analysis choices, fixed before looking at results
MIN_GROUP = 10                 # minimum tumour and normal samples for a cancer to be tested
FDR = 0.05
N_BOOT = 1000                  # bootstrap replicates for AUC confidence intervals
N_PERM = 1000                  # permutations for the optimal-cutpoint test
CUTPOINT_RANGE = (0.2, 0.8)    # quantiles scanned by the "optimal cutpoint" approach
GSEA_PERM = 1000
SEED = 20250522
