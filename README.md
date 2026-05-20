# TCGA-LUAD Differential Expression & Pathway Analysis

Bulk RNA-seq analysis comparing lung adenocarcinoma tumor vs. normal tissue.  
**60 samples** (30 tumor, 30 normal) from TCGA-LUAD · DESeq2 · fgsea Hallmark pathways

---

## Pipeline

```
TCGA download (barcodes filtered — targeted subset, not full dataset)
        ↓
QC (library size, PCA, sample distances)
        ↓
DESeq2 (raw counts, Wald test, tumor vs normal)
        ↓
fgsea (Hallmark pathways, ranked gene list)
```

---

## DESeq2 Results Summary

| | Count |
|---|---|
| Genes tested | 33,136 |
| Significant (padj < 0.05, \|log2FC\| > 1) | 9,641 |
| Upregulated in tumor | 6,193 |
| Downregulated in tumor | 3,448 |
| Outliers | 0 |

---

## Key Pathway Results (fgsea — Hallmark)

### Upregulated in Tumor (NES > 0)

| Pathway | NES | padj |
|---|---|---|
| HALLMARK_G2M_CHECKPOINT | +1.89 | 4.1e-12 |
| HALLMARK_E2F_TARGETS | +1.87 | 9.3e-12 |
| HALLMARK_GLYCOLYSIS | +1.66 | 9.4e-05 |
| HALLMARK_MITOTIC_SPINDLE | +1.64 | 1.6e-04 |
| HALLMARK_MYC_TARGETS_V1 | +1.55 | 2.3e-03 |
| HALLMARK_PI3K_AKT_MTOR_SIGNALING | +1.54 | 5.9e-03 |
| HALLMARK_SPERMATOGENESIS | +1.49 | 1.6e-02 |
| HALLMARK_MTORC1_SIGNALING | +1.41 | 3.6e-02 |
| HALLMARK_MYC_TARGETS_V2 | +1.48 | 4.9e-02 |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | +1.34 | 1.0e-01 |
| HALLMARK_DNA_REPAIR | +1.31 | 1.1e-01 |

### Downregulated in Tumor (NES < 0)

| Pathway | NES | padj |
|---|---|---|
| HALLMARK_ADIPOGENESIS | −1.85 | 6.3e-06 |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | −1.81 | 1.2e-05 |
| HALLMARK_INFLAMMATORY_RESPONSE | −1.70 | 2.5e-04 |
| HALLMARK_FATTY_ACID_METABOLISM | −1.51 | 8.2e-03 |
| HALLMARK_TGF_BETA_SIGNALING | −1.50 | 6.7e-02 |
| HALLMARK_HEME_METABOLISM | −1.50 | 9.5e-03 |
| HALLMARK_MYOGENESIS | −1.39 | 4.2e-02 |
| HALLMARK_CHOLESTEROL_HOMEOSTASIS | −1.40 | 1.0e-01 |
| HALLMARK_IL6_JAK_STAT3_SIGNALING | −1.30 | 1.6e-01 |

---

## Biological Interpretation

**Upregulated — Proliferation programs active:**
- E2F + G2M + Mitotic Spindle → cell cycle dysregulation, uncontrolled division
- MYC targets → oncogenic growth program active
- Glycolysis → Warburg effect (cancer cells shift to aerobic glycolysis)
- PI3K/AKT/mTOR → survival and growth signaling constitutively ON

**Downregulated — Normal lung identity lost:**
- Adipogenesis → stromal/tissue identity lost
- TNFA/NFkB + Inflammatory Response → immune signaling suppressed in tumor
- Fatty acid metabolism → normal lipid metabolism gone
- TGF-beta → complex role; suppression may reflect immune evasion

---

## LUAD vs BRCA — Comparison

| Pathway | LUAD | BRCA |
|---|---|---|
| E2F_TARGETS | ✅ Up | ✅ Up |
| G2M_CHECKPOINT | ✅ Up | ✅ Up |
| MYC_TARGETS | ✅ Up | ✅ Up |
| MTORC1_SIGNALING | ✅ Up | ✅ Up |
| GLYCOLYSIS | ✅ Up (stronger) | — |
| PI3K_AKT_MTOR | ✅ Up | — |
| ADIPOGENESIS | ✅ Down | ✅ Down |
| FATTY_ACID_METABOLISM | ✅ Down | ✅ Down |
| MYOGENESIS | ✅ Down | ✅ Down |
| TNFA/NFkB | ✅ Down (LUAD) | ✅ Up (BRCA) |
| INFLAMMATORY_RESPONSE | ✅ Down (LUAD) | — |

**Key difference:** TNFA/NFkB and Inflammatory Response are **downregulated in LUAD** but were **upregulated in BRCA** — reflecting different immune microenvironments between lung and breast cancer.

**Shared across both cancers:** Proliferation programs (E2F, G2M, MYC) are universally active — classic hallmarks of cancer regardless of tissue origin.

---

## Requirements

```r
BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "fgsea"
))

install.packages(c("msigdbr", "data.table"))
```

---

## Usage

```r
source("luad.R")
```

- `set.seed(42)` used throughout for reproducibility
- Raw counts (STAR - Counts, unstranded) — not TPM/FPKM
- Reference level: Solid Tissue Normal
- Filter: genes with ≥ 10 counts in at least 1 sample
- fgsea ranking: `log2FC × −log10(pvalue + 1e−300)`

---

## Project Structure

```
tcga-luad-deseq2/
├── luad.R               # full pipeline
├── README.md
└── .gitignore
```

---

