library(TCGAbiolinks)
library(SummarizedExperiment)

R.version$version.string

#Tcga data query 
query= GDCquery(project = "TCGA-LUAD",
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                workflow.type = "STAR - Counts",
                sample.type = c("Primary Tumor", "Solid Tissue Normal")
)


metadata=getResults(query)

colnames(metadata)
table(metadata$sample_type)

tumor_lung <- metadata[metadata$sample_type =="Primary Tumor", ]
normal_lung <- metadata[metadata$sample_type =="Solid Tissue Normal", ]


nrow(tumor_lung)
nrow(normal_lung)

set.seed(42)
tumor_30 <- tumor_lung[sample(nrow(tumor_lung), 30), ]
normal_30 <- normal_lung[sample(nrow(normal_lung), 30), ]

#combine both 
data <- rbind(tumor_30, normal_30)
selected_data <- data$cases
head(selected_data)

query_60 <- GDCquery(project = "TCGA-LUAD",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts",
                     sample.type = c("Primary Tumor", "Solid Tissue Normal"), 
                     barcode=selected_data) #filtering here at this step

GDCdownload(query_60, method="api", files.per.chunk=10)

se<- GDCprepare(query_60)

#counts
count_mat <- assay(se, "unstranded")
#extracting coldata
sample_info <- data.frame(sample_type=se$sample_type, row.names = colnames(se))
#rowinfo
gene_info <- rowData(se)

sample_info$sample_type <- factor(
  sample_info$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor")
)  
  
  #DESEQDataSet
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = sample_info, design=~sample_type)

#filter low counts
keep <- rowSums(counts(dds)>=10) >=1
dds <- dds[keep, ]

#Rundeseq2
dds <- DESeq2(dds)

res <- results(dds,
               contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"),
               alpha = 0.05
)


summary(res)
res_clean <- res[!is.na(res$padj), ]
nrow(res_clean)  # kitne genes bache

# Significant = padj < 0.05 AND |log2FC| > 1
res_sig <- res_clean[
  res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1,
]

nrow(res_sig)

# Tumor more expressed
upregulated <- res_sig[res_sig$log2FoldChange > 1, ]
upregulated <- upregulated[order(upregulated$log2FoldChange, decreasing = TRUE), ]

# Tumor less expressed  
downregulated <- res_sig[res_sig$log2FoldChange < -1, ]
downregulated <- downregulated[order(downregulated$log2FoldChange), ]

# Count
nrow(upregulated)
nrow(downregulated)

# Top genes dekho
head(upregulated, 10)
head(downregulated, 10)

colnames(gene_info) 
#mapping IDs to names
id_to_symbol <- setNames(
  gene_info$gene_name,
  rownames(gene_info)
)
dim(res_clean)
rownames(res_clean)
# Adding symbols to Res_clean 
res_clean$symbol <- id_to_symbol[rownames(res_clean)]

head(res_clean[, c("log2FoldChange", "padj", "symbol")])

#Pathway genes
library(fgsea)
library(msigdbr)

ranked_genes <- res_clean$log2FoldChange * -log10(res_clean$pvalue + 1e-300)

names(ranked_genes) <- res_clean$symbol

# removing NA symbols
ranked_genes <- ranked_genes[!is.na(names(ranked_genes))]

# Sorting necessary for fgsea
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Verify
head(ranked_genes, 5)
tail(ranked_genes, 5)
length(ranked_genes)


# Hallmark pathways
pathways_h <- msigdbr(
  species = "Homo sapiens",
  collection = "H"          # H = Hallmark
)

# fgsea format conversion— we need named list, not long format
pathway_list <- split(
  pathways_h$gene_symbol,
  pathways_h$gs_name
)

sum(duplicated(names(ranked_genes)))

# Duplicates removing 
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

# Verify
sum(duplicated(names(ranked_genes)))  # should turn out to be 0 


# Step 1: removing Duplicates 
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

# Step 2: fgsea
set.seed(42)
fgsea_res <- fgsea(
  pathways = pathway_list,
  stats    = ranked_genes,
  minSize  = 15,
  maxSize  = 500
)

# Step 3: Results 
head(fgsea_res[order(fgsea_res$padj), ], 20)


# Significant pathways
sig_pathways <- fgsea_res[fgsea_res$padj < 0.05, ]
sig_pathways <- sig_pathways[order(sig_pathways$NES, decreasing = TRUE), ]

# Upregulated pathways (NES > 0 = enriched in Tumor)
up_pathways   <- sig_pathways[sig_pathways$NES > 0, ]

# Downregulated pathways (NES < 0 = depleted in Tumor)
down_pathways <- sig_pathways[sig_pathways$NES < 0, ]

nrow(up_pathways)
nrow(down_pathways)

colnames(se)