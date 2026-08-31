#downloading transcriptome data from TCGA - HNSC firehose legacy
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)

clin <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
table(clin$tissue_or_organ_of_origin) 

#"Tongue, NOS", "Base of tongue, NOS", "Border of tongue","Ventral surface of tongue, NOS
tongue_samples <- clin[grep("tongue", clin$tissue_or_organ_of_origin, ignore.case = TRUE),]
table(tongue_samples$tissue_or_organ_of_origin)

query <- GDCquery(project = "TCGA-HNSC",data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",barcode = tongue_samples$submitter_id,
  workflow.type = "STAR - Counts",access = "open")

sample_info <- getResults(query)
nrow(sample_info)
GDCdownload(query = query, method = "api", files.per.chunk = 20) 
data <- GDCprepare(query = query)

genes <- rowData(data)

if ("gene_type" %in% colnames(genes)) {
  protein_coding <- genes$gene_type == "protein_coding"
  data_clean <- data[protein_coding, ]
} else {
  data_clean <- data
}

counts  <- assay(data_clean, "unstranded")
samples <- colData(data_clean)
genes1  <- rowData(data_clean)

save(data, data_clean, counts, samples, genes1, 
     file = "tcga_hnsc_tongue_clean.RData")

table(samples$sample_type)  # check sample_type — Primary Tumor vs Solid Tissue Normal are the usual two groups
genes <- rowData(data)
table(samples$sample_type) #has solid tissue normal from tongue

library(DESeq2)
library(apeglm)

# DESeq2 needs the counts and colData in the same sample order — TCGAbiolinks
# preserves this automatically since both come from the same data_clean object
dds <- DESeqDataSetFromMatrix(countData = counts, colData   = samples,
  design    = ~sample_type)

# remove genes with essentially no expression across all samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# set the reference level explicitly so log2FC direction is unambiguous
dds$sample_type <- relevel(factor(dds$sample_type), ref = "Solid Tissue Normal")

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
barplot(sizeFactors(dds), las = 2, main = "Size factors per sample") 
# quick sanity check — size factors shouldn't vary wildly if libraries are comparable

dds <- estimateDispersions(dds)
plotDispEsts(dds)

dds <- nbinomWaldTest(dds)
resultsNames(dds)

dds <- DESeq(dds) #go through the message on the console and connect it to the lecture
res <- results(dds)

summary(res)

# shrink log2FC for reliable ranking (needed for GSEA later)
res_shrunk <- lfcShrink(dds, coef = "sample_type_Primary.Tumor_vs_Solid.Tissue.Normal",
                        type = "apeglm")

# tidy into a data frame, map to gene symbols
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(as.data.frame(genes1) %>% select(gene_id, gene_name), by = "gene_id") %>%
  arrange(padj)

deg <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

#how to get up and down regulated gene here?

write.csv(res_df, "deseq2_results_full.csv", row.names = FALSE)
write.csv(deg, "deseq2_significant_genes.csv", row.names = FALSE)

#ORA - clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)

# ORA needs Entrez IDs — map from gene_name (SYMBOL)
deg_entrez <- bitr(deg$gene_name, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db, drop = TRUE)
#GO database
ego <- enrichGO(gene = deg_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

dotplot(ego, showCategory = 20)

#KEGG database
ekegg <- enrichKEGG(gene = deg_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
dotplot(ekegg, showCategory = 20)

#GSEA - using all genes
# build ranked gene list from the FULL results table (res_df), not `deg`
gsea_input <- res_df %>%
  filter(!is.na(padj), !is.na(gene_name)) %>%
  distinct(gene_name, .keep_all = TRUE)

# map to Entrez for org.Hs.eg.db-based GSEA
gsea_entrez <- bitr(gsea_input$gene_name, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db, drop = TRUE)

gsea_input <- gsea_input %>%
  inner_join(gsea_entrez, by = c("gene_name" = "SYMBOL"))

ranked_genes <- gsea_input$log2FoldChange
names(ranked_genes) <- gsea_input$ENTREZID
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

gsea_go <- gseGO(geneList = ranked_genes, OrgDb = org.Hs.eg.db,
                 ont = "BP", keyType = "ENTREZID",
                 pvalueCutoff = 0.05, verbose = FALSE)
dotplot(gsea_go, showCategory = 20)
#ridgeplot(gsea_go)
enrichplot::gseaplot2(gsea_go, geneSetID = 1)   # enrichment plot for the top term

gsea_kegg <- gseKEGG(geneList = ranked_genes, organism = "hsa",
                     pvalueCutoff = 0.05, verbose = FALSE)

dotplot(gsea_kegg, showCategory = 20)
ridgeplot(gsea_kegg)
gseaplot2(gsea_kegg, geneSetID = 1)



