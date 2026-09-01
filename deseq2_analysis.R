############ DO NOT RUN THE CODE - START ###############
#load the necessary libraries
library(TCGAbiolinks)
library(SummarizedExperiment)

#downloading transcriptome data from TCGA - HNSC firehose legacy
clin <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
table(clin$tissue_or_organ_of_origin) 
#Result - "Tongue, NOS", "Base of tongue, NOS", "Border of tongue","Ventral surface of tongue, NOS

#store the sample ids of all the patients with tongue as the tissue of origin
tongue_samples <- clin[grep("tongue", clin$tissue_or_organ_of_origin, ignore.case = TRUE),]
table(tongue_samples$tissue_or_organ_of_origin)

#fetch the rnaseq data for the selected patient ids
query <- GDCquery(project = "TCGA-HNSC",data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",barcode = tongue_samples$submitter_id,
                  workflow.type = "STAR - Counts",access = "open")
sample_info <- getResults(query)
nrow(sample_info)
GDCdownload(query = query, method = "api", files.per.chunk = 20) 
data <- GDCprepare(query = query)
genes <- rowData(data)

#select and retain only protein coding genes
if ("gene_type" %in% colnames(genes)) {
  protein_coding <- genes$gene_type == "protein_coding"
  data_clean <- data[protein_coding, ]
} else {
  data_clean <- data
}
#from the assaydata store the counts and sample metadata
counts  <- assay(data_clean, "unstranded")
samples <- colData(data_clean)
genes1  <- rowData(data_clean)

save(counts, file = "data/tcga_hnsc_counts.RData")
save(samples, file = "data/tcga_hnsc_samples.RData")
save(genes1, file = "data/tcga_hnsc_genes.RData")

############ DO NOT RUN THE CODE - END ###############
#load necessary packages

library(dplyr)
library(tibble)
library(DESeq2)
library(apeglm)

#load the necessary data files
load("data/tcga_hnsc_counts.RData")
load("data/tcga_hnsc_samples.RData")
load("data/tcga_hnsc_genes.RData")

table(samples$sample_type)  # check sample_type - Primary Tumor vs Solid Tissue Normal are the usual two groups

#DESeq2 needs the counts and colData in the same sample order — TCGAbiolinks
#preserves this automatically since both come from the same data_clean object
#create the S4 DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples,
                              design = ~sample_type)

#remove genes with essentially no expression across all samples
#exercise 1: Explain what is stored in the variable 'keep'
keep <- rowSums(counts(dds) >= 10) >= 3 
#exercise 2: How many genes have been dropped?
dds <- dds[keep, ] 

#set the reference level explicitly
dds$sample_type <- relevel(factor(dds$sample_type), ref = "Solid Tissue Normal")

#estimate the size factors
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
barplot(sizeFactors(dds), las = 2, main = "Size factors per sample") 

#estimate dispersions
dds <- estimateDispersions(dds)
plotDispEsts(dds) #exercise 3: Explain the plot

#QC: PCA to check sample segregation
# variance-stabilizing transformation for QC — blind=TRUE ignores the design formula 
#the PCA reflects raw structure in the data
vsd <- vst(dds, blind = TRUE)

# quick built-in PCA plot
plotPCA(vsd, intgroup = "sample_type")
#exercise 4: Do tumor and normal samples separate cleanly along PC1? 

dds <- nbinomWaldTest(dds)
resultsNames(dds)

#exercise 5: Go through the message on the console and connect it to the lecture
dds <- DESeq(dds) 
res <- results(dds)
summary(res)

#shrink log2FC for reliable ranking (needed for GSEA later)
#if you have problems with apeglm library - do not worry and skip this step
res_shrunk <- lfcShrink(dds, coef = "sample_type_Primary.Tumor_vs_Solid.Tissue.Normal",
                        type = "apeglm")

#tidy into a data frame, map to gene symbols
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(as.data.frame(genes1) %>% select(gene_id, gene_name), by = "gene_id") %>%
  arrange(padj)

deg <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
#exercise 5: How to get up and down regulated genes here?

write.csv(res_df, "deseq2_results_full.csv", row.names = FALSE)
write.csv(deg, "deseq2_significant_genes.csv", row.names = FALSE)

############### Enrichment analysis ###############
#Over representation analysis - clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ORA needs Entrez IDs - map from gene_name (SYMBOL)
deg_entrez <- bitr(deg$gene_name, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db, drop = TRUE)
# use GO-BP database as reference
ego <- enrichGO(gene = deg_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

dotplot(ego, showCategory = 20)

#use KEGG database as reference
ekegg <- enrichKEGG(gene = deg_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
dotplot(ekegg, showCategory = 20)

# GSEA - using all genes
# build ranked gene list from the FULL results table (res_df), not deg
gsea_input <- res_df %>% filter(!is.na(padj), !is.na(gene_name)) %>%
  distinct(gene_name, .keep_all = TRUE)

# map to Entrez for org.Hs.eg.db-based GSEA
gsea_entrez <- bitr(gsea_input$gene_name, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db, drop = TRUE)

gsea_input <- gsea_input %>% inner_join(gsea_entrez, by = c("gene_name" = "SYMBOL"))

ranked_genes <- gsea_input$log2FoldChange
names(ranked_genes) <- gsea_input$ENTREZID
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

#GO database
gsea_go <- gseGO(geneList = ranked_genes, OrgDb = org.Hs.eg.db,
                 ont = "BP", keyType = "ENTREZID",
                 pvalueCutoff = 0.05, verbose = FALSE)
#dotplot(gsea_go, showCategory = 20)
dotplot(gsea_go, showCategory = 20, split = ".sign") + 
  facet_grid(. ~ .sign)   # "activated" = NES > 0, "suppressed" = NES < 0
ridgeplot(gsea_go) + labs(x = "log2 fold change distribution")

enrichplot::gseaplot2(gsea_go, geneSetID = 1)   # enrichment plot for the top term
# top 3 by adjusted p-value
top_ids <- gsea_go@result %>% arrange(p.adjust) %>% slice_head(n = 3) %>% pull(ID)
enrichplot::gseaplot2(gsea_go, geneSetID = top_ids, pvalue_table = TRUE)

#KEGG database
gsea_kegg <- gseKEGG(geneList = ranked_genes, organism = "hsa",
                     pvalueCutoff = 0.05, verbose = FALSE)

#dotplot(gsea_kegg, showCategory = 20)
dotplot(gsea_go, showCategory = 20, split = ".sign") + 
  facet_grid(. ~ .sign)   # "activated" = NES > 0, "suppressed" = NES < 0
ridgeplot(gsea_go) + labs(x = "log2 fold change distribution")

enrichplot::gseaplot2(gsea_kegg, geneSetID = 1)   # enrichment plot for the top term
# top 3 by adjusted p-value
top_ids <- gsea_kegg@result %>% arrange(p.adjust) %>% slice_head(n = 3) %>% pull(ID)
enrichplot::gseaplot2(gsea_kegg, geneSetID = top_ids, pvalue_table = TRUE)
