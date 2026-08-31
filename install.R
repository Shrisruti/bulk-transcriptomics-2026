#install required packages for RNA-seq tutorial

# Install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")}

# Bioconductor packages
bioc_packages <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", "enrichplot", "org.Hs.eg.db", "clusterProfiler")

# CRAN packages
cran_packages <- c("tibble", "ggplot2", "dplyr")

installed <- rownames(installed.packages())

for (pkg in bioc_packages) {
  if (!pkg %in% installed) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

for (pkg in cran_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}
