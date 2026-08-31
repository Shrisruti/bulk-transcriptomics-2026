#install required packages for RNA-seq tutorial

packages <- c("TCGAbiolinks","SummarizedExperiment","tibble","ggplot2","dplyr","DESeq2")

installed <- rownames(installed.packages())

for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(DESeq2)

# ==========================================================
# Install required packages for the scRNA-seq tutorial
# ==========================================================

cran_packages <- c(
  "Seurat",
  "patchwork",
  "ggplot2",
  "dplyr"
)

# Install missing packages
installed <- rownames(installed.packages())

for (pkg in cran_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load packages
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)

cat("\n")
cat("=========================================\n")
cat("All required packages are installed.\n")
cat("You are ready to run the tutorial.\n")
cat("=========================================\n")
