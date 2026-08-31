#install required packages for RNA-seq tutorial

packages <- c("TCGAbiolinks","SummarizedExperiment","tibble",
              "ggplot2","dplyr","DESeq2", "clusterProfiler")

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
library(ggplot2)
library(clusterProfiler)
