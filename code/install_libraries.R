install.packages(c("workflowr", "dplyr", "ggplot2", "rmarkdown", "DT"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer", "airway"))


install.packages(c("ggfortify", "ggpubr", "ggplot2", "ggsci", "corrplot", "sva", "pheatmap"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c( "GEOquery", "DESeq2", "Glimma", "AnnotationDbi", "org.Hs.eg.db",
                        "clusterProfiler", "DOSE", "enrichplot", "topGO", "GeneTonic" ))


library(GEOquery)
library(DESeq2)
library(ggfortify)
library(ggpubr)
library(ggsci)
library(Glimma)
library(corrplot)
library(sva)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("clusterProfiler")
library(pheatmap)
library(DOSE)
library(enrichplot)
library(topGO)
library(GeneTonic)
