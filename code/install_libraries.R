install.packages(c("workflowr", "dplyr", "ggplot2", "rmarkdown", "DT"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer", "airway"))
