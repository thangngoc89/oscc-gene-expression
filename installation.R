if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("plotly", quietly = TRUE))
  install.packages("plotly")

if (!require("listviewer", quietly = TRUE))
  install.packages("listviewer")

if (!require("RSQLite", quietly = TRUE))
  install.packages("RSQLite")