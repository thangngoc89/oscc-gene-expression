library(DESeq2)
library(jsonlite)

dds2 <- readRDS("./dds2.RDS")

gene_name <- "FER1L6"
# Count
plotCounts <- plotCounts(dds2, gene= gene_name, intgroup =c("cluster2","SampleID"), returnData=TRUE)
# log(count)
plotCounts$count <- log(plotCounts$count)
# Turn cluster2 levels into integer
plotCounts$cluster2 <- as.numeric(levels(plotCounts$cluster2))
# Rename columns (for smaller JSON size)
colnames(plotCounts) <- c('y','x','n')
json <-toJSON(plotCounts, rownames = FALSE, digits=NA)

message(json)


