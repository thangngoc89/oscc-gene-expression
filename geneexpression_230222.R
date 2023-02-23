# Set working directory
setwd("/Users/khoa/OneDrive - UMP/Đề tài Oral Cancer/11.Data/reproduce_Khoa")
options(stringsAsFactors = FALSE)

#CPU core
options(mc.cores = parallel::detectCores())

library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(plotly)
library(listviewer)
# Set color
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

#Load data
dds2 <- readRDS("./dds2.RDS")

# Plot count
gene <- c("TNFSF12", "TNFRSF25", "TNFRSF12A", "WNT5A", "CTHRC1", "SFRP2")

plot_gene_expression <- function(gene_name) {
  plot <- tryCatch(
    {
      plotCounts <- plotCounts(dds2, gene= gene_name, intgroup =c("cluster2","SampleID"), returnData=TRUE)
      plot <- ggplot(plotCounts, aes(x=cluster2, y=log(count), color =SampleID)) +
        geom_point(size =4) +
        labs(title = gene_name) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    },
    error=function(cond) {
      return(NA)
    },
    warning=function(cond) {
      return(NA)
    }
  )    
  
  return(plot)
}

for (i in 1:length(gene)) {
  plot <- plot_gene_expression(gene[[i]])
  if (any(is.na(plot))) {
    message(paste("Gene", gene[[i]], "not available"))
  }
  else {
    png(paste("plot_", gene[[i]], ".png", sep = ""))
    print(plot)
    dev.off()  
    message(paste("Extracted gene", gene[[i]]))
  }
}

plot <- plot_gene_expression("TNFSF12")
json <- plotly_json(plot, pretty = FALSE)
write(json[["x"]][["data"]], "TNFSF12.json")

# library(readr)
# readr::read_csv_chunked("genesymbol.csv", function(chunk) {
#   message(chunk)
# })
symbols <- read.csv("genesymbol.csv")
symbols <- symbols[,2]

library(RSQLite)
conn <- dbConnect(RSQLite::SQLite(), "geneexpression.db")
# Create toy data frames
car <- c('Camaro', 'California', 'Mustang', 'Explorer')
make <- c('Chevrolet','Ferrari','Ford','Ford')
df1 <- data.frame(car,make)
car <- c('Corolla', 'Lancer', 'Sportage', 'XE')
make <- c('Toyota','Mitsubishi','Kia','Jaguar')
df2 <- data.frame(car,make)
# Add them to a list
dfList <- list(df1,df2)
# Write a table by appending the data frames inside the list
for(k in 1:length(dfList)){
  dbWriteTable(conn,"Cars_and_Makes", dfList[[k]], append = TRUE)
}
conn <- dbConnect(RSQLite::SQLite(), "gene_expr.db")
for(i in 1:length(symbols)) {
  name <- c(symbols[[i]])
  plot <- plot_gene_expression(name)
  if (any(is.na(plot))) {
    expr <- NA
    df <- data.frame(name, expr)
    dbWriteTable(conn,"gene_expr", df, append = TRUE)
  }
  else {
    json <- plotly_json(plot, pretty = FALSE)
    expr <- toString(json[["x"]][["data"]])
    df <- data.frame(name, expr)
    dbWriteTable(conn,"gene_expr", df, append = TRUE)
  }
  message(i)
}

dbDisconnect(conn)