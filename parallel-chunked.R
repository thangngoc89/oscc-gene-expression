library(doParallel)
library(parallel)

options(stringsAsFactors = FALSE)

total_cores <- parallel::detectCores()
registerDoParallel(total_cores)

foreach(core=1:total_cores, .combine = c) %dopar% {
  library(RColorBrewer)
  library(DESeq2)
  library(ggplot2)
  library(plotly)
  library(listviewer)
  library(RSQLite)
  
  symbols <- read.csv("genesymbol.csv")
  symbols <- symbols[,2]
  
  chunked <- split(symbols,             
                   cut(seq_along(symbols),
                       total_cores,
                       labels = FALSE))
  symbols <- chunked[[core]]
  
  # Set color
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  heat_colors <- brewer.pal(n = 6, name = "YlOrRd")
  
  #Load data
  dds2 <- readRDS("./dds2.RDS")
  
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
  
  json_gene_expression <- function (gene_name) {
    plot <- plot_gene_expression(gene_name)
    if (any(is.na(plot))) {
      return(NA)
    }
    else {
      json <- plotly_json(plot, pretty = FALSE)
      if (interactive()) {
        return(toString(json[["x"]][["data"]]))
      }
      
      return(toString(json))
    }
  }
  
  
  conn <- dbConnect(RSQLite::SQLite(), paste("gene_expr-", core, ".db", sep=""))
  
  
  
  # Batch writes into 100 genes
  
  chunk_length <- 20
  chunked <- split(symbols,             # Applying split() function
                   ceiling(seq_along(symbols) / chunk_length))
  
  for(i in 1:length(chunked)) {
    name <- as.vector(chunked[[i]])
    expr <- lapply(name, json_gene_expression)
    expr <- as.vector(expr)
    nested <- list(name=name, expr=expr)
    df <-  as.data.frame(do.call(cbind, nested))
    
    dbWriteTable(conn,"gene_expr", df, append = TRUE)
    message(paste(core,i*chunk_length))
  }
  dbDisconnect(conn)
}

