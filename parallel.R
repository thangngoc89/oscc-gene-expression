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
  
  conn <- dbConnect(RSQLite::SQLite(), paste("gene_expr-", core, ".db", sep=""))
  
  is_interactive <- interactive()
  
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
      if (is_interactive) {
        expr <- toString(json[["x"]][["data"]])  
      }
      else {
        expr <- toString(json)  
      }

      df <- data.frame(name, expr)
      dbWriteTable(conn,"gene_expr", df, append = TRUE)
    }
    message(paste(core,i))
  }
  dbDisconnect(conn)
}

