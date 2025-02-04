library(doParallel)
library(parallel)

options(stringsAsFactors = FALSE)

total_cores <- parallel::detectCores()
registerDoParallel(total_cores)

dir.create("parallel-temp")

foreach(core=1:total_cores, .combine = c) %dopar% {
  library(DESeq2)
  library(jsonlite)
  library(RSQLite)
  
  symbols <- read.csv("genesymbol.csv")
  symbols <- symbols[,2]
  
  chunked <- split(symbols,             
                   cut(seq_along(symbols),
                       total_cores,
                       labels = FALSE))
  symbols <- chunked[[core]]
  
  #Load data
  dds2 <- readRDS("./dds2.RDS")
  
  plot_gene_expression <- function(gene_name) {
    plot <- tryCatch(
      {
        # Count
        plotCounts <- plotCounts(dds2, gene= gene_name, intgroup =c("cluster2","SampleID"), returnData=TRUE)
        # log(count)
        plotCounts$count <- log(plotCounts$count)
        # Turn cluster2 levels into integer
        # plotCounts$cluster2 <- as.numeric(levels(plotCounts$cluster2))
        # Rename columns (for smaller JSON size)
        colnames(plotCounts) <- c('y','x','n')
        json <-toJSON(plotCounts, rownames = FALSE, digits=NA)
        return(json)
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
  
  conn <- dbConnect(RSQLite::SQLite(), paste("parallel-temp/gene_expr-", core, ".db", sep=""))
  
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
      expr <- toString(plot)  
      df <- data.frame(name, expr)
      dbWriteTable(conn,"gene_expr", df, append = TRUE)
    }
    message(paste(core,i))
  }
  dbDisconnect(conn)
}

source("./merge-db.R")
unlink("parallel-temp")