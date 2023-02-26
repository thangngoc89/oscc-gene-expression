library(RSQLite)
library(parallel)

symbols <- read.csv("genesymbol.csv")
symbols <- length(symbols[[1]])

if (file.exists("merged.db")) {
  file.remove("merged.db")
}

conn <- dbConnect(RSQLite::SQLite(), "merged.db")

total_cores <- parallel::detectCores()

get_temp_db_name <- function(core) {
  return(paste("parallel-temp/gene_expr-", core, ".db", sep=""))
}

for(core in 1:total_cores) {
  core_conn <- dbConnect(RSQLite::SQLite(), get_temp_db_name(core))
  core_result <- dbReadTable(core_conn, "gene_expr")
  dbWriteTable(conn, "gene_expr", core_result, append = TRUE)
  dbDisconnect(core_conn)
  message(paste("Done", core))
}

total_rows <- dbGetQuery(conn, "SELECT count(name) as length FROM gene_expr")
if (total_rows$length != symbols) {
  message("Total rows not equal")
  message(paste("Symbols:", symbols))
  message(paste("Total rows:", total_rows$length))
} else {
  message("All checked !")
  message("Remove temporary databases")
  
  for(core in 1:total_cores) {
    db_name <- get_temp_db_name(core)
    if (file.exists(db_name)) {
      file.remove(db_name)
    }
  }
}

dbDisconnect(conn)