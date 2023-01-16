## to readout genelist in Szakasits_et_al-2009-The_Plant_Journal paper (Table S1)

library(pdftools)
library(stringr)
library(tabulizer)
library(dplyr)

setColbyWidth <- function(table, setCols, colwidth){
  w <- apply(table, 2, function(x) unname(unlist(sapply(x, nchar))))
  w <- apply(w, 2, max)
}
## to verifiy columns of the table
colValidate <- function(table, setCols, colwidth = NA) {
  
  wssub <- "##"
  if(is.null(colnames(table))) {
    if(any(setCols %in% table[1,])) {
      colnames(table) <- table[1,]
      table <- table[2:nrow(table),]
    } else if (ncol(table) == length(setCols)) {
      colnames(table) <- setCols
    } else {
      if(is.na(colwidth)) stop("please set column width")
      table <- setColbyWidth(table, setCols, colwidth)
    }
  }
  colnames(table) <- trimws(colnames(table), which="both")
  descol <- which(grepl("Description", colnames(table)))
  table[, descol] <- gsub(" ", wssub, table[,descol])
  table[, descol] <- sub("(\\d{1,5})##(AT[1-9,M,C]G\\d{5}.*)", 
                         "\\1 \\2 ", 
                         table[, descol],
                         ignore.case = T)
  nongenerow <- !grepl("\\d{1,5} AT[1-9,M,C]G\\d{5}", table[, descol], ignore.case = T)
  table[nongenerow, descol] <- sub("^", "0 ", table[nongenerow, descol])
  
  df <- data.frame(matrix(ncol = length(setCols), nrow = nrow(table)))
  colnames(df) <- setCols
  
  for(n in setCols) {
    if(n %in% colnames(table)) {
       df[,n] <- table[,n]
    } else {
       ninx <- which(grepl(n, colnames(table), ignore.case = F))
       
       if(length(ninx) == 0) {
         df[,n] <- NA
         next
       } else {
         orig.colnames <- trimws(unlist(str_split(colnames(table)[ninx], " ")), which = "both")
         npos <- which(n == orig.colnames)
         
         dat <- sapply(table[, ninx], function(x) unlist(str_split(x, " "))[npos])
         df[,n] <- unname(dat)
       }
       
    }
  }
  
 df$Description <- gsub(wssub, " ", df$Description)
 return(df)
}

desConsolidate <- function(table, delim = " ") {
  
  pointer <- 0
  for(r in 1:nrow(table)) {
    if(table[r, "Rank"] > 0) pointer <- r
    else {
      if(pointer == 0) next
      table[pointer, "Description"] <- paste(table[pointer, "Description"], 
                                          table[r, "Description"], sep = delim)
    }
  }
  table <- filter(table, Rank > 0)
  
  return(table)
}

getGeneID <- function(string) {
  geneid.pattern <- "AT[1-9,M,C]G\\d{5}"
  string <- trimws(string, which = "both")
  geneID <- sub(sprintf("^(%s)\\.\\d.*", geneid.pattern), "\\1", string, ignore.case = T)
  geneID <- trimws(geneID, which="both")
  return(geneID)
}

getGeneName <- function(strings) {
  
  sym <- str_extract(strings, "Symbol:.*?\\|")
  alias <- str_extract(strings, "Aliases:.*")
  sym <- paste(sym, alias, sep = ",")
  sym <- gsub("/|;", ",", sym)
  sym <- gsub("Symbol:|Aliases:|NA|\\|", "", sym)
  sym <- unname(unlist(sapply(sym, function(x) unlist(str_split(x, ","))[1])))
  sym <- trimws(sym, which = "both")
  return(sym)
  
}

## main
proj.path <- "syncytia2009"
#list.pdf <- file.path(proj.path, "ref", "TableS1.Syncytium.vs.Root(all).pdf")
ts5 <- file.path(proj.path, "ref", "TableS5.HighlyExpressed.in.Syncytium.pdf") 
pages <- pdf_info(ts5)$pages

knownCols <- c("Rank", "Description", "Sync", "Root", "M", "t")#, "adj.q", "B")
num.cols <- c("Sync", "Root", "M", "t")#, "adj.q", "B")
results <- list()
for(p in 2:pages) {
  cat("processing page", p, "\n")
  label <- sprintf("page_%d", p)
  tables <- extract_tables(ts5, pages = p)
  tables <- lapply(tables, function(t) colValidate(t, knownCols))
  table <- bind_rows(tables, .id = "table_in_page")
  table <- desConsolidate(table, delim = " ")
  results[[label]] <- table
}

result <- bind_rows(results, .id = "page")
ppct <- which(grepl("%", result$adj.q))
result[ppct, "adj.q"] <- sub("%", "", result[ppct, "adj.q"])
result[ppct, "adj.q"] <- as.numeric(result[ppct, "adj.q"])/100
result$adj.q <- as.numeric(result$adj.q)
result[num.cols] <- sapply(result[num.cols], as.numeric)

result$GeneID <- getGeneID(result$Description)
result$GeneName <- getGeneName(result$Description)
result$Description <- gsub(",", ";", result$Description)

## write the result table to file
write.csv(result, file.path(proj.path, "tableS5.csv"), quote = F)
## to read the result
res2 <- read.csv(file.path(proj.path, "tableS5.csv"), 
                 header = T, row.names = 1, stringsAsFactors = F)

