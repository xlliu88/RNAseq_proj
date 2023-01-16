  
### import Results

importResults.old <- function(resultfiles, pattern = NA) {
  result <- list()
  
  if(!is.na(pattern[1])) {
    compare.names <- compare.names[sapply(pattern, function(x) which(grepl(pattern, names(compare.names))))]
  }
  
 
  for (name in names(compare.names)) {
    if(name %in% names(resultfiles)){
      file_name <- resultfiles[[name]]
      sep <- ifelse(grepl(".txt$", file_name), "\t", ",")
      res <- read.table(file.path(resultfiles$path, resultfiles[[name]]), header = T, sep = sep, row.names = 1)
      res$Gene_id <- row.names(res)
      res$padj <- ifelse(is.na(res$padj), 1, res$padj)
      row.names(res) <- res$Gene_id
      res <- res[order(row.names(res)),]
      result[[name]] <- res
    }
  }
  return(result)
}

importResults <- function(resultfiles, pattern = NA) {
  result <- list()
  compare.names <- compare.names[names(compare.names) %in% names(resultfiles)]
  if(!is.na(pattern[1])) {
    compare.names <- compare.names[sapply(pattern, function(x) {
      which(grepl(pattern, names(compare.names)))
      })]
  }
  
  result <- lapply(names(compare.names), function(x) {
    file_name <- resultfiles[[x]]
    sep <- ifelse(grepl(".txt$", file_name), "\t", ",")
    qot <- "\"" #ifelse(grepl(".txt$", file_name), "\"", "")
    res <- read_delim(file.path(resultfiles$path, resultfiles[[x]]),
                      col_types = cols(),
                      delim = sep, 
                      quote = qot)  %>% 
      mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
      mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue)) %>%
      mutate(Gene_name = ifelse(is.na(Gene_name), "", Gene_name)) %>% 
      rename("Gene_id" = `...1`) %>%
      arrange(Gene_id) %>%
      dplyr::select(Gene_id, Gene_name, everything())
    return(res)
  })
  names(result) <- names(compare.names)
  return(result)
}

##
### consolidate LCM BCN comparisons to one dataframe
consolidateData <- function(data, treatment = c("BCN", "CLE")){
    if(treatment == "BCN") {
      trts <- c("wbc", "vbc", "vwc", "vwb", "clvBCN")
    } else if (treatment == "CLE") {
      trts <- c("wec", "vec", "vwc", "vwe", "clvCLE")
    } else {
      stop("wrong treatment type, should be one of (BCN, CLE)")
    }
    
    nr <- nrow(data[[1]])
    res <- data.frame(matrix(nrow=nr, ncol=0))
    for (n in trts) {
      df <- data[[n]]
      geneInfo <- select(df, Gene_name, Gene_type, Gene_description)
      df <- select(df, log2FoldChange, padj)
      df <- round(df, 4)
      colnames(df) <- paste(c("lfc", "p"), n, sep = ".")
      res <- cbind(res, df)
    }
    
    res$dLFC <- res[[paste("lfc", trts[1], sep = ".")]] - res[[paste("lfc", trts[2], sep = ".")]]
    res$Gene_id <- row.names(res)
    res <- cbind(res, geneInfo)
    res <- res[order(res$dLFC),]
    
    return(res)
}   
