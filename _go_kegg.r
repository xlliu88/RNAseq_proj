## general functions
## take a gene list and ontology of interest
## genelist should be a named vector of p values
## whichGO is one of "BP", "CC", "MF"
## three steps to get statistical result
##  build GO Object
##  runTest (to get topGOresult object) 
##  GenTable  to get statistical result
## return a data frame with GO terms and p-value
## if not out.path==NA, will out put a graph in pdf format

## outFormat: all, combined
goTair <- function(genes, whichGO,
                   pthreshold = 0.05, 
                   nodeSize = 10, 
                   topNodes = NA, 
                   algorithm = "weight01", 
                   stat = c("Fisher", "KS"), 
                   outFormat = "combined", 
                   combMethod = "gmean",
                   out.path = NA, 
                   fn.prefix = NA) {
  

  
  fn <- sprintf("%s.%s.%s.%s", fn.prefix, whichGO, 
                paste0(algorithm, collapse = "-"), 
                paste0(stat, collapse = "-"))
  
  GOobj <- new("topGOdata",
               ontology = whichGO,
               allGenes = genes,
               geneSel = function(x) x < pthreshold ,
               nodeSize = nodeSize, 
               annot = annFUN.org,
               mapping = "org.At.tair.db")
  
  par <- expand.grid(algorithm, stat, stringsAsFactors = F)
  
  result <- lapply(1:nrow(par), 
                   function(i) runTest(GOobj, 
                                       algorithm = par[i, 1], 
                                       statistic = par[i, 2]))
  names(result) <- apply(par, 1, function(x) str_c(x, collapse = ""))

  if(nrow(par) == 1) outFormat <- "all"
  if(tolower(outFormat) == "all") {
      n <- length(score(result[[1]]))
      res.table <- do.call("GenTable", c(list(object = GOobj),
                                         result,
                                         list(topNodes = n))) %>% 
        mutate(Term = str_replace(Term, ",", ";"))
      
      if(nrow(par) == 1) colnames(res.table)[ncol(res.table)] <- "pvalue"
      
      if (!is.na(out.path)) {
        for(i in 1:length(result)) {
          graph.name <- file.path(out.path, paste0(fn.prefix, whichGO, names(result)[i], sep = "."))
          sigGO <- sum(score(result[[i]]) < 0.01)
          topNodes <- ifelse(is.na(topNodes), sigGO, topNodes)
          printGraph(GOobj, result[[i]], firstSigNodes = topNodes, 
                     fn.prefix = graph.name, pdfSW=T, useInfo="all")
        }
        
        suffix <- ifelse(length(result) == 1, "", "all")   
        table.name <- sprintf("%s.%s.csv", fn, suffix)
        write_csv(res.table, file.path(out.path, table.name))
     }
      
  } else if (tolower(outFormat) == "combined") {
      res <- combineGoResultsInList(result, method = combMethod)
      res.table <- GenTable(GOobj, pvalue = res, topNodes = length(score(res))) %>% 
        arrange(pvalue) %>% 
        mutate(Term = str_replace_all(Term, ",", "; "))

      if(!is.na(out.path)) {
          sigGO <- sum(score(res) < 0.01)
          graph.name <- file.path(out.path, paste0(fn, "combined", sep = "."))
          topNodes <- ifelse(is.na(topNodes), sigGO, topNodes)
          printGraph(GOobj, res, firstSigNodes = topNodes, 
                     fn.prefix = graph.name, 
                     pdfSW=T, useInfo="all")
          
          table.name <- paste(fn, "combined.csv", sep = ".")
          write_csv(res.table, file.path(out.path, table.name))
      }
  } 
  return(res.table)
}


combineGoResultsInList <- function(result, method = "gmean") {
    
    ## takes a list of topGOresult objects and combined them to one dataframe.
    ## method = c("gmean", "min", "median")
    
     if(!class(result)=="list") {
       stop("Input error: only takes list of topGOresult objects as Input!")
     }

     if(length(result) == 1) {
        return(result[[1]])
     } else {
        res <- combineResults(result[[1]], result[[2]], method = method)
        result[[1]] <- res
        result[[2]] <- NULL
       
        result <- combineGoResultsInList(result, method = method)
     }
  
}

keggGeneSet <- function(KEGGid, species = "ath") {
  
  # list genes in given KEGG pathway
  # only works for arabidopsis for now; other species not tested
  # if KEGGid not provided, will return all pathways
  
  # return a named list: with KEGGis as names; and geneid as values
  
  if(missing(KEGGid)) {
    kegg.ath <- keggList("pathway", species)
    names(kegg.ath) <- sub("path:", "", names(kegg.ath))
    KEGGid <- names(kegg.ath)
  }
  
  result <- list()
  for (id in KEGGid) {
      pw <- keggGet(id)
      genes <- pw[[1]]$GENE
      if (is.null(genes)) {
        result[[id]] <- "NA"
        next
      }
      pw2 <- genes[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
      pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
      result[[id]] <- pw2
  } 
  
  return(result)
}

goInKEGG <- function(KEGGid, species = "ath"){
  
  # list associated go terms in given KEGG pathway
  # only works for arabidopsis for now; other species not tested
  # if KEGGid not provided, will return all pathways
  
  if(missing(KEGGid)) {
    kegg.ath <- keggList("pathway", species)
    names(kegg.ath) <- sub("path:", "", names(kegg.ath))
    KEGGid <- names(kegg.ath)
  }
  
  result <- list()
  for (id in KEGGid) {
    pw <- keggGet(id)
    dblinks <- pw[[1]]$DBLINKS
    go <- dblinks[grep("GO:", dblinks)]
    if (is.null(go)) {
      result[[id]] < "NA"
      next
    }
    go <- sub("GO: ", "", go)
    golist <- unlist(str_split(go, pattern = " "))
    result[[id]] <- golist
  } 
  
  return(result)
}

keggTest <- function(KEGGid, geneList, gSet, p.thd = 0.05, method = c("Fisher", "Wilcox")) {
  
  fun <- get(sprintf("keggTest%s", method))
    df <- sapply(KEGGid, function(id) fun(id, geneList, gSet[[id]], p.thd = p.thd)) %>% 
      t() %>% 
      as_tibble() %>% 
      mutate(Set.size = as.numeric(Set.size),
             Significant = as.numeric(Significant),
             p.value = as.numeric(p.value)) %>% 
     rename_with(function(x) str_replace(x, "KEGG_name.*", "KEGG_name"), starts_with("KEGG_name")) 
  
    return(df)
    
}

## test for one pathway with Wilcox method
keggTestWilcox <- function(kegg_id, geneList, gSet, p.thd = 0.05) {
   #pw_genes_set <- keggGeneSet(KEGGid = kegg_id) %>% unlist()
   #cat(sprintf("id:\t%s\n", kegg_id))
   genes_in_pw <- intersect(names(geneList), gSet)
   genes_out_pw <- setdiff(names(geneList), genes_in_pw)
   scores_in <- geneList[genes_in_pw]
   scores_out <- geneList[genes_out_pw]
   if (length(scores_in) > 0){
     p.value <- wilcox.test(scores_in, scores_out, alternative = "less")$p.value
   } else{
     p.value <- NA
   }
   genes_sig <- geneList[geneList < p.thd] %>% names
   sig_in_pw <- intersect(genes_sig, genes_in_pw)
   kegg_name <- keggList("pathway", "ath")[sprintf("path:%s", kegg_id)]
   return(c(KEGG_id = kegg_id, 
            KEGG_name = kegg_name,
            Set.size = length(genes_in_pw), 
            Significant = length(sig_in_pw), 
            p.value = p.value))
  
  
}

keggTestFisher <- function(kegg_id, geneList, gSet, p.thd = 0.05) {
     #pw_genes_set <- keggGeneSet(KEGGid = kegg_id) %>% unlist()
     #cat(sprintf("id:\t%s\n", kegg_id))
     genes_in_pw <- intersect(names(geneList), gSet)
     genes_sig <- geneList[geneList < p.thd] %>% names
     sig_in_pw <- intersect(genes_sig, genes_in_pw)
     
     K <- length(genes_in_pw)
     M <- length(genes_sig)
     X <- length(sig_in_pw)
     N <- length(geneList)
      
      ## Fisher's exact test
      cont.table <- matrix(c(X, K-X, M-X, N-K-M+X),
                           nrow = 2, ncol = 2)
      ftest <- fisher.test(cont.table, alternative = "greater")
      kegg_name <- keggList("pathway", "ath")[sprintf("path:%s", kegg_id)]
      return(c(KEGG_id = kegg_id, 
               KEGG_name = kegg_name,
               Set.size = length(genes_in_pw), 
               Significant = length(sig_in_pw), 
               p.value = ftest$p.value))
      
}

pathwayStats <- function(dds, gset, ref, samp, direction, same.dir = T) {
  
  ## to return a dataframe of pathview statistics
  
  cnts <- counts(dds, normalized = T)
  ref.col <- name2col(dds, ref)
  sample.col <- name2col(dds, samp)
  
  kegg.test <- gage(cnts, gset, ref = ref.col, samp = sample.col, same.dir = same.dir)
  
  if (tolower(direction) == "up") {
    pw.df <- as.data.frame(kegg.test$greater)
    pw.df$direction <- "UP"
  } else {
    pw.df <- as.data.frame(kegg.test$less)
    pw.df$direction <- "DOWN"
  }
  
  pw.df$q.val <- ifelse(is.na(pw.df$q.val), 1, pw.df$q.val)

  pw.df$KEGG.ID <- row.names(pw.df)
  pw.names <- importPathway()
  pw.df <- merge(pw.df, pw.names, by = "KEGG.ID")
  row.names(pw.df) <- NULL
  pw.df <- select(pw.df, KEGG.ID, KEGG.Name, set.size, direction, p.geomean, stat.mean, p.val, q.val)
  pw.df <- pw.df[order(pw.df$q.val),]
  return(pw.df)
  
}

name2col <- function(dds, keys) {
  # to return column number based on key
  
  fcts <- ddsGetFactors(dds)
  
  results <- list()
  for (f in fcts) {
    x <- sapply(keys, function(k) which(dds[[f]] == k))
    results[[f]] <- unlist(x)
  }
  
  result <- Reduce(intersect, results)
  return(result)
}

pathwayGraph <- function(df, pathway.id, pthreshold = 0.05, same.layer = F, 
                         img.path = ".", out.prefix = "", out.suffix = "pathview", 
                         file.type = "pdf") {
  
  # takes a dataframe of gene expression data
  # should have log2FoldChage and padj columns
  # rownames of the dataframe should be gene id (AGI)
  
  pws <- importPathway()
  pw.name <- pws$KEGG.Name[which(pws$KEGG.ID == pathway.id)]
  filename <- sprintf("%s.%s_%s.%s", out.prefix, pathway.id, pw.name, file.type) %>% 
    gsub("/", "-", .) %>% 
    gsub(";", "-", .)
  old.name <- sprintf("%s.%s.%s", pathway.id, out.suffix, file.type) 
  
  lfc <- ifelse(df$padj < pthreshold, df$log2FoldChange, 0)
  
  if(is_tibble(df)) {
    names(lfc) <- df$Gene_id
  } else {
    names(lfc) <- row.names(df)
  }
  
  pathview(lfc, 
           pathway.id = pathway.id, 
           species = "ath", 
           gene.idtype = "tair", 
           same.layer = same.layer,
           out.suffix = out.suffix,
           kegg.native = ifelse(file.type == "pdf", F, T))
  
  # in default, graph saved in the current working directory wit default name
  # to rename and move pathway graph 
  cat(sprintf("rename file to:\n\t%s\n",filename))
  file.rename(old.name, filename)
  cat(sprintf("copy file to:\n\t%s\n", img.path))
  file.copy(filename, file.path(img.path, filename))
  
  # to remove old files
  files2rm <- list.files(pattern = paste0(pathway.id, "*"))
  file.remove(files2rm)

}

pathwayGeneCluster <- function(dds, KEGGid, expr, ntop = 50, direction = "up", imgpath = NA, imgname = NA) {
  ## to plot normalized count heatmap of genes from a given KEGG pathway
  ## dds: a dds object
  ## KEGGid: in form of athxxxxxx
  ## expr: DEseq2 result table; have: geneid, log2foldchange, padj
  
  if (tolower(direction) == "up") {
      gccol <- colorRampPalette(brewer.pal(9, "YlOrRd"))(255)  ## color for genecluster
      expr <- expr[with(expr, padj < 0.05 & log2FoldChange > 0),]

  } else if(tolower(direction) == "down") {
      gccol <- colorRampPalette(brewer.pal(9, "Blues"))(255)  ## color for genecluster
      expr <- expr[with(expr, padj < 0.05 & log2FoldChange < 0),]

  } else {
      gccol <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))  ## color for genecluster
      direction <- "Up/Down"
  }
  

  gene.pw <- genesInKEGG(KEGGid)[[KEGGid]]
  expr <- expr[row.names(expr) %in% gene.pw, ]
  expr <- expr[order(abs(expr$log2FoldChange), decreasing = T),]
  n <- min(ntop, nrow(expr))
  expr <- expr[1:n,]
  expr <- expr[order(expr$log2FoldChange, decreasing = T),]
  
  dds <- ddsCombineFactors(dds)
  dds <- vst(dds)
  dds <- dds[row.names(dds) %in% row.names(expr), ]
  
  row.names(dds) <- ifelse(mcols(dds)$Gene_name == "", row.names(dds), paste0(row.names(dds), "_", mcols(dds)$Gene_name))
  
  img <- list(type = "tiff", w = 10, h = 8, unit = "in", dpi = 300, px = 20, lsz = 1.5, psz = 6)
  pw.name <- keggGet(KEGGid)[[1]]$NAME
  pw.name <- unlist(str_split(pw.name, " - Arabidopsis"))[1]
  title <- sprintf('Top %d %s Regulated Genes in KEGG pathway "%s"', n, direction, pw.name)
  title <- tools::toTitleCase(title)
  mat <- as.matrix(assay(dds))
  colnames(mat) <- colData(dds)$Group
 
  pheatmap(mat, 
           scale="row", 
           trace="none", 
           col = gccol, 
           cellwidth = img$px * 1.6,
           #cellheight = img$px * 0.8,
           border_color = NA,
           angle_col = 90, 
           fontsize = img$px * 0.6,
           show_rownames = T,  
           main = title, 
           filename = ifelse(is.na(imgname), NA, file.path(imgpath, imgname))) 
           #width = img$w, height = img$w)
  
}

goVisBar <- function(data, y = "pvalue", pthreshold = 0.01, maxGO = 20, 
                     col = NA, showID = F, img.path = NA, fn.prefix = NA) {
  value <- as.numeric(data[[y]])
  if(showID) {
    names(value) <- paste0(data$Term, " (", data$GO.ID, ")")
  } else {
    names(value) <- data$Term
  }
  value <- value[order(value)]
  value <- value[value < pthreshold]
  if(length(value)==0) {
    message("NO significantly Enriched GO Term")
    return()
  }
  value <- value[1:min(maxGO, length(value))]
  text <- rev(paste0(value))
  text <- ifelse(text == "1e-30", "< 1e-30", text)
  value <- rev(-log10(value))
  # file.name <- paste("EnrichedGO", fn.prefix, paste0("top(", length(value), ")"), "png", sep = ".")
  fn <- sprintf("EnrichedGO_%s_(top%d).pdf", fn.prefix, length(value))
  paretoh(value, label = text, xlabel = "Enrichment (-logP)", col = col, 
          addline.v = -log10(pthreshold), 
          img.path = img.path, img.name = fn)
  
}

## takes a list of GO-enrichment results and plot p-vallue as a bar plot
goVisBar2 <- function(data, y = "pvalue", pthreshold = 0.01, showThreshold = T, 
                      cols = NA, cat.names = NA, select.method = "allSig", maxGO = 20, 
                      showLegend = T, showID = F, title = "", img.path = NA, img.name = NA) {
  
  if(!class(data)=="list") stop("Need a list of GO Results")
  
  temp <- data
  
  temp <- lapply(temp, function(df) df[order(df$GO.ID),])
  ps <- lapply(temp, function(df) as.numeric(df[[y]]))
  ps <- as.matrix(as.data.frame(ps))
  if(showID) {
    row.names(ps) <- paste0(temp[[1]]$Term, " (", temp[[1]]$GO.ID, ")")
  } else {
    row.names(ps) <- temp[[1]]$Term
  }
  
  anysig <- apply(ps, 1, function(x) any(x < pthreshold))
  allsig <- apply(ps, 1, function(x) all(x < pthreshold))
  if (select.method == "allSig") {
    ps <- ps[allsig, ,drop = F]
  } else {
    ps <- ps[anysig, ,drop = F]
  }
  if (nrow(ps) == 0) {
     message("NO item for plot")
     return()
  } else if (nrow(ps) < 4) {
     showLegend <- F
  }
  
  p.prods <- apply(ps, 1, prod)
  ps <- ps[order(p.prods),, drop = F]
  ps.log <- -log10(ps)
  ps.top <- ps.log[1:min(maxGO, nrow(ps.log)), , drop = F]
  ps.top <- ps.top[c(nrow(ps.top):1), , drop = F]
  

  ## set color
  if(is.na(cols[1])) {
      cols <- brewer.pal(12,"Set3")[1:ncol(ps.top)]
  } else {
      cols <- cols
  }
  
  if(is.na(cat.names)[1]) cat.names <- colnames(ps.top)
  textlength <- max(nchar(row.names(ps.top)))
  
  ## set image file
  if(!is.na(img.path)) {
    # names.split <- sapply(colnames(ps.top), function(v) str_split(v, "\\."))
    # names.df <- t(as.data.frame(names.split))
    # names.ele <- apply(names.df, 2, unique)
    # names.ele2 <- unlist(lapply(names.ele, function(x) paste0(x, collapse = "_")))
    pdf(file.path(img.path, img.name), 
        width = textlength * 0.2, 
        height = nrow(ps.top) * 0.3 + 1.5) 
        #res = 300, unit = "in")
  }
  
  rmar <- max(nchar(colnames(ps.top))) * 0.5
  par(mar = c(6, textlength * 0.45, 1,1))#
  bp <- barplot(t(ps.top), beside = T, horiz = T, 
               # width = 0.4, #space = 0.2, 
                border = NA, axes = F,
                cex.axis = 1.5,
                las = 1, col = cols,
                yaxt = "n")
  if (showLegend) {
        legend("bottomright",
              legend = cat.names,
              fill = cols,
              border = NA, 
              bty = "n",
              xpd = T,
              cex = 1.1,
              inset = c(0,0))
  }
  ax2 <- axis(side = 2, at = bp[seq(2, length(bp), ncol(ps.top))], labels = row.names(ps.top), tick = F, las= 2, xpd = T, col.axis = "gray15", col = "gray15", mgp=c(3,0,0), cex.axis=1.2)
  ax1 <- axis(side = 1, labels= T, col.axis = "gray15", col = "gray15", cex.axis = 1.2, las=1, mgp=c(1,1,0))
  title(xlab = "Enrichment", line = 2.2, cex.lab = 1.5)
  abline(v=0, col = "gray15", lwd = 1.5)
  if(showThreshold) abline(v = -log10(pthreshold), lty = 2, col = "gray90", lwd = 1.5)
  if (!is.na(img.path)) {
     dev.off()
  }
  
}


goVisHM <- function(data, y = "pvalue", keys, include = "all", pthreshold = 0.01, maxGO = 30, showID = F, title = "", img.path = NA, img.name = NA) {
  
  # take a list of go results and produce a heatmap
  # keys: keys to filter results. keys should be substring in the names(data)
  if(!class(data)=="list") stop("Need a list of GO Results")
  
  data.use <- list()
  if (!include == "all") {
    for(n in names(data)) {
      includeExp <- sapply(include, function(x) grepl(x,n))
      if(any(includeExp)) {
        data.use[[n]] <- data[[n]]
      }
    }
  } else {
    data.use <- data
  }
  
  temp <- list()
  for(n in names(data.use)) {
    keyInName <- sapply(keys, function(x) grepl(x, n))
    if (all(keyInName)) {
      temp[[n]] <- data.use[[n]]
    }
  }
  
  whichgo <- unique(str_extract(names(temp), "BP|CC|MF"))
  if(is.na(whichgo)) {
    stop("NO GO type found in names(data)\n")
  } else if(length(whichgo) > 1) {
    stop("Different GO types in data\n")
  } else if(!whichgo %in% c("BP", "MF", "CC")) {
    stop("unknown GO type\n")
  }
  
  temp <- lapply(temp, function(df) df[order(df$GO.ID),])
  ps <- lapply(temp, function(x) as.numeric(x[[y]]))
  ps <- as.matrix(as.data.frame(ps))
  if(showID) {
    row.names(ps) <- paste0(temp[[1]]$Term, " (", temp[[1]]$GO.ID, ")")
  } else {
    row.names(ps) <- temp[[1]]$Term
  }
  
  anysig <- apply(ps, 1, function(x) any(x < pthreshold))
  allsig <- apply(ps, 1, function(x) all(x < pthreshold))
  ps <- ps[anysig,]
  pVars <- rowVars(ps)
  ps <- ps[order(pVars, decreasing = T),]
  ps <- -log10(ps)
  ps <- ps[1:min(nrow(ps), maxGO),]
  
  hmcol <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))
  pheatmap(ps, 
           scale="row", 
           trace="none", 
           col = hmcol, 
           border_color = NA,
           angle_col = 90, 
           fontsize = 16,  
           cellwidth = 30,
           show_rownames = T,  
           filename = ifelse(is.na(img.name), NA, file.path(img.path, img.name)), 
           #width = img$w, height = img$w,
           main=title)
  
  
}

readGO <- function(path, patterns) {
  
  ## read go results to a list
  result <- list()
  
  files <- list.files(path = path, pattern = "*.csv")
  for (f in files) {
    pattern.match <- sapply(patterns, function(x) grepl(x, f))   
    if(!all(pattern.match)) next
    res.name <- paste0(unlist(str_split(f, "\\."))[1:4], collapse = ".")
    res.df <- read.csv(file.path(path,f), header = T, stringsAsFactors = F)
    result[[res.name]] <- res.df
  }
  
  return(result)
}

readPathwayRes <- function(path, patterns) {
  
  ## read go results to a list
  result <- list()
  
  files <- list.files(path = path, pattern = "PathWayEnrichment.*.csv")
  for (f in files) {
    pattern.match <- sapply(patterns, function(x) grepl(x, f))   
    if(!all(pattern.match)) next
    res.name <- paste0(unlist(str_split(f, "\\."))[2:4], collapse = ".")
    res.df <- read.csv(file.path(path,f), header = T, stringsAsFactors = F)
    result[[res.name]] <- res.df
  }
  
  return(result)
}

dfIntersect <- function(data, by.col) {
   
  # input: a list of dataframes with same structure but different number of rows
  # output: return the a list of dataframes based on the intersect column defined by by.col
  if(by.col == 0) {
    data <- lapply(data, function(df) { df[["tempcol"]] <- row.names(df); df})
    by.col <- "tempcol"
  }
  
  x <- lapply(data, function(df) df[[by.col]])
  comm <- Reduce(intersect, x)
  result <- lapply(data, function(df) df[df[[by.col]] %in% comm, ])
  
  result <- lapply(result, function(df) df[,!colnames(df) == "tempcol"])
  return(result)
}

dfUnion <- function(data, by.col, carryover.cols = NA) {
  
  # takes a list of dataframs with same columns
  # return d list of dataframs based on the union of by.col
  # carryover.cols will be filled
  # numeric columns will be filled by NA
  
  comm.dfs <- dfIntersect(data, by.col = by.col)
  if(by.col == 0) {
    data <- lapply(data, function(df) {df[["tempcol"]] <- row.names(df); df})
    comm.dfs <- lapply(comm.dfs, function(df) {df[["tempcol"]] <- row.names(df); df})
    by.col <- "tempcol"
    carryover.cols <- c(carryover.cols, by.col)
  }
  
  ncol <- unique(unlist(lapply(data, ncol)))
  names <- unique(unlist(lapply(data, colnames)))[1:ncol]
  if (length(ncol) > 1) {
    message("different number of columns\nreturned NA")
    return()
  }
  
  uniq.dfs <- lapply(data, function(df) df[!df[[by.col]] %in% comm.dfs[[1]][[by.col]],])
  df <- data.frame(matrix(nrow = 0, ncol = ncol))
  colnames(df) <- names
  for (n in names(uniq.dfs)) { 
     colnames(uniq.dfs[[n]]) <- names
     df <- rbind(df,uniq.dfs[[n]])
  }
  
  numcols <-unlist(sapply(df, function(x) ifelse(any(is.na(x))|suppressMessages(all(!is.na(as.numeric(as.character(x))))), T, F)))
 
  blank <- colnames(df)[!colnames(df) %in% carryover.cols]
  df[, blank] <- ""
  df[, numcols] <- NA
  
  result <- list()
  for (n in names(data)) {
     add <- df[!df[[by.col]] %in% data[[n]][[by.col]], ]
     result[[n]] <- rbind(data[[n]], add)
  }
  
  result <- lapply(result, function(df) df[,!colnames(df) == "tempcol"])
  return(result)
                     
}

go2Gene <- function(go.df, exp.df, go.pthreshold = 0.01) {
  
  # takes a dataframe of go statistics (from goTair), and a dataframe of gene expression
  # return genes that in significant enriched GO terms
  
  go.sig <- filter(go.df, pvalue < go.pthreshold)
  goTerms <- go.df[, c("GO.ID", "Term")]
  genes <- AnnotationDbi::select(org.At.tair.db, keys = go.sig$GO.ID, columns = "TAIR", keytype = "GO")
  genes <- merge(genes, go.df, by.x="GO", by.y="GO.ID", all.x=T)
  genes <- select(genes, TAIR, GO, Term, ONTOLOGY, EVIDENCE)
  
  res <- merge(genes, geneInfo, by.x="TAIR", by.y="Gene_id", all=F)
  
  return(res)
}

importGO <- function(reload = F) {
  
   # to load a dataframe of go terms
   # including 4 columns: GO.ID, Ontology, name of the term, Definition
  
   go.v <- package.version("GO.db")
   goterm.path <- "./genomes/goterms"
   go.filename <- paste("goterms", go.v, "csv", sep = ".")
   go.file <- file.path(goterm.path, go.filename)
   
   if(file.exists(go.file) & !reload) {
     go.table <- read.csv(go.file, header = T)
   } else {
     go.table <- data.frame(Term = Term(GOTERM), Ontology = Ontology(GOTERM), Definition = Definition(GOTERM))
     go.table$GO.ID <- row.names(go.table)
     go.tabel$Term <- gsub(",", ";", go.table$Term)
     go.tabel$Definition <- gsub(",", ";", go.table$Definition)
     
     go.table <- select(go.table, GO.ID, Ontology, Term, Definition)
     row.names(go.table) <- NULL
     write.table(go.table, go.file, row.names = F, col.names = T, sep = ",", quote = F)
   }
   
   go.table$Ontology <- as.factor(go.table$Ontology)
   return(go.table)
   
}

importPathway <- function(species = "ath", reload = F) {
    
    # to load a dataframe of KEGG terms
    # including 4 columns: KEGG.ID, KEGG.NAME, genes, Definition
    
   # go.v <- package.version("GO.db")
    pw.path <- "./genomes/goterms"
    pw.filename <- paste("KEGG", species, "csv", sep = ".")
    pw.file <- file.path(pw.path, pw.filename)
    
    if(file.exists(pw.file) & !reload) {
      pw.table <- read.csv(pw.file, header = T)
    } else {
      pw.ath <- keggList("pathway", species)
      names(pw.ath) <- sub("path:", "", names(pw.ath))
      pw.ath <- sapply(pw.ath, function(k) str_split(k, " - Arabidopsis")[[1]][1])
      pw.ath <- sub(",", ";", pw.ath)
      pw.ath <- as.data.frame(pw.ath)
      pw.ath$ID <- row.names(pw.ath)
      
      pw.genes <- keggGeneSet()
      pw.g <- unlist(lapply(pw.genes, function(g) paste0(g, collapse = "; ")))
      pw.gtable <- as.data.frame(pw.g)
      pw.gtable$ID <- row.names(pw.gtable)
      
      pw.table <- merge(pw.ath, pw.gtable, by="ID")
      
      colnames(pw.table) <- c("KEGG.ID", "KEGG.Name", "KEGG.genes")
      row.names(pw.table) <- NULL
      write.table(pw.table, pw.file, row.names = F, col.names = T, sep = ",", quote = F)
    }
    
    return(pw.table)
}

goKeys <- function() {
  
  key <- list()
  key$hormones <- "hormone|auxin|indoleacetic.*acid|cytokinin|ethylene|Brassinosteroid|abscisic.*acid|gibberellin|salicylic.*acid|jasmonic.*acid|Strigolactone"
  key$development <- "root|shoot|petal|sepal|floral|morphogenesis|development|polarity|pattern formation|cell fate|differentiatioin"
  key$cellfate <- "meristem|stem cell|differentiation|specification|pattern formation|cell fate"
  key$vascular <- "xylem|phloem|cambium|meristem|embryo"
  key$defense <- "bacterium|bacterial|chitin|insect|fungus|fungal|virus|defense|resistance|immune|cell death|autophagy"
  key$abiotic <- "response to|stimulus"
  key$chromsome <- "DNA methylation|DNA recombination|histone|chromosom|recombination|DNA conformation"
  key$transport <- "transport"
  key$deve2 <- paste(key$development, key$cellfate, sep = "|")
  key$expreg <- "translation|transcription|DNA methylation|RNA methylation|histone.*methylation|mRNA splicing|histone acetylation"
  
  return(key)
}
