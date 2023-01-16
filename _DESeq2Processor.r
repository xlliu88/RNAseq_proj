## DESeq2 processing functions


ddsFromFeatureCounts <- function(x, conditions, design, source, gffa3 = NA, 
                                 gene_info = NA, normalize = T) {
  
  # read count table and conditions
  # first column of conditions table should be the same as column names of count matrix
  # value: normalized dds object
  
  if (!tolower(source) %in% c("featurecount", "stringtie", "tximport")) stop("unknown source of count table")
  
  cond <- read.table(conditions, header = T, sep = "\t", stringsAsFactors = F)
  
  if (tolower(source) == "featurecount") {
    fcnts <- read.table(x, header = TRUE, skip = 1, row.names = "Geneid", stringsAsFactors = F)
    fcnts$Chr <- substring(fcnts$Chr,1,1)
    fcnts$Strand <- substring(fcnts$Strand,1,1)
    geneinfo <- fcnts[,1:5]
    counttable <- fcnts[,6:ncol(fcnts)]

    samples <- gsub("([[:digit:]])\\.([[:digit:]])", "\\1-\\2", names(counttable))
    samples <- unlist(sapply(samples, function(s) lapply(str_split(s, "\\."), function(x) x[length(x)-1])))
    names(counttable) <- vgsub(samples, c("\\.", "_ss", "_sorted"), c("-","",""))
    
  } else if(tolower(source) == "stringtie") {
    counttable <- read.table(x, header = T, row.names = "gene_id", sep = ",", stringsAsFactors = F)
    colnames(counttable) <- sub("\\.", "-", colnames(counttable))
    geneinfo <- NA
    
  } else if(tolower(source) == "tximport") {
    counttable <- read.table(x, header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
    counttable <- round(as.matrix(counttable))
    mode(counttable) <- "integer"
    geneinfo <- NA
  }
  
  if (all(colnames(counttable) == cond[,1])) {
    conditions <- cond[, 2:ncol(cond), drop = F]
    rownames(conditions) <- cond[,1]
  } else {
    stop("please check condition settings.")
  }
  
  ## set condition factor levels. levels will be the same as it's appearance order in conditon file
  for (n in 1:ncol(conditions)) {
    lvls <- unique(conditions[,n])
    conditions[,n] <- factor(conditions[,n], levels = lvls)
  }
  counttable <- as.matrix(counttable)
  dds <- DESeqDataSetFromMatrix(counttable, DataFrame(conditions), design = design)
  mcols(dds) <- geneinfo
  
  if (!is.na(gffa3)) {
    geneDes <- extractGeneDesFromGFF3(gffa3)
    dds <- addGeneDescription(dds, geneDes)
  }
   
  if (!is.na(gene_info)) {
    geneDes <- read.table(gene_info, header = T, row.names = 1, sep = "\t")
    dds <- addGeneDescription(dds, geneDes)
  }
  
  if(normalize) dds <- estimateSizeFactors(dds)
  
  return(dds)
}

DEseq2contrast <- function(dds, outpath, ref_levels = NA, countThreshold = 0, alpha = 0.1, overwrite=F) {
  
  if(!file.exists(outpath)) dir.create(outpath)
  
  dds <- dds[which(rowSums(counts(dds)) >= countThreshold), ]
  fct.names <- ddsGetFactors(dds)
  
  if (!all(is.na(ref_levels))) {
      dds <- ddsFactorRelevel(dds, ref_levels = ref_levels)
   }
  
  #dds <- ddsCombineFactors(dds)
  if("Gene_description" %in% colnames(mcols(dds))) {
      geneInfo <- mcols(dds)[, c("Gene_name", "Gene_type", "Gene_description")]
      geneInfo <- geneInfo[order(row.names(geneInfo)),]
  } 
  
  dds2<- DESeq(dds)
  rn <- resultsNames(dds2)
  int.names <- rn[which(containsAll(rn, fct.names))]
  simple.names <- rn[which(containsAny(rn,fct.names) & !containsAll(rn,fct.names))]

  source <- as.character(unique(dds$Source))
  #for (r in resultsNames(dds2)) {
  for(r in c(simple.names, int.names)){
    comp.name <- contrast2Comparison(dds2, r)
    filename <- paste(source, comp.name, "csv", sep = ".")
    res <- results(dds2, name = r, alpha = alpha)
    res <- res[order(row.names(res)), ]
    if (exists("geneInfo")) res <- cbind(res, geneInfo)
   
    if (file.exists(file.path(outpath, filename)) & !overwrite) {
       cat("   !!! output file exist. No data was written.\n")
    } else {
       cat("   writing results... \n")
       cat("       resultsNames: ", r, "\n")
       cat("       output file: ",filename,"\n")
       write.csv(res, file.path(outpath,filename))
       cat("\n")
    }
    
  }
  
  comb.contrast <- expand.grid(simple.names, int.names, stringsAsFactors = F)
  for (r in 1:nrow(comb.contrast)) {
      cont.name <- list(as.character(comb.contrast[r,]))
      comp.name <- contrast2Comparison(dds2, cont.name)
      if(is.na(comp.name)) next
        
      filename <- paste(source, comp.name, "csv", sep = ".")
      
      res <- results(dds2, cont.name)
      res <- res[order(row.names(res)), ]
      if (exists("geneInfo")) res <- cbind(res, geneInfo)
      
      if (file.exists(file.path(outpath, filename)) & !overwrite) {
        cat("   !!! output file exist. No data was written.\n")
      } else {
        cat("   writing results... \n")
        cat("       combined resultsNames: ", paste0(unlist(cont.name), collapse = ", "), "\n")
        cat("       output file: ",filename,"\n")
        write.csv(res, file.path(outpath,filename))
        cat("\n")
      }
      
  }
  
  
}

ddsFactorRelevel <- function(dds, ref_levels) {
  
  fcts <- ddsGetFactors(dds)
  for(f in fcts) {
    ref <- ref_levels[which(ref_levels %in% levels(dds[[f]]))]
    
    if (length(ref)==0) {
       ref <- levels(dds[[f]])[1]
    } else if (length(ref) > 1) {
       msg <- sprintf("more then one ref level found for Factor %s", f) 
       stop(msg)
    }  
    
    dds[[f]] <- relevel(dds[[f]], ref = ref)
    }
  
  return(dds)
}

counts2 <- function(dds, format = "long", normalize = T) {
  # to return normalized counts of a dds
  # only support melted format for now
  cnts <- counts(dds, normalized = normalize)
  cond <- as.data.frame(colData(dds))
  fcts <- ddsGetFactors(dds)
  
  if (!format == "long") return(cnts)
  
  cnts.melt <- melt(cnts, value.name = "Counts")
  colnames(cnts.melt) <- c("Gene_id", "Label", "Counts")
  for (f in fcts) {
    cnts.melt[[f]] <- cond[match(cnts.melt[,2], row.names(cond)), f]
  }
  
  cnts.melt <- cnts.melt[, c("Gene_id", "Label", fcts, "Counts")]
  return(cnts.melt)
}

DEseq2pairwise <- function(dds, outpath, q = F, overwrite=F) {
  
  if(!file.exists(outpath)) dir.create(outpath) 
  ## combine all factors into one
  dds <- ddsCombineFactors(dds, combined_name = "Group")
  
  ## calculate and write pairwise comparison results.
  source <- as.character(unique(dds$Source))
  idx <- t(combn(levels(dds$Group),2))
  for (i in 1:nrow(idx)) {
    subdds <- ddsSlice(dds, factorName="Group", factorLevel = idx[i,]) 
    subdds$Group <- relevel(subdds$Group, ref = idx[i,1])
    design(subdds) <- ~ Group
    
    subdds <- DESeq(subdds, quiet=q)
    res <- results(subdds)
    
    res <- res[order(row.names(res)), ]
    gd <- as.data.frame(mcols(dds)[, c("Gene_name", "Gene_type", "Gene_description")])
    gd <- gd[order(row.names(gd)),]
    res <- cbind(res, gd)
    #res[["Gene_name"]] <- gd$Gene_name
    #res[["Gene_type"]] <- gd$Gene_type
    #res[["Gene_description"]] <- gd$Gene_description
    
    comparison <- paste0(rev(idx[i,]), collapse = ".vs.")
    output_name <- paste0(source, ".", comparison, ".csv")
    
    cat("writing results for ", comparison, "\n")
    cat("   file: ", output_name, "\n")
    if (file.exists(file.path(outpath, output_name)) & !overwrite) {
       cat("   !!! output file exist. No data was written.\n")
    } else {
      cat("   writing results... \n")
      write.csv(res, file.path(outpath,output_name))
      cat("\n")
    }
  }

}

extractGeneDesFromGFF3 <- function(gff3, overwrite=F) {
  
  geneDes_file <- sub("\\.gff3",".annotation", gff3)
  if(file.exists(geneDes_file) & overwrite==F) {
    geneDes <- read.table(geneDes_file, header = T, row.names = 1, sep = "\t", quote="", stringsAsFactors = F)
    return(geneDes)
  }
  
  gff <- read.table(gff3, header = F, sep = "\t", quote = "", fill = T, skip = 15)
  names(gff) <- c("Chr", "Source","Feature", "Start", "End", "Score", "Strand","Frame", "Attribute")
  gff <- filter(gff, Feature %in% c("gene","ncRNA_gene"))
  
  gff$Gene_id <- str_extract(gff$Attribute, "ID=gene:.*?;")
  gff$Gene_id <- vgsub(gff$Gene_id, c("ID=gene:", ";"), rep("",2))
  gff$Gene_type <- str_extract(gff$Attribute, "biotype=.*?;")
  gff$Gene_type <- vgsub(gff$Gene_type, c("biotype=", ";"), rep("",2))
  gff$Gene_name <- str_extract(gff$Attribute, "Name=.*?;")
  gff$Gene_name <- vgsub(gff$Gene_name, c("Name=", ";"), rep("",2))
  gff$Gene_name[which(is.na(gff$Gene_name))] <- ""
  
  
  gff$Gene_description <- vgsub(gff$Attribute, c("ID=gene:.*?;", "Name=", "biotype=", "description=", "\\[Source:.*araport11", "gene_id=.*?;", "logic_name=araport11"), rep("", 7))
  gff$Gene_description <- vgsub(gff$Gene_description, c("\\%3B", "\\%2C"), rep(";", 2))
  row.names(gff) <- gff$Gene_id
  geneDes <- select(gff, Gene_id, Gene_name, Gene_type, Gene_description)
  
  
  write.table(geneDes, geneDes_file, row.names = F, col.names = T, quote = F, sep = "\t")
  return(geneDes)
}

addGeneDescription <- function(dds, gene_description){
  #browser()
  #dds <- dds[order(row.names(dds)),]
  #gene_description <- gene_description[order(row.names(gene_description)), ]
  match.id <- match(row.names(dds), row.names(gene_description))
  
  if(any(is.na(match.id))) {
    stop("Gene_description file not complete!")
  }
  
  df <- cbind(mcols(dds), gene_description[match.id, ])
  na_cols <- apply(df, 2, function(x) all(is.na(x)))
  mcols(dds) <- df[, !na_cols]
  # if (all(row.names(mcols(dds)) == row.names(gene_description))) {
  #   mcols(dds) <- cbind(mcols(dds), gene_description)
  #   return(dds)
  #   
  # } else if (all(row.names(dds) %in% row.names(gene_description))){
  #   message("dds object has less genes then gene_description file\n")
  #   mcols(dds) <- simpleMerge(mcols(dds), gene_description, by=0, all.x=T)
  #   
  # } else if (all(row.names(gene_description) %in% row.names(dds))) {
  #   message("dds object has more genes then gene_description file\n Some genes will not has description information")
  #   mcols(dds) <- simpleMerge(mcols(dds), gene_description, by=0, all.x=T)
  #   
  # } else {
  #   stop("dds and gene_description have different sets of row names")
  # }
  return(dds) 
}

ddsAddGeneInfo <- function(dds, stgtf) {
  
  ## extract gene info from gtf assemblied by StringTie
  
  df.gtf <- read.table(stgtf, skip=2, header=F, sep = "\t",stringsAsFactors = F)
  names(df.gtf) <- c("Chr", "Source","Feature", "Start", "End", "Score", "Strand","Frame", "Attribute")
  df.gtf <- filter(df.gtf, Feature == "transcript")
  
  geneid <- str_extract(df.gtf$Attribute, "gene_id\\s\\w*\\.\\d*\\;")
  geneid <- vgsub(geneid, c("gene_id\\s", ";"), rep("",2))
  
  refgeneid <- str_extract(df.gtf$Attribute, "ref_gene_id\\s\\w*\\;")
  refgeneid <- vgsub(refgeneid, c("ref_gene_id\\s", ";"), rep("",2))
  df.gtf$stID <- geneid
  df.gtf$AGI <- refgeneid
  df.gtf$stID[is.na(df.gtf$stID)] <- df.gtf$AGI[is.na(df.gtf$stID)]
  
  res <- df.gtf[!duplicated(df.gtf$stID), ]
  noAGI <- which(is.na(res$AGI))
  res$AGI[noAGI] <- ""
  row.names(res) <- res$stID
  res <- select(res, Chr, Source, Start, End, Score, Strand, AGI)
  
  res <- res[order(row.names(res)), ]
  dds <- dds[order(row.names(dds)), ]
  if(all(row.names(dds)==row.names(res))) {
      #row.names(dds) <- row.names(res)
      mcols(dds) <- res
    } else {
      stop("row.names of dds don't match with new ids")
    }
  
  return(dds)
  
}

ddsCombineFactors <- function(dds, combined_name = "Group", collapse = "_") {
  
  ## to combine all factors in dds design into group
  
  fts <- ddsGetFactors(dds)
  cols <- match(fts, colnames(colData(dds)))
  condtable <- as.data.frame(colData(dds)[, cols]) ## as.data.frame: convert S4 DataFrame to S3 "data.frame"
  if(length(fts) == 1) {
     dds[[combined_name]] <- dds[[fts]]
     return(dds)
  } 
  fctlvls <- list()
  for(f in fts){
    fctlvls[[f]] <- levels(condtable[[f]])
  }
  tmp <- expand.grid(fctlvls)
  comblvl <- apply(tmp, 1, paste, collapse = collapse)
  
  dds[[combined_name]] <- apply(condtable, 1, paste, collapse = collapse)
  dds[[combined_name]] <- factor(dds[[combined_name]], levels = comblvl)
  return(dds)
  
}

ddsGetFactors <- function(dds) {
  
  # to get factors in dds design, return factors as a vector
  
  d <- as.character(dds@design)[2]
  fts <- unlist(strsplit(d, "\\+"))
  fts <- fts[which(!grepl(":", fts))]
  fts <- trimws(fts)
  
  return(fts)
}

ddsSlice <- function(dds, dim=2, factorName, factorLevel) {
  
  if (dim %in% c(1, "row")) {
    cat(sprintf("slice by row. Filtering with factor: %s\n", factorName))
    message("under construction")
  } else if ( dim %in% c(2, "col", "column")) {
    cat(sprintf("slice by column. Filtering with factor: %s\n", factorName))
    dds <- dds[, dds[[factorName]] %in% factorLevel]
    
    for (c in names(colData(dds))) {
      
      if ( is.factor(dds[[c]])) {
        dds[[c]] <- droplevels(dds[[c]])
      }
    }
  }
  
  return(dds)
}

ddsRemoveSample <- function(dds, samples) {
  c <- colData(dds)
  if (!all(samples %in% row.names(c))) stop("not all samples to remove found in assay")
  
  keep <- which(!row.names(c) %in% samples)
  dds <- dds[, keep]
  
  return(dds)
  
}

resultDiff <- function(path1, results1, path2,results2, result_names, cat.names = c("Result1", "Result2"), lfc_cutoff = 0, padj_cutoff = 0.05, genetype = "all", outpath, overwrite = F) {
  
  ##  to compare DEseq2 results from two methods.
  ##  the comparison based on padj
  ##  outputs
  ##    a file summarize number of overlapped de genes, overlapped de.up, overlapped de.down, result1.up, result1.down, result2.up, result2.down
  ##    venn diagram of overlapped de genes,  and genes only found in one method
  ##    venn diagram of overlapped up/down genes, and up/down genes only found in one method
  ##    a csv file with list of overlapped genes, overlapped up, overlapped down; 
  ##    
  ##  returns a data frame with 13 cols
  ##    geneID, baseMean.ref, log2FC.ref, padj.ref, 
  ##            baseMean.trt, log2FC.trt, padj.trt, 
  ##            trt, comparison, difftype, gene_name, gene_type, gene_description
  
  
  if (!all(names(results1)==names(results2))) stop("results list have different names")
  if (!all(names(results1)==names(result_names))) stop("different names in results and result_names")
  
  cutoffs <- paste0("[", "padj-", padj_cutoff, " ", "lfc-", lfc_cutoff, "]")
  venndir <- paste("venn", cutoffs, genetype, sep=".")
  vennpath <- file.path(outpath, venndir)
  if (!file.exists(outpath)) dir.create(outpath)
  if (!file.exists(vennpath)) dir.create(vennpath)
  
  df.diff <- data.frame()
  de.summary <- data.frame()
  rowname_summ <- c()
  for (n in names(results1)) {
    result_name <- result_names[[n]]
    cols <- brewer.pal(8, "Set3")[c(1,4)]
    #label <- paste("de", result_name)
    
    res1 <- read.csv(file.path(path1, files[[n]]), header = T, stringsAsFactors = F)
    res2 <- read.csv(file.path(path2, files[[n]]), header = T, stringsAsFactors = F)
    row.names(res1) <- res1$X
    res1 <- res1[order(row.names(res1)),]
    row.names(res2) <- res2$X
    res2 <- res2[order(row.names(res2)),]
    ## differentially expressed genes
    res1.de <- filter(res1, !is.na(padj) & padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
    res2.de <- filter(res2, !is.na(padj) & padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
    
    ## protein genes
    if (!genetype=="all") {
      res1.de <- filter(res1.de, Gene_type %in% genetype)
      res2.de <- filter(res2.de, Gene_type %in% genetype)
    }
    
    ## upregulated and downregulated
    res1.up <- filter(res1.de, log2FoldChange > 0)
    res1.down <- filter(res1.de, log2FoldChange < 0)
    res2.up <- filter(res2.de, log2FoldChange > 0)
    res2.down <- filter(res2.de, log2FoldChange < 0)
    
    ## venn plot
    vennlist <- list(res1.de$X, res2.de$X)
    vennlist.up <- list(res1.up$X, res2.up$X)
    vennlist.down <- list(res1.down$X, res2.down$X)
    names(vennlist) <- names(vennlist.up) <- names(vennlist.down) <- cat.names
    
    o <- calculate.overlap(vennlist)
    o.up <- calculate.overlap(vennlist.up)
    o.down <- calculate.overlap(vennlist.down)
    if (all(o$a1 %in% o$a3) & all(o$a2 %in% o$a3)) next; # if not difference, skip to next comparison
    
    venn.name <- file.path(vennpath, paste0("de.", result_name, ".tiff"))
    venn.name.up <- file.path(vennpath, paste0("up.", result_name, ".tiff"))
    venn.name.down <- file.path(vennpath, paste0("down.", result_name, ".tiff"))
    
    venn.diagram(vennlist, filename = venn.name, 
                 fill = cols, col = cols, 
                 main = paste(result_name,"(DEG)"), main.cex = 1.5,
                 cat.cex = 1.5, cex = 1)
    
    venn.diagram(vennlist.up, filename = venn.name.up, 
                 fill = cols, col=cols, 
                 main = paste(result_name,"(DEG-up)"), main.cex = 1.5,
                 cat.cex = 1.5, cex = 1)
    
    venn.diagram(vennlist.down, filename = venn.name.down, 
                 fill = cols, col=cols, 
                 main=paste(result_name,"(DEG-down)"), main.cex = 1.5,
                 cat.cex = 1.5, cex = 1)
    
    colnames(res1) <- paste(colnames(res1), "ref", sep=".")
    colnames(res2) <- paste(colnames(res2), "trt", sep=".")
    
    if(all(row.names(res1)==row.names(res2))) {
      res <- cbind(res1, res2)
    } else {
      res <- simpleMerge(res1, res2, by.x = "X.ref", by.y = "X.trt", all=T)
    }
    
    res$Diff_type <- NA
    res$Diff_type[which(row.names(res1) %in% o$a3)] <- "Overlapped"
    res$Diff_type[which(row.names(res1) %in% o$a1 & !row.names(res1) %in% o$a3)] <- "Missed"
    res$Diff_type[which(row.names(res2) %in% o$a2 & !row.names(res2) %in% o$a3)] <- "New"
    res.diff <- filter(res, !is.na(Diff_type))
    
    res.diff$Trt <- cat.names[2]
    res.diff$Comparison <- result_name
    
    res.diff <- select(res.diff, X.ref, baseMean.ref, baseMean.trt, log2FoldChange.ref, log2FoldChange.trt, padj.ref, padj.trt, 
                       Trt, Comparison, Diff_type, 
                       Gene_name.ref, Gene_type.ref, Gene_description.ref) 
    names(res.diff) <- c("GeneID","baseMean.ref", "baseMean.trt", "LFC.ref", "LFC.trt", "padj.ref", "padj.trt", 
                         "Trt", "Comparison", "Diff_type", 
                         "Gene_name", "Gene_type", "Gene_description")
    
    df.diff <- rbind(df.diff, res.diff)
    
    summ <- c(length(o$a1), length(o.up$a2), length(o$a3), 
              length(o.up$a1), length(o.up$a2), length(o.up$a3), 
              length(o.down$a1), length(o.down$a2), length(o.down$a3))
    de.summary <- rbind(de.summary, summ)
    rowname_summ <- c(rowname_summ, result_name)
    cat("comparison: ", result_name, "\n")
    cat("    # of total genes in ", cat.names[1], ": ", nrow(res1), "\n")
    cat("      of which diff expressed: ", nrow(res1.de), "\n")
    cat("    # of total genes in ", cat.names[2], " result: ", nrow(res2), "\n")
    cat("      of which diff expressed: ", nrow(res2.de), "\n")
    cat("    # overlap of XL & Bill de: ", length(o$a3), "\n")
    cat("============\n")
    
  }
  
  names(de.summary) <- c(paste("de", c(cat.names,"overlap"), sep = "."),
                         paste("up", c(cat.names,"overlap"), sep = "."),
                         paste("down", c(cat.names,"overlap"), sep = "."))
  row.names(de.summary) <- rowname_summ
  
  base_name <- paste(paste0(cat.names, collapse = ".vs."), cutoffs, genetype, sep = "_")
  sumfile.name <- paste(base_name, "compare.summary.csv", sep = "_")
  df.diff.name <- paste(base_name, "different-DEGs.csv", sep = "_")
  write.csv(de.summary, file.path(outpath, sumfile.name), quote=F, row.names = T)
  write.table(df.diff, file.path(outpath, df.diff.name), quote=F, row.names = F, col.names = T, sep = ",")
  
  do.call(file.remove, list(list.files(vennpath, pattern = "*.log", full.names = T)))
  
  ## prepare df.diff
  df.diff$Trt <- as.factor(df.diff$Trt)
  df.diff$Comparison <- as.factor(df.diff$Comparison)
  df.diff$Diff_type <- as.factor(df.diff$Diff_type)
  return(df.diff)
}

getOverlapGenes <- function(result.files, outpath, result.names = NA, padj_cutoff = 0.05, lfc_cutoff = 1, dds = NA, plotcounts = NA) {
  
  # a function to find overlap genes in two or more comparisons
  # output:
  #   venn diagram
  #   list of genes overlapped
  #
  
  if (length(result.files) < 2) stop("more then 2 results needed for comparison")
  if (length(result.names) < 2) {
    result.files <- paste0("res", 1:length(result.files))
  } else if (!length(result.files) == length(result.names)) {
    n <- length(result.files)
    m <- length(result.names)
    if (n > m) result.names <- c(result.names, paste0("res", (m+1):n))
    if (n < m) result.names <- result.names[1:n]
  }
  
  ## basic setup
  if(!file.exists(file.path(outpath,"plots"))) dir.create(file.path(outpath, "plots"))
  cutoffs <- paste0("[padj-",padj_cutoff," lfc-",lfc_cutoff, "]")
  label <- paste0(result.names, collapse = ".vs.")
  
  ## get ovlapped genes
  results <- list()
  res.df <- data.frame()
  for (f in result.files) {
    
    idx <- which(result.files == f)
    name <- result.names[idx]
    
    res <- read.csv(file.path(outpath, f), row.names=1, stringsAsFactors = F, header = T)
    res$padj <- ifelse(is.na(res$padj), 1, res$padj)
    res$GeneID <- row.names(res)
    res.sign <- filter(res, padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
    res.sign <- res.sign[order(res.sign$padj), ]
    results[[name]] <- res.sign$GeneID
  }
  
  ## set up for venn plot
  venn.cols <- brewer.pal(8, "Set2")[1:length(result.files)]
  venn.name <- paste0(label, "_", cutoffs, "_venn", ".tiff")
  o <- calculate.overlap(results)
  venn.diagram(results, file.path(outpath, "plots", venn.name), 
               fill = venn.cols, col = venn.cols, 
               category.names = names(results),
               main = paste0("Overlapped DEGs - ", label))
  do.call(file.remove, list(list.files(file.path(outpath, "plots"), pattern = "*.log", full.names = T)))
  
  ## plotCounts
  if (is.na(plotcounts)) {
    return(results)
  } else if (!class(dds) == "DESeqDataSet") {
    msg <- sprintf(dds, "is not a DESeqDataSet Object")
    stop(msg)
  }
  
  plt.title <- paste0("Overlapped-DEGs: ", label)
  plt.label <- paste(label, cutoffs, "Overlabpped-DEGs", sep = "_")
  genes <- o[[length(o)]]
  
  if (length(genes) > 210) {
    cat("gene list too long, first 200 will be ploted\n")
    genes <- genes[1:210]
  }
  
  ng <- length(genes)
  step <- 30
  for (i in 1:ceiling(ng/step)) {
    start <- step*(i-1) + 1
    end <- step * i
    cat("plot genes ", start, " to ", min(end, ng))
    cat("\n")
    imagename <- paste0(plt.label, "_", i, ".tiff")
    
    subgenes <- genes[start:min(end, ng)]
    ncol = min(length(subgenes), 6)
    plotCounts2(dds.LCM, subgenes, ncol = ncol, main = plt.title)
    ggsave(file.path(outpath, "plots", imagename), device = "tiff", scale = 1, width = 12, height = 8, units = "in", dpi = 150)   
  }
  
  return(o)
  
}

simpleMerge <- function(df1, df2, by=NULL, by.x=NULL, by.y=NULL, all.x=F, all.y=F, all=F) {
  
  if (is.null(c(by, by.x))) stop("Merge Column for df1 is missing")
  if (is.null(c(by, by.y))) stop("Merge Column for df2 is missing")
  if (is.null(by.x)) by.x <- by
  if (is.null(by.y)) by.y <- by
  
  if (is.character(by.x) & !by.x %in% colnames(df1)) stop("Merge column for df1 is not found")
  if (is.character(by.y) & !by.y %in% colnames(df2)) stop("Merge column for df2 is not found")
  
  if (is.null(by)) {
    row.names(df1) <- df1[[by.x]]
    row.names(df2) <- df2[[by.y]]
    by=0
  } 
  if (ncol(df1) == 1) df1$Row.names <- row.names(df1)
  if (ncol(df2) == 1) df2$Row.names <- row.names(df2)
  
  if (by==0) {
      df1 <- df1[order(row.names(df1)),]
      df2 <- df2[order(row.names(df2)),]
      
      df1.overlap <- df1[which(row.names(df1) %in% row.names(df2)),]
      df2.overlap <- df2[which(row.names(df2) %in% row.names(df1)),]
      df1.left <- df1[which(!row.names(df1) %in% row.names(df2)),]
      df2.left <- df2[which(!row.names(df2) %in% row.names(df1)),]
      
      df.overlap <- cbind(df1.overlap[order(row.names(df1.overlap)),], 
                          df2.overlap[order(row.names(df2.overlap)),])
      if("Row.names" %in% colnames(df.overlap)) df.overlap <- df.overlap[, !colnames(df.overlap)=="Row.names"]
      
      if(all.x==T) {
        df2.left <- data.frame(matrix("", nrow = nrow(df1.left),ncol=ncol(df2.left)))
        names(df2.left) <- names(df2)
        df.left <- cbind(df1.left, df2.left)
        
      } else if (all.y==T) {
        df1.left <- data.frame(matrix("", nrow = nrow(df2.left),ncol=ncol(df1.left)))
        names(df1.left) <- names(df1)
        df.left <- cbind(df1.left, df2.left)
        
      } else if ((all.x==T & all.y==T) | all==T) {
        df.left <- merge(df1.left,df2.left,by = 0, all=all)
        row.names(df.left) <- df.left$Row.names
        keep <- colnames(df.left)[!grepl("Row.names", colnames(df.left))]
        df.left <- select(df.left, keep)
      } else {
        df.left <- data.frame(matrix(nrow = 0, ncol = ncol(df.overlap)))
      }
   }
  
  if("Row.names" %in% colnames(df.left)) df.left <- df.left[, !colnames(df.left)=="Row.names"]
  colnames(df.left) <- colnames(df.overlap)
  merged <- rbind(df.overlap, df.left)
  merged <- merged[, colSums(is.na(merged)) < nrow(merged)]
  
  return(merged) 
}

vgsub <- function(x, v_pattern, v_replacement, fixed = TRUE) {

  if(length(v_replacement) == 1) v_replacement <- rep(v_replacement, times = length(v_pattern))
  if (!length(v_pattern) == length(v_replacement)) {
    stop("Vector of pattern and Vector of replacement have different length")
  }
  
  for (i in 1: length(v_pattern)){
    x <- gsub(v_pattern[i], v_replacement[i], x, fixed = F)
  }
  return(x)
}

containsAll <- function(v1, v2) {
   ## to return elements in V1 that contains all elements in V2 as substring
  
  e <- rep(T, length(v1))
  for (v in v2) {
    ev <- grepl(v, v1)
    e <- e & ev
  }
  
  return(e)
}

containsAny <- function(v1, v2) {
  ## to return elements in V1 that contains all elements in V2 as substring
  
  e <- grepl(paste(v2, collapse = "|"), v1)

  return(e)
}

contrast2Comparison <- function(dds, resName) {
  
  # to convert resultsName / contrast list to comparison
  # dds should be a dds after running DESeq() 
  # factors should be releveled
  # for now only works for upto 2 factor design.
  
  fct.names <- ddsGetFactors(dds)
  fcts <- list()
  for (f in fct.names) {
    fcts[[f]] <- levels(dds[[f]])
  }
  
  # "Condition_BCN_vs_Ctrl": WT.BCN.vs.WT.Ctrl
  # "Genotype_II_vs_I": WT.ctrl.vs.clv.ctrl
  if(!class(resName)=="list" & length(resName)==1) {
      if(containsAll(resName, fct.names) & length(fct.names > 1)) {
         ## Interaction comparison
         temp <- vgsub(resName, fct.names, rep("", length(fct.names)))
         comp.name <- paste("interaction", temp, sep = "_")
      } else {
         ## simple comparison
         comp.fct <- fct.names[which(sapply(fct.names, function(x) grepl(x, resName)))]
         resName.split <- unlist(str_split(resName, "_"))
         comp.lvls <- resName.split[which(resName.split %in% fcts[[comp.fct]])]
        
         left.fct <- fct.names[which(!fct.names == comp.fct)]
        
         ## to compile comparison
         comp.name <- ""
         for (n in names(fcts)) {
            if(n == comp.fct) {
               comp.name <- paste(comp.name, comp.lvls, sep = "_")
            } else {
               comp.name <- paste(comp.name, fcts[[n]][1], sep = "_")
            }
         }
         comp.name <- paste0(substring(comp.name, 2), collapse = ".vs.")
      }
    
  } else if (class(resName)=="list") {
      #browser()
      #ctrstNames <- unlist(resName)
      main.eff <- as.character(resName[[1]][1])
      int.eff <- as.character(resName[[1]][2])
      
      comp.fct <- fct.names[which(sapply(fct.names, function(x) grepl(x, main.eff)))]
      main.split <- unlist(str_split(main.eff, "_"))
      comp.lvls <- main.split[which(main.split %in% fcts[[comp.fct]])]
      
      int.lvls <- vgsub(int.eff, fct.names, rep("", length(fct.names)))
      int.lvls <- unlist(str_split(int.lvls, "\\."))
      
      if(!(any(int.lvls %in% comp.lvls))) {
         msg <- sprintf("invalid combination of resultsNames:\n\t%s \n\treturn: NA\n", 
                        paste0(unlist(resName), collapse = ", "))
         warning(msg)
         return(NA)
      }
      
      int.lvl <- int.lvls[which(!int.lvls %in% comp.lvls)]
      
      comp.name <- ""
      for(n in names(fcts)) {
        if(n == comp.fct) {
          comp.name <- paste(comp.name, comp.lvls, sep = "_")
          cat(comp.name, "\n")
        } else if (int.lvl %in% fcts[[n]]) {
          comp.name <- paste(comp.name, int.lvl, sep = "_")
          cat(comp.name, "\n")
        }
      }
      comp.name <- paste0(substring(comp.name, 2), collapse = ".vs.")
  } else {
    ## 
    comp.name <- NA
  }
  
  return(comp.name)
  
}

# pattern to filter file name
# filter to filter with metadata
# for now, it only merge single page PDFs
mergePDFs <- function(path, pattern = NA, filter = list(), 
                     order_by = c("name", "ct", "mt"), keep = T, out.name = NA) {
  
  if(!class(filter) == "list") stop("filter has to be a list")
  
  pdfs <- list.files(path = path, pattern = "*.pdf$")
  pages <- unname(sapply(file.path(path, pdfs), pdf_length))
  pdfs <- pdfs[which(pages == 1)]
  pdfs <- file.path(path, pdfs)
  if(is.na(out.name)) out.name <- "combined.pdf"
  combined <- pdf_combine(pdfs, output = file.path(path, out.name))
  
  if(!keep) {
    file.remove(file.path(path, pdfs))
  }
}
  