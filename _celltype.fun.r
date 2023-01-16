## To Calculated cell specific genes;
## cty: a given cell type
## LM: a linear model
## queries: a matrix of 0s and 1s to define comparison pairs; 1: do compare, 0: don't do compare
## value: return a list, where
##        each slot contains genes expressing higher in the given cell-type than the compared one
##        and the slot name is the comparison . 

CellSpecificGenes <- function(cty, LM, queries, svobj = NULL, 
                              addFC.thd = NA, addp.thd = NA) {
  if(!cty %in% row.names(queries)) stop("cell type doesn't exist")
  if(!is.na(addFC.thd)) addp.thd <- ifelse(is.na(addp.thd), 1, addp.thd)
  if(!is.na(addp.thd)) addFC.thd <- ifelse(is.na(addFC.thd), 1, addFC.thd)

  ## build a contrast matrix
  # in the model ~ marker, since the first factor is set as the base level,
  # everything else is relative to this level
  # so when compare this factor to any other factors, you don't have to substract
  # this level again. thus the first row of contrast matrix is set to 0 
  qlist <- paste0(cty, "-", colnames(queries)[queries[cty,]==1]) 
  contrast.matrix <- makeContrasts(contrasts = qlist, levels = LM$design)
  contrast.matrix[1, ] <- 0
  # 
  # if(!is.null(svobj)) {
  #   if (!class(svobj) == "list") stop("wrong sva object")
  #   if (is.null(svobj$n.sv)) stop("wrong sva object items")
  #   
  #   contrast.matrix <- rbind(contrast.matrix, 
  #                            matrix(rep(0, svobj$n.sv * length(qlist)), 
  #                                   ncol = length(qlist)))
  # }
  
  ## build a linear model
  fit2 <- contrasts.fit(LM, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  ## get list of genes that has higher expression in given cell type 
  ## compare to other cell types
  pw.result <- list()
  for (q in qlist) {
    res <- topTable(fit2, coef = q, number = nrow(fit2))
    g1 <- row.names(res)[res$logFC >= log2(fc.thd) & res$adj.P.Val <= p.thd]
    g2 <- NULL
    if(!is.na(addFC.thd)) {
      idx <- which(res$logFC >= log2(addFC.thd) & res$adj.P.Val <= addp.thd)
      g2 <- row.names(res)[idx]
    }
    pw.result[[q]] <- unique(c(g1,g2))
  }
  
  return(pw.result)
}

getCellSpecificGenes2 <- function(eset, design, queries, p2g, sva = F, svobj = NULL,
                                 out.file = NA, addFC.thd = NA, addp.thd = NA) {
    
    pdata <- phenoData(eset)@data
    mks <- factor(pdata$Marker)
    queries <- queries[row.names(queries) %in% mks, colnames(queries) %in% mks]
    if(!("Marker" %in% colnames(pdata))) {
      stop("Marker column not found in targets!")
    }
    if(!all(c("array_element_name", "locus") %in% colnames(p2g))) {
      stop("check probe to gene ID conversion table")
    }
    
    if(sva & is.null(svobj)) stop("sva object missing")
    if(!sva) svobj <- NULL
    # if(sva) {
    #   design <- model.matrix( ~ mks)
    #   colnames(design) <- c("Intercept", levels(mks))
    #   n.sv <- num.sv(exprs(eset), design, method = "leek")
    #   svobj <- sva(exprs(eset), design, n.sv=n.sv)
    #   design <- cbind(design, svobj$sv)
    # } else {
    #   design <- model.matrix( ~ mks)
    #   colnames(design) <- c("Intercept", levels(mks))
    #   svobj <- NULL
    # }

    fit <- lmFit(eset, design)
    
    cat("  calculate enriched genes in each cell type...\n\t")
    cty.probes <- list()
    for (cty in levels(mks)) {
      cat(sprintf("%s, ", cty))
      temp <- CellSpecificGenes(cty, fit, queries = queries, svobj = svobj, 
                                addFC.thd = addFC.thd, addp.thd = addp.thd)
      cty.probes[[cty]] <- calculate.overlap2(temp)
    }
    cat("\n")
    
    ## convert probeID to geneID
    cat(" Convert probeID to geneID...\n")
    cty.genes <- lapply(cty.probes, function(x) {
      p2g[match(x, p2g$array_element_name), ]$locus})
    cty.genes <- lapply(cty.genes, function(x) {
      unlist(sapply(x, function(g) unlist(str_split(g, ";"))))})
    cty.genes <- lapply(cty.genes, function(x) x[!x == "no_match"])
    cty.genes <- lapply(cty.genes, unique)
    emptymarker <- which(sapply(cty.genes, is.null))
    if(any(emptymarker)) cty.genes[[emptymarker]] <- " "
    
    ## convert result to a dataframe and save to a file
    cat(" Save cell type specific genes to file...\n")
    if(is.na(out.file)) {
       enrdfile <- sprintf("EnrdGeneList.%s.fc%1.1f.p%s.csv", 
                           mthd, fc.thd, 
                           formatC(p.thd, format = "e", digits = 0))
       out.file <- file.path(out.path, enrdfile)
    }
    n.obs <- sapply(cty.genes, length)
    seq.max <- seq_len(max(n.obs))
    mat <- sapply(cty.genes, "[", i = seq.max)
    mat[is.na(mat)] <- ""
    df.cty.genes <- as.data.frame(mat)
    cat("   Write result in: ", out.file, "\n")
    write.table(df.cty.genes, out.file, sep = ",", 
                row.names = F, col.names = T, quote = F)
    
    return(cty.genes)
}

getCellSpecificGenes <- function(eset, targets, queries, p2g, sva = F, 
                                 out.file = NA, addFC.thd = NA, addp.thd = NA) {
    
    if(!("Marker" %in% colnames(targets))) {
      stop("Marker column not found in targets!")
    }
    if(!all(c("array_element_name", "locus") %in% colnames(p2g))) {
      stop("check probe to gene ID conversion table")
    }
    
    mks <- factor(targets$Marker)
    queries <- queries[row.names(queries) %in% mks, colnames(queries) %in% mks]
    
    if(sva) {
      design <- model.matrix( ~ mks)
      colnames(design) <- c("Intercept", levels(mks))
      n.sv <- num.sv(exprs(eset), design, method = "leek")
      svobj <- sva(exprs(eset), design, n.sv=n.sv)
      design <- cbind(design, svobj$sv)
    } else {
      design <- model.matrix( ~ mks)
      colnames(design) <- c("Intercept", levels(mks))
      svobj <- NULL
    }
    
    fit <- lmFit(eset, design)
    
    cat("  calculate enriched genes in each cell type...\n\t")
    cty.probes <- list()
    for (cty in levels(mks)) {
      cat(sprintf("%s, ", cty))
      temp <- CellSpecificGenes(cty, fit, queries = queries, svobj = svobj, 
                                addFC.thd = addFC.thd, addp.thd = addp.thd)
      cty.probes[[cty]] <- calculate.overlap2(temp)
    }
    cat("\n")
    
    ## convert probeID to geneID
    cat(" Convert probeID to geneID...\n")
    cty.genes <- lapply(cty.probes, function(x) {
      p2g[match(x, p2g$array_element_name), ]$locus})
    cty.genes <- lapply(cty.genes, function(x) {
      unlist(sapply(x, function(g) unlist(str_split(g, ";"))))})
    cty.genes <- lapply(cty.genes, function(x) x[!x == "no_match"])
    cty.genes <- lapply(cty.genes, unique)
    emptymarker <- which(sapply(cty.genes, is.null))
    if(any(emptymarker)) cty.genes[[emptymarker]] <- " "
    
    ## convert result to a dataframe and save to a file
    cat(" Save cell type specific genes to file...\n")
    if(is.na(out.file)) {
       enrdfile <- sprintf("EnrdGeneList.%s.fc%1.1f.p%s.csv", 
                           mthd, fc.thd, 
                           formatC(p.thd, format = "e", digits = 0))
       out.file <- file.path(out.path, enrdfile)
    }
    n.obs <- sapply(cty.genes, length)
    seq.max <- seq_len(max(n.obs))
    mat <- sapply(cty.genes, "[", i = seq.max)
    mat[is.na(mat)] <- ""
    df.cty.genes <- as.data.frame(mat)
    cat("   Write result in: ", out.file, "\n")
    write.table(df.cty.genes, out.file, sep = ",", 
                row.names = F, col.names = T, quote = F)
    
    return(cty.genes)
}

calculate.overlap2 <- function(lst) {
  if(length(lst) < 0) stop("empty input")
  
  if(length(lst) == 1) {
    names(lst) <- NULL
    return(unname(unlist(lst)))
  }
  
  if(length(lst) > 1) {
    itst <- intersect(lst[[1]], lst[[2]])
    lst[[2]] <- itst
    names(lst)[2] <- "Overlapped"
    lst[[1]] <- NULL
    lst <- calculate.overlap2(lst)
  }
}

ma.PCA <- function(madata, ntop = 1000, title = "") {
  
  img <- list(px = 12)
  mat4pca <- as.matrix(exprs(madata))
  colnames(mat4pca) <- sampleNames(madata)
  
  Pvars<- rowVars(mat4pca)
  select <- order(-Pvars)[seq_len(min(ntop, length(Pvars)))]
  PCA <- prcomp(t(mat4pca[select,]), scale = F)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  dataGG <- data.frame(PC1=PCA$x[,1], PC2=PCA$x[,2], 
                       PC3=PCA$x[,3], PC4=PCA$x[,4],
                       sampleNo = colnames(mat4pca),
                       sample = phenoData(madata)@data$Cell.type)
  cols <- c(brewer.pal(9, "Set1"),brewer.pal(12, "Set3"), 
            brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
  cols <- cols[1:nlevels(as.factor(dataGG$sample))]
  
  PCAplot <- ggplot(dataGG, aes(PC1, PC2, label = sampleNo)) +
    ggtitle(title) + 
    geom_point(aes(col=factor(sample)),size=5*img$px/14) +
    scale_colour_manual(values = cols) + 
    geom_text_repel(size = 3*img$px/14,
                    color = "gray60",
                    hjust="inward", vjust = "inward", 
                    show.legend = F) +
    labs(title = title,
         x = paste0("PC1, VarExp:", round(percentVar[1], 4)),
         y = paste0("PC2, VarExp:", round(percentVar[2], 4))) + 
    
    theme(#text = element_text(size = img$px),
      axis.title.x = element_text(colour = "black",size = img$px),   # hide x axis title.
      axis.ticks.x = element_line(size = img$lsz * 0.8, color = "black"), #hide x axis ticks.
      axis.title.y = element_text(color = "black", size = img$px),
      axis.ticks.y = element_line(size = img$lsz * 0.8, color = "black"),
      axis.text.x = element_text(vjust = 1, hjust = 1, color = "black", size = img$px),
      axis.text.y = element_text(face = "plain", color = "black", size = img$px),
      
      legend.background = element_blank(),
      legend.key = element_rect(fill="white", color="white"),
      legend.key.size = unit(img$px*0.8,"points"),
      legend.text = element_text(size = img$px * 0.8, color = "black"),
      legend.spacing.y = unit(img$px/5, "points"),
      legend.title = element_blank(),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = img$lsz),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0, size = img$px))
  
  return(PCAplot)
  
}

readAffy2 <- function(path, targets, p2g, default.rownames = F) {
 # targets <- targets[order(targets$Marker), ]
  files <- file.path(path, targets$File.name) 
  cat(sprintf("reading microarray data: %d files\n", length(files)))
  ab.raw <- ReadAffy(filenames = files)
  
  if(!default.rownames) {
    if(!is.null(targets$Treatment) & nlevels(factor(targets$Treatment)) > 1) {
      row.names(targets) <- paste(targets$Marker, targets$Treatment, targets$Rep, 
                                  sep = ".")  
    } else {
      row.names(targets) <- paste0(targets$Marker, ".", targets$Rep)
    }
  }
                               
  mt <- match(sampleNames(ab.raw), targets$File.name)
  phenoData(ab.raw) <- new("AnnotatedDataFrame", data = targets[mt, ])
  protocolData(ab.raw) <- new("AnnotatedDataFrame", data = targets[mt, ])
  row.names(p2g) <- p2g$array_element_name
  fmt <- match(row.names(ab.raw), row.names(p2g))
  fData(ab.raw) <- p2g[fmt,]
  
  return(ab.raw)
}

cty.init <- function(path, targets) {
  
  if(!file.exists(path)) dir.create(path)
  
  ## make sub-directories
  #subdirs <- c("vsBrady", "vsEIS", "vsLCM", "vsEnL", "logs", "qc.path")
  subdirs <- c("vsBrady", "logs", "qc.path")
  for (d in 1:length(subdirs)) {
    assign(subdirs[d], file.path(path, subdirs[d]), envir = .GlobalEnv)
    if(!file.exists(get(subdirs[d]))) dir.create(get(subdirs[d]))
  }

  assign("logfile", file.path(logs, "log.txt"), envir = .GlobalEnv)
  write(sprintf("\n## Running date: %s", Sys.Date()), logfile, append = T)
  
  # prepare syncytium genes
  EIS.wbc.up <- EIS.wbc %>%
    filter(log2FoldChange > log2(fc.thd), padj < p.thd) %>%
    droplevels()
  
  LCM.wbc.up <- LCM.wbc %>%
    filter(log2FoldChange > log2(fc.thd), padj < p.thd) %>%
    droplevels()
  
  wbc.up <- list(EIS = EIS.wbc.up$Gene_id,
                 LCM = LCM.wbc.up$Gene_id)
 # wbc.up$EnL <- intersect(wbc.up$EIS, wbc.up$LCM)
  
  EIS.vbc.up <- EIS.vbc %>%
    filter(log2FoldChange > log2(fc.thd), padj < p.thd) %>%
    droplevels()
  LCM.vbc.up <- LCM.vbc %>%
    filter(log2FoldChange > log2(fc.thd), padj < p.thd) %>%
    droplevels()
  
  vbc.up <- list(EIS = EIS.vbc.up$Gene_id,
                 LCM = LCM.vbc.up$Gene_id)
  #vbc.up$EnL <- intersect(vbc.up$EIS, vbc.up$LCM)
  
  ## gene universe
  ma.genes <- probe2gene$locus
  ma.genes <- ma.genes[!ma.genes == "no_match"]
  ma.genes <- sapply(ma.genes, function(x) unlist(str_split(x, ";")))
  ma.genes <- unname(unlist(ma.genes))
  gene_universe <- intersect(intersect(EIS.wbc$Gene_id, 
                                       LCM.wbc$Gene_id), ma.genes)
  
  assign("wbc.up", wbc.up, envir = .GlobalEnv)
  assign("vbc.up", vbc.up, envir = .GlobalEnv)
  assign("gene_universe", gene_universe, envir = .GlobalEnv)
}

comp2Brady <- function(genelist, Bradylist) {
  if(!file.exists(vsBrady)) dir.create(vsBrady)
  vsBrady.stats.file <- file.path(vsBrady, sprintf("overlap.rate%s.csv", suffix))
  
  if(!("ARR15" %in% names(genelist))) {
    genelist[["ARR15"]] <- c(genelist[["ARRI"]], 
                             genelist[["ARRII"]], 
                             genelist[["ARRIII"]])
  }
  genelist <- lapply(genelist, unique)
  cols <- brewer.pal(8, "Set2")[1:2]
  
  sumlist <- list()
  for (marker in names(genelist)){
    if (!(marker %in% names(Bradylist))) next
    b <- Bradylist[[marker]]
    a <- genelist[[marker]]
    vlist <- list(XL = a, BSM = b)
    o <- calculate.overlap(vlist)
    o.n <- unlist(sapply(o, length))
    sumlist[[marker]] <- c(o.n[1],o.n[2],o.n[3], round(100*o.n[3]/o.n[1],1))
    
    title <- sprintf("%s [XL.vs.BSM]", marker)
    fn <- sprintf("%s.png", marker)
    venn.diagram(vlist, 
                 cat.pos = 180, 
                 fill = cols,
                 main = title,
                 file.path(vsBrady, fn))
  }
  logs <- list.files(path = vsBrady, pattern = "*.log")
  file.remove(file.path(vsBrady, logs))
  
  cat("writing result...\n")
  vsBrady.rate <- t(as.data.frame(sumlist))
  colnames(vsBrady.rate) <- c("XL", "BSM", "Overlap", "O.R")
  write.table(vsBrady.rate, vsBrady.stats.file, sep = ",", 
              row.names = T, col.names = NA, quote = F)
}

# comp2Syncytium: to compare cell type enriched gene list to Syncytium enrich genes
#   genelist: list; one entry to one cell type
#   Syn.up: list; EIS, LCM, EnL
#   gene_universe: 
#   outpath: optional; if not assigned, will use the globe variable
# return: list of data frames
#         1. sum:    summary of overlapped genes
#         2. genes:  list of overlapped genes (EIS&LCM&Celltype-enriched)
#         3. fisher: p.values of fisher's exact test
comp2Syncytium <- function(genelist, Syn.up, gene_universe, outpath = NA) {
  
  if(!is.na(outpath)) {
    out.path <- outpath
    if(!file.exists(out.path)) dir.create(out.path)
  }
  
  e.path <- file.path(out.path, "euler")
  if(!file.exists(e.path)) dir.create(e.path)
  # go.path <- file.path(out.path, "go")
  ## consulidate ARRs to ARR15
  # if(!("ARR15" %in% names(genelist))){
  #   arr <- names(genelist)[which(grepl("ARRI", names(genelist)))]
  #   genelist[["ARR15"]] <- unique(unlist(genelist[arr]))
  # }
       
  genelist <- lapply(genelist, function(x) x[x %in% gene_universe])
  genelist <- lapply(genelist, unique)
  Syn.up <- lapply(Syn.up, function(x) x[x %in% gene_universe])
  
  cols <- brewer.pal(8, "Set2")[1:2]
  n.overlap <- list()
  # n.overlap[["Total"]] <- unlist(lapply(genelist, length))
  
  cat("fisher test and euler plot...\n\t")
  fisher.p <- list()
  overlapped.genes <- list()
  for(marker in names(genelist)) {
    cat(marker, ", ")
    n.overlap[[marker]]["Total"] <- length(genelist[[marker]])
    
    for (sample in names(Syn.up)) {
      vlist <- list(Syncytium = Syn.up[[sample]],
                    x = genelist[[marker]])
      names(vlist)[2] <- marker
      
      o <- calculate.overlap(vlist)
      
      ## Fisher's exact test
      cont.table <- data.frame(matrix(nrow = 2, ncol = 2))
      cont.table[1,1] <- length(o$a3)
      cont.table[1,2] <- length(o$a1) - length(o$a3)
      cont.table[2,1] <- length(o$a2) - length(o$a3)
      cont.table[2,2] <- length(gene_universe) - length(o$a2) - cont.table[1,2]
      ftest <- fisher.test(cont.table, alternative = "greater")
      fisher.p[[marker]][sample] <- ftest$p.value
      
      # overlapped gene summary
      n.overlap[[marker]][sample] <- length(o$a1)
      n.overlap[[marker]][sprintf("in%s", sample)] <- length(o$a3)
      n.overlap[[marker]][sprintf("perIn%s", sample)] <- 
        sprintf("%1.1f%%", length(o$a3)/length(o$a2) * 100)
    }
    
    # euler plot
    title <- sprintf("%s.vs.Syncytium", marker)
    e.fn <- file.path(e.path, sprintf("%s.pdf", title))
    elist <- c(Syn.up[c(1:2)], list(m = genelist[[marker]]))
    
    names(elist)[length(elist)] <- marker
    
    e.col <- brewer.pal(length(elist), "Set2")
    eplt <-  plot(euler(elist), 
           fill = colRamp(e.col,1.5),
           edges = colRamp(e.col,1),
           quantities = T, 
           labels = T )#identical(legend, T))
    
    # has to plot() again, otherwise would plot in a loop
    pdf(e.fn)
    plot(eplt)
    graphics.off()
    
    overlapped.genes[[marker]] <- calculate.overlap2(elist)
  } # end of outter loop (markers)
  cat("\n")
  cat("saving overlapped genes - vsEISvsLCM\n")
  n.obs <- sapply(overlapped.genes, length)
  seq.max <- seq_len(max(n.obs))
  mat <- sapply(overlapped.genes, "[", i = seq.max)
  mat[is.na(mat)] <- ""
  df.overlapped.genes <- as.data.frame(mat)
  write.table(df.overlapped.genes, 
              file.path(out.path, sprintf("[List].OverlappedGenes_%s&MK.csv", 
                                          paste0(names(Syn.up), collapse = "&"))), 
              sep = ",", row.names = F, col.names = T, quote = F)
  
  cat("Saving fisher test results...\n")
  df.sum <- t(as.data.frame(n.overlap))
  df.fisher <- t(as.data.frame(fisher.p))
  
  sumfile <- "OverlappedGenesSummary.csv"
  fisherfile <- "FishersExactTest.csv"
  write.table(df.sum, file.path(out.path, sumfile), 
              row.names = T, col.names = NA, sep = ",", quote = F)
  write.table(df.fisher, file.path(out.path, fisherfile), 
              row.names = T, col.names = NA, sep = ",", quote = F)
 
  
  return(list(fisher = df.fisher, sum = df.sum, genes = df.overlapped.genes))
}

# comp2Syncytium: to compare cell type enriched gene list to Syncytium enrich genes
# genelist: list; one entry to one cell type
# Syn.up: list; EIS, LCM, EnL
# gene_universe: 
# outpath: optional; if not assigned, will use the globe variable
comp2Syncytium.old <- function(genelist, Syn.up, gene_universe, outpath = NA) {
  
  if(!is.na(outpath)) {
    out.path <- outpath
    if(!file.exists(out.path)) dir.create(out.path)
  }
  
  ## consulidate ARRs to ARR15
  if(!("ARR15" %in% names(genelist))){
    arr <- names(genelist)[which(grepl("ARRI", names(genelist)))]
    genelist[["ARR15"]] <- unique(unlist(genelist[arr]))
  }
       
  genelist <- lapply(genelist, function(x) x[x %in% gene_universe])
  genelist <- lapply(genelist, unique)
  Syn.up <- lapply(Syn.up, function(x) x[x %in% gene_universe])
  
  cols <- brewer.pal(8, "Set2")[1:2]
  n.overlap <- list()
  n.overlap[["Total"]] <- unlist(lapply(genelist, length))
  fisher.p <- list()
  for(syn in names(Syn.up)) {
    if(is.na(outpath)) {
      path <- get(sprintf("vs%s", syn))
    } else {
      path <- file.path(out.path, sprintf("vs%s", syn))
    }
    if(!file.exists(path)) dir.create(path)
    
    n.overlap[[syn]] <- rep(length(Syn.up[[syn]]), length(genelist))
    name.o <- sprintf("In%s", syn)
    
    for (marker in names(genelist)) {
      vlist <- list(Syn.up = Syn.up[[syn]], x = genelist[[marker]])
      names(vlist) <- c("Syncytium(up)", marker)
      
      title <- sprintf("%s.vs.Syncytium(%s)", marker, syn)
      vname <- file.path(path, sprintf("%s.png", title))
      venn.diagram(vlist, 
                   cat.pos = 180, 
                   fill = cols,
                   main = title,
                   filename = vname)
      
      ## summary 
      o <- calculate.overlap(vlist)
      n.overlap[[name.o]][marker] <- length(o$a3)

      ## Fisher's exact test
      cont.table <- data.frame(matrix(nrow = 2, ncol = 2))
      cont.table[1,1] <- length(o$a3)
      cont.table[1,2] <- length(o$a1) - length(o$a3)
      cont.table[2,1] <- length(o$a2) - length(o$a3)
      cont.table[2,2] <- length(gene_universe) - length(o$a2) - cont.table[1,2]
      ftest <- fisher.test(cont.table, alternative = "greater")
      fisher.p[[syn]][marker] <- ftest$p.value
    }
    pct <- sprintf("%1.1f%%", n.overlap[[name.o]]/n.overlap[["Total"]] * 100)
    n.overlap[[sprintf("perIn%s",syn)]] <- pct
    
    logs <- list.files(path = path, pattern = "*.log")
    file.remove(file.path(path, logs))
  }
  ## euler venn
  cat("euler venn plot...\n")
  glist <- list()
  e.path <- file.path(out.path, "vsEISvsLCM")
  if(!file.exists(e.path)) dir.create(e.path)
  for(mk in names(genelist)) {
    # mk <- mks[1]
    vlist <- list(EIS = Syn.up[["EIS"]],
                  LCM = Syn.up[["LCM"]],
                  m = genelist[[mk]])
    names(vlist)[3] <- mk
    
    fn <- file.path(e.path, sprintf("%s.png", mk))
    png(fn, width = 480, height = 480)
    plot(euler(vlist), 
         fill = colRamp(brewer.pal(3, "Set2"),1.5),
         edges = colRamp(brewer.pal(3, "Set2"),1),
         quantities = T, 
         labels = T )#identical(legend, T))
    dev.off()
    #graphics.off()
    glist[[mk]] <- calculate.overlap(vlist)$a5
  }
  cat("saving overlapped genes - vsEISvsLCM\n")
  n.obs <- sapply(glist, length)
  seq.max <- seq_len(max(n.obs))
  mat <- sapply(glist, "[", i = seq.max)
  mat[is.na(mat)] <- ""
  ov3.genes <- as.data.frame(mat)
  write.table(ov3.genes, file.path(e.path, "[List].Overlapped.Gene.csv"), 
              sep = ",", row.names = F, col.names = T, quote = F)
  
  cat("Saving fisher test and hypergeomic test results...\n")
  df.sum <- as.data.frame(n.overlap)
  df.fisher <- as.data.frame(fisher.p)
  
  sumfile <- sprintf("Overlapped.Genes.Summary%s.csv", suffix)
  fisherfile <- sprintf("FishersExactTest%s.csv", suffix)
  write.table(df.sum, file.path(out.path, sumfile), 
              row.names = T, col.names = NA, sep = ",", quote = F)
  write.table(df.fisher, file.path(out.path, fisherfile), 
              row.names = T, col.names = NA, sep = ",", quote = F)
  
  return(list(sum = df.sum, fisher = df.fisher))
}

fisherPlot <- function(df, cols = NA) {
  
  if(any(is.na(cols))) cols <- brewer.pal(ncol(df), "Set2")
  scores <- -log10(as.matrix(df))
  rprod <- apply(scores, 1, function(x) prod(x) ** (1/length(x)))
  ord <- order(-rprod)
  res <- melt(scores)
  names(res) <- c("Marker", "Sample", "Score")
  res$Cell <- mk2cty$Map0[match(res$Marker, mk2cty$Marker)]
  res$Cell <- sprintf("%s (%s)", res$Cell, res$Marker)
  #res$score <- -log10(res$fisher.p)
  #res <- res[order(res$score, decreasing = T),]
  res$Cell <- factor(res$Cell, levels = unique(res$Cell)[ord])
  
  p <- ggplot(res, aes(x=Cell, y=Score, group = Sample, fill = Sample)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_y_continuous(expand=c(0, 0)) +
    xlab("Cell Types") + 
    ylab("-log10(pvalue)") + 
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") + 
    annotate("text", x = nlevels(res$Cell), y=2, hjust = -0.1,
             label = "p = 0.01", 
             size = 5,
             angle = 0, 
             col = colRamp("red", 1.5)) +
    scale_fill_manual(values = cols) + 
    coord_flip() +
    theme_bw()
  
  return(p)
}

# save final targets and querymatrix
saveSettings <- function() {
  
  path <- logs
  dt <- gsub("[-: ]", "", Sys.time())
  if(!file.exists(path)) dir.create(path)
  write.table(cel.targets, file.path(path, "targets_full.txt"), sep = "\t", quote = F)
 # write.table(tar, file.path(path, "targets.txt"),sep = "\t", quote = F)
 # write.table(query.mat, file.path(path, "queries.txt"), sep = "\t", quote = F)
  
  # copy all souced files to logs folder
  sapply(source_files, function(x) file.copy(file.path("./scripts", x), 
                                             file.path(path, paste0(dt, x))))
  # copy master script to logs folder
  file.copy(file.path("./scripts", thisfile), 
            file.path(path, paste(dt, thisfile, sep = "_")))
  # write running info file
  assign("readme", file.path(path, ".readme.txt"), envir = .GlobalEnv)
  apd <- ifelse(file.exists(readme), T, F)
  write(sprintf("\n## Running Date: %s", Sys.Date()), 
          file = readme, append = apd)
  
  # variables
  gse <- paste0(levels(as.factor(tar$GSE)),collapse = ", ")
  mks <- paste0(levels(as.factor(tar$Marker)),collapse = ", ")
  write(sprintf("Dataset:\t%s", gse), 
        file = readme, append = T)
  write(sprintf("Markers:\t%s", mks), 
        file = readme, append = T)
  write(sprintf("Normalization Method:\t%s", mthd) , 
        file = readme, append = T)
  write(sprintf("Foldchange threshold:\t%.2f", fc.thd), 
        file = readme, append = T)
  write(sprintf("adj.p threshold:\t%s", formatC(p.thd, format="e", digits = 0)), 
        file = readme, append = T)
  write(sprintf("fisher_test p threshold:\t%s", 
                formatC(fisher.p.thd, format="e", digits = 0)), 
        file = readme, append = T)
  
}

writeSessionInfo <- function(file) {
  #readme <- file.path(path, ".readme.txt")
  write("\n## SessionInfo ================================================",
        file = file, append = T)
  
  ssinfo <- sessionInfo()
  for(i in names(ssinfo)) {
    
    if(i %in% c("otherPkgs", "loadedOnly")) {
       title <- ifelse(i == "otherPkgs", 
                       "other attached packages:", 
                       "loaded via a namespace (and not attached):")
       write(title, file = file, append = T)
       v <- unlist(lapply(ssinfo[[i]], function(x) x[["Version"]]))
       pkgInfo <- sprintf("%s_%s", names(v), v)
       for(n in 1:ceiling(length(pkgInfo)/5)) {
          sublist <- pkgInfo[1:5+5*(n-1)]
          write(sprintf("\t%s", paste0(sublist, collapse = "\t")),
                file = file, append = T) 
       }
       
    } else if (i == "basePkgs") {
       bpkgs <- paste0(ssinfo[[i]], collapse = "\t")
       write("attached base packages: ", file = file, append = T)
       write(bpkgs, file = file, append = T)
    } else if(i == "R.version") {
       write(sprintf("%s: %s", i,ssinfo[[i]][["version.string"]]),
             file = file, append = T)
    } else {
       write(sprintf("%s:\t%s", i, ssinfo[[i]]), file = file, append = T)
    }
    write("", file = file, append = T)
  }
}
