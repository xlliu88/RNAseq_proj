## code need clean up

## Cell type Signature. master script
## cell specific genes after sva fitting
## major restructure of files
## added PCA with pcromp function
## add go
## change venn to euler
## will merge qc pdfs to single pdfs
## svd use to centerized data

## re-run 20200208_ARR15[gcrma_fc1.2-p1e03] set
## add set for non-overlapping cell-types


## import files
source_files <- c("_loadPackages.r", 
                  "_DESeq2Processor.r", 
                  "_plotters.r", 
                  "_importResult.bill.r",
                  "_bill.vars.r", 
                  "_celltype.vars.r", 
                  "_celltype.fun.r", 
                  "_go_kegg.r")

sapply(file.path("./scripts", source_files), source)

###### to override variables  ###################

assign("p.thd", 0.001, envir = .GlobalEnv)
assign("fc.thd", 1.0, envir = .GlobalEnv)
assign("mthd", "gcrma", envir = .GlobalEnv)
prefix <- sprintf("Enrd.SVA3.[%s_fc%1.1f-p%sx]", mthd, fc.thd, formatC(p.thd, format = "e", digits = 0))
prefix <- gsub("-", "", prefix)
out.path <- file.path(proj.path, prefix)
if(!file.exists(out.path)) dir.create(out.path)

## modified 2/28/2021
qc_path <- file.path(out.path, "qc2")
if(!file.exists(qc_path)) dir.create(qc_path)

## set up additional variables
# GSE.rm <- c("GSE7631") #, "GSE21553", "GSE61408", "GSE35580") #, "GSE64253")
MK.rm <- c("A8", "None")

tar.tmp <- cel.targets %>%
  filter(!Marker %in% MK.rm, DuplicateOf == 0)%>%
  droplevels()  

tar.tmp$Substrate <- sub("sucrose", "Sucr", tar.tmp$Substrate)
tar.tmp$Substrate <- sub(" \\(3mM\\)", "", tar.tmp$Substrate)
  
## set plot colors
cols <- colorGenerator(50, "qual")
hmcol <- colorGenerator(255, type = "seq", theme = "YlOrRd")
ntop.hm <- 500
ntop.dend <- 500 

## nrow(expr.gc.full)
## read ma data and normalization
cat("Reading microarray data...\n")
ab.raw.full <- readAffy2(cel.path, tar.tmp, probe2gene)
eset.raw.full <- computeExprSet(ab.raw.full, pmcorrect.method="pmonly",
                           summary.method="avgdiff")
expr.raw.full <- log2(exprs(eset.raw.full) + 0.1)

  
# cty.init(out.path)
  dataset <- "Mock"
  keep <- which(tar.tmp$Treatment == "Mock")
  cat("Mocks dataset\n")
  suffix <- sprintf("[%s_before_sva]", dataset)
 
  start_time <- Sys.time()
  tar <- tar.tmp[keep,]
  mks <- levels(as.factor(tar$Marker))
  qmat <- query.mat[colnames(query.mat) %in% mks, row.names(query.mat) %in% mks]
  ab.raw <- ab.raw.full[, keep]
  expr.raw <- expr.raw.full[, keep]
  eset.gc <- call.exprs(ab.raw, algorithm = mthd)
  expr.gc <- exprs(eset.gc)
  
  sn <- sub("Mock.", "", colnames(expr.raw))
  sampleNames(eset.gc) <- sn
  colnames(expr.raw) <- sn
  colnames(expr.gc) <- sn
  
  pdata <- phenoData(eset.gc)@data
  marker <- factor(pdata$Marker)
  gse <- factor(pdata$GSE)
  cty <- pdata$Cell.type
  
  write.csv(expr.raw, 
            file.path(out.path, sprintf("Expression_matrix_%s_raw.csv", dataset)), 
            quote = F)
  write.csv(expr.gc, 
            file.path(out.path, sprintf("Expression_matrix_%s_%s.csv",dataset, mthd)), 
            quote = F)
  
  write.table(tar, file.path(logs, sprintf("targets_%s.txt", dataset)), sep = "\t", quote = F)
  write.table(qmat, file.path(logs, sprintf("queries_%s.txt", dataset)), sep = "\t", quote = F)
  pdf(file.path(logs, sprintf("queries_%s.pdf", dataset)), width = 7, height = 7)
  image(1:nrow(qmat), 1:ncol(qmat), z = qmat, col = hmcol, 
        xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  abline(v = 0:nrow(qmat) + 0.5, col = "gray")
  abline(h = 0:nrow(qmat) + 0.5, col = "gray")
  axis(1, 1:nrow(qmat), colnames(qmat), las = 2, col = "gray", cex.axis = 0.75)
  axis(2, 1:nrow(qmat), row.names(qmat), las = 2, col = "gray", cex.axis = 0.75)
  dev.off()
  
  ##################################
  message("Ploting expression data before SVA ...\n")
  ## boxplot of raw
  ti <- sprintf("Expression_raw.[%s]", dataset)
  fn <- paste0("a0.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 5, pointsize = 8)
  par(mar = c(10.5, 5.1, 2.1, 5.1))
  boxplot(expr.raw, las = 2, ylab = "log2(expression)", 
          cex = 0.5,
          cex.axis = 0.5, 
          outline = F,
          #horizontal = T,
          col = cols[gse],
          lwd = 0.5)
  #par(xpd = T)
  legend("topright",#x = nrow(tar) + 7, y = mean(expr.raw), 
         levels(gse), 
         col = cols[1:nlevels(gse)], 
         pch = 15, cex = 0.5, bg = "white")
  title(ti, line = -1, cex = 0.3)
  dev.off()
  
  ## boxplot of normalized expression
  ti <- sprintf("Expression_%s.%s", mthd, suffix)
  fn <- paste0("a1.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 5, pointsize = 8)
  par(mar = c(10.5, 5.1, 2.1, 5.1))
  boxplot(expr.gc, las = 2, ylab = "log2(expression)", 
          cex = 0.5,
          cex.axis = 0.5, 
          outline = F,
          col = cols[gse],
          lwd = 0.5)
  legend("topright",#x = nrow(tar) + 7, y = mean(expr.raw), 
         levels(gse), 
         col = cols[1:nlevels(gse)], 
         pch = 15, cex = 0.5, bg = "white")
  title(ti, line = -1, cex = 0.3)
  dev.off()
  
  ## sample correlation density
  corrR<- cor(expr.raw)
  corrN <- cor(expr.gc)

  # dendrogram and heatmap using sample correlation
  cat("plot dendrogram and heatmap, before SVA...\n")
  topv <- order(-rowVars(expr.gc))
  topv.r <- order(-rowVars(expr.raw))
  
  corrR <- cor(expr.raw[topv.r[1:ntop.dend], ])
  corrN <- cor(expr.gc[topv[1:ntop.dend], ])
  rn <- cty[match(row.names(corrR), row.names(pdata))]
  #row.names(corrR) <- rn
  #row.names(corrN) <- rn
    
  dend_labR <- sprintf("%s [%s]", colnames(corrR), rn)
  dendRcorr <- as.dendrogram(hclust(as.dist(1-corrR))) %>%
    set("labels_cex", 500/nrow(tar)/8) 
  dendRcorr.ord <- order.dendrogram(dendRcorr)
  brcol <- cols[as.numeric(marker[dendRcorr.ord])]
  dendRcorr <- dendRcorr %>%
    set("labels", dend_labR[dendRcorr.ord]) %>%
    set("labels_colors", brcol)
    
  dend_labN <- sprintf("%s [%s]", colnames(corrN), rn)
  dendNcorr <- as.dendrogram(hclust(as.dist(1-corrN))) %>%
    set("labels_cex", 500/nrow(tar)/8) 
  dendNcorr.ord <- order.dendrogram(dendNcorr)
  brcol <- cols[as.numeric(marker[dendNcorr.ord])]
  dendNcorr <- dendNcorr %>%
    set("labels", dend_labN[dendNcorr.ord]) %>%
    set("labels_colors", brcol)
  
  # plot dendrogram
  ti <- "dendrogram_from_corr_raw"
  fn <- paste0("b0.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 5, pointsize = 8)
  plot(dendRcorr, main = ti)#, horiz = T)
  dev.off()
  
  ti <- sprintf("dendrogram_from_corr_%s.%s", mthd, suffix)
  fn <- paste0("b1.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 5, pointsize = 8)
  plot(dendNcorr, main = ti)
  dev.off()
  
  ## heatmap of sample correlation
  row.names(corrR) <- rn
  ti  <- "heatmap(corr)_raw_all"
  fn <- paste0("ca0.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 8.5, pointsize = 8)
  heatmap.2(1-corrR, trace = "none", col = rev(hmcol),
            Colv = dendRcorr, Rowv = dendRcorr,
            keysize = 1, key.title = "", key.xlab = "", key.ylab = "",
            density.info = "none", main = ti)
  dev.off()
  
  row.names(corrN) <- rn
  ti  <- sprintf("heatmap(corr)_%s_all.%s", mthd, suffix)
  fn <- paste0("ca1.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 8.5, pointsize = 8)
  heatmap.2(1-corrN, trace = "none", col = rev(hmcol),
            Colv = dendNcorr, Rowv = dendNcorr,
            keysize = 1, key.title = "", key.xlab = "", key.ylab = "",
            density.info = "none",
            main = ti)
  dev.off()
  
  ## heatmap of top genes clustering
  ti <- sprintf("genecluster(top%d)_from_corr_raw", ntop.hm)
  fn <- paste0("ct0.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 8.5, pointsize = 8) 
  hmR <- heatmap.2(expr.raw[topv[1:ntop.hm],],
                   Colv = dendRcorr, 
                   trace = "none", 
                   dendrogram = "column", 
                   col = hmcol,
                   ColSideColors = cols[marker],
                   keysize = 1,
                   key.title = "",
                   key.ylab = "",
                   key.xlab = "",
                   density.info = "none",
                   main = ti)
  dev.off()
  
  ti <- sprintf("genecluster(top%d)_%s_top", ntop.hm, mthd)
  fn <- paste0("ct1.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 8.5, pointsize = 8) 
  hmN <- heatmap.2(expr.gc[topv[1:ntop.hm],],
                   Colv = dendNcorr, 
                   trace = "none", 
                   dendrogram = "column", 
                   col = hmcol,
                   ColSideColors = cols[marker],
                   keysize = 1,
                   key.title = "",
                   key.ylab = "",
                   key.xlab = "",
                   density.info = "none",
                   main = ti)
  dev.off()
  
  
  
  # fit without sva ---------------------------------------------------------
  message("linear fit before SVA ... \n")
  X <- model.matrix( ~ marker)
  fit <- lmFit(eset.gc, X)
  beta <- fit$coefficients
  
  ## residues of simple linear model
  Yhat <- t(X %*% t(beta))
  errN <- expr.gc - Yhat
  
  ti <- sprintf("residues %s.%s", mthd, suffix)
  fn <- paste0("ar1.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 5, pointsize = 8)
  par(mar = c(8.5, 2.1, 2.1, 2.1))
  boxplot(errN, las = 2, ylab = "residues", cex = 0.5, outpch = ".",
          col = cols[factor(tar$GSE)],
          outline = F,
          lwd = 0.5)
  legend("topright", levels(gse), 
         col = cols[1:nlevels(gse)], 
         pch = 15, 
         cex = 0.5, 
         bg = "white")
  title(ti, line = -1, cex.main = 0.7)
  dev.off()
  
  ## % variances of PC
  cat("Explorering source of variances... \n")
  s <- svd(expr.gc - rowMeans(expr.gc))
  D <- s$d
  var_explained <- D^2/sum(D^2)
  ti = sprintf("Variance_Explained_%s %s", mthd, suffix)
  pdf(file.path(qc.path, paste0("d1.", ti, ".pdf")), width = 7, height = 4, pointsize = 10)
  plot(var_explained[1:20] * 100, pch = 16, 
       main = ti,
       xaxt = "n",
       ylim = c(0, max(var_explained) * 110),
       ylab = "% variance explained",
       xlab = "Principal component")
  axis(side = 1, at = 1:20)
  text(x = 1:5, y = var_explained[1:5] * 100, 
       labels = sprintf("%.1f%%", var_explained[1:5] * 100), cex = 0.7,
       adj = c(0.5, -0.8))
  dev.off()
  
  
  ## source of variances
  var_source <- c("Marker", "Date", "GSE", "Treatment", "Age", "Substrate", "Ecotype", "Lab")
  for (vs in var_source) {
    lvls <- nlevels(factor(tar[[vs]]))
    if(lvls == 1) {
      var_source <- var_source[!var_source == vs]
      next
    }
    tar[[vs]] <- as.factor(tar[[vs]])
  }
   var4fit <- tar %>% 
      select(all_of(var_source))
   a <- vector(mode = "list", length = 4)
  
  pltn <- length(var_source) + 1
  pltr <- ifelse(pltn > 3, 2, 1)
  pltc <- ifelse(pltr == 1, pltn, ceiling(pltn/pltr))
  for(i in 1:4) {
    message(sprintf("fitting pc%i with variables", i))
    var4fit$pc <- s$v[, i]
    fml <- sprintf("pc ~ 1 + %s", str_c(var_source, collapse = "+")) %>% as.formula()
    a[[i]] <- aov(fml, data = var4fit)
    summary(a[[i]]) %>% print()
    
    fn <- sprintf("e1.PC%i_var_source_%s %s.pdf", i, mthd, suffix)
    pdf(file.path(qc.path, fn), width = pltc * 2.5, height = pltr * 2.5, pointsize = 8)
    par(mfrow = c(pltr, pltc), mar = c(5.1, 4.1, 1.1, 1.1))
    for(vs in var_source) {
      boxplot(split(s$v[, i], factor(tar[[vs]])), las = 2, outline = F, outpch = "")
      stripchart(split(s$v[, i], factor(tar[[vs]])),  method = "jitter", 
                 add = T, vertical = T, cex = 0.5, pch = 19, col = "gray30")
      title(vs, line = -1)
    }
    plot(var_explained[1:20] * 100, pch = 16,
         xaxt = "n",
         ylab = "% variance explained",
         xlab = "Principal component")
    axis(side = 1, at = 1:20)
    title(sprintf("PC%i (%.1f%%)", i, D[i]^2/sum(D^2) * 100), line = -1)
    points(i, var_explained[i] * 100, pch = 16, col = "red")
    dev.off()
  }
  
  ## pca colored by variables
  pltc <- ceiling((pltn-1)/pltr)
  xlab <- sprintf("PC1 (%.1f%% Var)", s$d[1]^2/sum(s$d^2) * 100)
  ylab <- sprintf("PC2 (%.1f%% Var)", s$d[2]^2/sum(s$d^2) * 100)
  fn <- sprintf("f1a.PCAIsvd)_by_factors_%s.%s.pdf", mthd, suffix)
  pdf(file.path(qc.path, fn), width = 2.5 * pltc, height = 2.5 * pltr)
  par(mfrow = c(pltr,pltc), mar = c(4.1,2.1,1.1,1.1))
  for(vs in var_source) {
    fac <- factor(tar[[vs]])
    legcol <- ifelse(nlevels(fac) > 15, 2, 1)
    plot(s$d[1]*s$v[, 1], s$d[2]*s$v[, 2], pch = 16, 
         col = cols[fac], cex = 0.5, cex.axis = 0.7, xlab = "", ylab = "")
         #xlab = xlab, ylab = ylab)
    legend("topright", levels(fac), pch = 19, xjust = 1, cex = 0.5,
           col = cols[1:nlevels(fac)], ncol = legcol)
    title(main = sprintf("colored_by_%s", vs), 
         line = 0.2, font.main = 1, cex.main = 0.8)
    title(xlab = xlab, ylab = ylab, cex.lab = 0.7, line = 2)
  }
  dev.off()
  
  ## pca plot from prcomp
  pca.fn <- sprintf("f1b.PCA(prcomp)_all_%s.%s.pdf", mthd, suffix)
  pca1 <- PCAplot(expr.gc, 
                  factor(tar$Marker), 
                  cols = cols, 
                  circle = T, 
                  title_suffix = suffix)
  ggsave(file.path(qc.path, pca.fn), width = 6.5, height = 4.5)
  pca.fn <- sprintf("f1c.PCA(prcomp)_top%d_%s.%s.pdf",ntop.dend, mthd, suffix)
  pca2 <- PCAplot(expr.gc, 
                  factor(tar$Marker), 
                  ntop = ntop.dend, 
                  cols = cols, 
                  circle = T, 
                  title_suffix = sprintf("%s[top%s]", suffix, ntop.dend))
  ggsave(file.path(qc.path, pca.fn), width = 6.5, height = 4.5)
  
  
  
  
# ==============================================================================  
# +++ sva fit ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ==============================================================================
  message("Linear fitting with SVA ...")
  suffix <- sub("before", "after", suffix)
  
  X <- model.matrix(~ marker)
  colnames(X) <- sub("marker", "", colnames(X))
  colnames(X)[which(colnames(X) == "(Intercept)")] <- "Intercept"
  svaobj <- sva(expr(eset.gc), X) #, mod0)   # same result with or without mod0
  message(sprintf("Number of surogate variables: %d", svaobj$n.sv))
  
  svaX <- cbind(X, svaobj$sv)
  colnames(svaX) <- c(levels(marker), paste0("sv", 1:ncol(svaobj$sv)))
  svafit <- lmFit(eset.gc, svaX)
  
  beta.s <- svafit$coefficients[, 1:ncol(X)]  # coef for original model
  signal <-t(svaX[, 1:ncol(X)] %*% t(beta.s))
  alpha.s <- svafit$coefficients[, (ncol(X) + 1):ncol(svaX)] # coef for batch effects
  batch <- t(svaX[, (ncol(X) + 1):ncol(svaX)] %*% t(alpha.s))
  
  err.sva <- expr.gc - signal - batch
  expr.sva <- expr.gc - batch
  expr.sva[which(expr.sva < 0)] <- 0
  
  message("saving expression matrix after sva...")
  write.csv(expr.sva, 
            file.path(out.path, sprintf("Expression_matrix_%s_%s_sva.csv", dataset, mthd)),
            quote = F)
  
  
  ## boxplot of expression after SVA
  cat("Ploting Expression, after SVA ...")
  ti <- sprintf("Expression_%s.%s", mthd, suffix)
  fn <- paste0("a2.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 5, pointsize = 8)
  par(mar = c(10.5, 5.1, 2.1, 5.1))
  boxplot(expr.sva, las = 2, ylab = "log2(expression)", 
          cex = 0.5,
          cex.axis = 0.5, 
          outline = F,
          col = cols[gse],
          lwd = 0.5)
  legend("topright",#x = nrow(tar) + 7, y = mean(expr.raw), 
         levels(gse), 
         col = cols[1:nlevels(gse)], 
         pch = 15, cex = 0.5, bg = "white")
  title(ti, line = -1, cex = 0.3)
  dev.off()
  
  ## residues boxplot
  ti <- sprintf("residues %s.%s", mthd, suffix)
  fn <- paste0("ar2.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 5, pointsize = 8)
  par(mar = c(8.5, 2.1, 2.1, 2.1))
  boxplot(err.sva, las = 2, ylab = "residues", cex = 0.5, outpch = ".",
          col = cols[factor(tar$GSE)],
          outline = F,
          lwd = 0.5)
  legend("topright", levels(gse), col = cols[1:nlevels(gse)], pch = 15, cex = 0.5, bg = "white")
  title(ti, line = -1, cex.main = 0.7)
  dev.off()

  ## residue correlation
  corr.err <- sapply(1:ncol(errN), function(x) cor(errN[, x], err.sva[, x]))
  pdf(file.path(qc.path, "ar3.residue_correlation_befor.vs.after_sva.pdf"), 
      width = 8, height = 5, pointsize = 8)
  par(mar = c(5.5, 3.5, 2.1, 1))
  plot(marker, corr.err, ylim = c(1, 1), las = 2, outline = F, 
       cex.axis = 0.9, cex.main = 1, font.main = 1,
       main = "correlation of residues befor and after sva fit",
       ylab = "", xlab = "")
  stripchart(split(corr.err, marker), method = "jitter", 
             vertical = T, add = T, pch = 19, col = "green4")
  dev.off()
  
  ### variance explained by PC
  ss <- svd(expr.sva - rowMeans(expr.sva))
  var_xpln.sva <- ss$d^2/sum(ss$d^2)
  ti <- sprintf("d2.Variance_Explained_%s.%s", mthd, suffix)
  pdf(file.path(qc.path, paste0(ti, ".pdf")), width = 6, height = 4, pointsize = 10)
  plot(var_xpln.sva[1:20] * 100, pch = 16, 
       main = ti,
       ylim = c(0, max(var_xpln.sva) * 110),
       ylab = "% variance explained",
       xlab = "Principal component")
  text(x = 1:5, y = var_xpln.sva[1:5] * 100, 
       labels = sprintf("%.1f%%", var_xpln.sva[1:5] * 100), cex = 0.7,
       adj = c(0.5, -0.8))
  dev.off()
  
  ## source of variances  
  var_source <- c("Marker", "Date", "GSE", "Treatment", "Age", "Substrate", "Ecotype", "Lab")
  for (vs in var_source) {
     lvls <- nlevels(factor(tar[[vs]]))
     if(lvls == 1) {
       var_source <- var_source[!var_source == vs]
       next
     }
     tar[[vs]] <- as.factor(tar[[vs]])
  }
  var4fit <- tar %>% 
    select(all_of(var_source))
  a <- vector(mode = "list", length = 4)
  pltn <- length(var_source) + 1
  pltr <- ifelse(pltn > 3, 2, 1)
  pltc <- ifelse(pltr == 1, pltn, ceiling(pltn/pltr))
  for(i in 1:4) {
    var4fit$pc <- ss$v[, i]
    message(sprintf("fitting pc%s with variables...", i))
    fml <- sprintf("pc ~ 1 + %s", str_c(var_source, collapse = "+")) %>% as.formula()
    # fit <- lm(fml, data = var4fit )
    a[[i]] <- aov(fml, data = var4fit )
    print(summary(a[[i]]))
    fn <- sprintf("e2.PC%i_var_source_%s %s.pdf", i, mthd, suffix)
    pdf(file.path(qc.path, fn), width = pltc * 2.5, height = pltr * 2.5, pointsize = 8)
    par(mfrow = c(pltr, pltc), mar = c(5.1, 4.1, 1.1, 1.1))
    
    for(vs in var_source) {
      boxplot(split(ss$v[, i], factor(tar[[vs]])), las = 2, outline = F, outpch = "")
      stripchart(split(ss$v[, i], factor(tar[[vs]])),  method = "jitter", 
                 add = T, vertical = T, cex = 0.5, pch = 19, col = "gray30")
      title(vs, line = -1)
    }
    plot(var_xpln.sva[1:20] * 100, pch = 16,
         xaxt = "n",
         ylab = "% variance explained",
         xlab = "Principal component")
    axis(side = 1, at = 1:20)
    title(sprintf("PC%i (%.1f%%)", i, ss$d[i]^2/sum(ss$d^2) * 100), line = -1)
    points(i, var_xpln.sva[i] * 100, pch = 16, col = "red")
    dev.off()
  }
  
  
  ## pca colored by variables
  pltc <- ceiling((pltn-1)/pltr)
  xlab <- sprintf("PC1 (%.1f%% Var)", ss$d[1]^2/sum(ss$d^2) * 100)
  ylab <- sprintf("PC2 (%.1f%% Var)", ss$d[2]^2/sum(ss$d^2) * 100)
  fn <- sprintf("f2a.PCA(svd)_by_factors_%s.%s.pdf", mthd, suffix)
  pdf(file.path(qc.path, fn), width = 3 * pltc, height = 3 * pltr)
  par(mfrow = c(pltr,pltc), mar = c(4.1,3.1,1.1,1.1))
  for(vs in var_source) {
    fac <- factor(tar[[vs]])
    legcol <- ifelse(nlevels(fac) > 15, 2, 1)
    plot(ss$d[1]*ss$v[, 1], ss$d[2]*ss$v[, 2], pch = 16, 
         col = cols[fac], cex = 0.8, cex.axis = 1, xlab = "", ylab = "")
    if(!vs == "Date"){
      legend("topright", levels(fac), pch = 19, xjust = 1, cex = 0.5,
             col = cols[1:nlevels(fac)], ncol = legcol)
    }
    title(main = sprintf("colored_by_%s", vs), 
         line = 0.2, font.main = 1, cex.main = 1)
    title(xlab = xlab, ylab = ylab, cex.lab = 1, line = 2)
  }
  dev.off()
# 
# PCA plot using prcomp
  pca.fn <- sprintf("f2b.PCA(prcomp)_all_%s.%s.pdf", mthd, suffix)
  pca1 <- PCAplot(expr.sva, factor(tar$Marker), cols = cols, circle = T, title_suffix = suffix)
  ggsave(file.path(qc.path, pca.fn), width = 6.5, height = 4.5)
  pca.fn <- sprintf("f2c.PCA(prcomp)_top%d_%s.%s.pdf",ntop.dend, mthd, suffix)
  pca2 <- PCAplot(expr.sva, factor(tar$Marker), ntop = ntop.dend, cols = cols, circle = T, title_suffix = suffix)
  ggsave(file.path(qc.path, pca.fn), width = 6.5, height = 4.5)
  
  
  
  ## heatmap and dendrogame using correlation
  topv.sva <- order(-rowVars(expr.sva))  
  corrS <- cor(expr.sva[topv.sva[1:ntop.dend], ])
  # row.names(corrS) <- rn
  dend_labS <- sprintf("%s [%s]", colnames(corrS), rn)
  dendScorr <- as.dendrogram(hclust(as.dist(1-corrS))) %>%
    set("labels_cex", 500/nrow(tar)/8)
    
  dendScorr.ord <- order.dendrogram(dendScorr)
  brcol <- cols[as.numeric(marker[dendScorr.ord])]
  dendScorr <- dendScorr %>%
    set("labels", dend_labS[dendScorr.ord]) %>%
    set("labels_colors", brcol)
  
  # plot dendrogame
  ti <- sprintf("dendrogram_from_corr_%s %s", mthd, suffix)
  fn <- paste0("b2.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 5, pointsize = 8)
  plot(dendScorr, main = ti)#, horiz = T)
  dev.off()
  
  # heatmap of sample correlation
  row.names(corrS) <- rn
  ti <- sprintf("heatmap(corr)_%s %s.pdf", mthd, suffix)
  fn <- paste0("ca2.", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 8.5, pointsize = 8) 
  hmS <- heatmap.2(1-corrS, trace = "none", col = rev(hmcol),
                   Colv = dendScorr, Rowv = dendScorr,
                   keysize = 1,
                   key.title = "",
                   key.ylab = "",
                   key.xlab = "",
                   density.info = "none",
                   main = ti)
  dev.off()
  
  # heatmap of gene clustering
  ti <- sprintf("genecluster(top%d)_%s %s.pdf", ntop.hm, mthd, suffix)
  fn <- paste0("ct2", ti, ".pdf")
  pdf(file.path(qc.path, fn), width = 8.5, height = 8.5, pointsize = 8) 
  hmS <- heatmap.2(expr.sva[topv.sva[1:ntop.hm],],
                   Colv = dendScorr, 
                   trace = "none", 
                   dendrogram = "column", 
                   col = hmcol,
                   ColSideColors = cols[marker],
                   keysize = 1,
                   key.title = "",
                   key.ylab = "",
                   key.xlab = "",
                   density.info = "none",
                   main = ti)
  dev.off()
  
  
  # before and after sva comparison -----------------------------------------
  ## correlation density
  ti <- sprintf("Sample Correlation Density %s", suffix)
  fn <- sprintf("as.SampleCorrelationDensity %s - %s.pdf", mthd, suffix)
  corrR <- cor(expr.raw)
  corrN <- cor(expr.gc)
  corrS <- cor(expr.sva)
  
  corr.d.max <- max(c(density(corrN)$y, 
                      density(corrR)$y, 
                      density(corrS)$y))

  pdf(file.path(qc.path, fn), width = 7, height = 7)
  plot(density(corrR), 
       main = ti, 
       cex = 0.5,
       xlim = c(0.4, 1), 
       ylim = c(0, corr.d.max), 
       lty = 3, lwd = 2)
  lines(density(corrN), col = "black", lty = 2, lwd = 1)
  lines(density(corrS), col = "black", lwd = 1)
  legend("topleft", legend = c("raw data", "normalized", "sva fitted"), 
         lty = c(3, 2, 1), lwd = c(2,1,1)) 
         # col = c("black", "blue", "red"))
  dev.off()
  
  ## hclust comparison using correlation
  dendNcorr <- dendNcorr %>%
    set("labels_cex", 600/nrow(tar)/8) %>%
    set("branches_lwd", 1)
  dendScorr <- dendScorr %>%
    set("labels_cex",600/nrow(tar)/8) %>%
    set("branches_lwd", 1)
  dlcor <- dendlist(dendNcorr, dendScorr) %>%
    untangle(method = "step1side")
  
  pdf(file.path(qc.path, sprintf("b3.dendrogram_from_corr_%s_before.vs.after_sva %s.pdf", mthd, suffix)), 
      width = 8.5, height = 8.5, pointsize = 8)
  tanglegram(dlcor, margin_inner = 8, fast = T, lwd = 1.5, 
             main_left = "Before SVA", 
             main_right = "After SVA")
  dev.off()
  
  
  ## merge qc files
  pdfs <- list.files(qc.path, pattern = "*.pdf$")
  pages <- unname(sapply(file.path(qc.path,pdfs), pdf_length))
  pdfs <- pdfs[which(pages == 1)]
  qc.fn <- sprintf("qc_plots_%s.pdf", dataset)
  pdf_combine(file.path(qc.path, pdfs), output = file.path(qc.path, qc.fn))
  for(f in pdfs) file.remove(file.path(qc.path, f)) #, file.path(qc.path, sprintf("%sx",f)))
  
  
  
# ################################################################################
# ### cell specific genes from microarray results ################################
# ################################################################################
  
  enrdGeneList.fn <- sprintf("[list].Cell-type_ENRD_Genes.%s%s.csv", prefix, dataset)
  cty.genes <- getCellSpecificGenes2(eset.gc, svaX, queries = query.mat, 
                               p2g = probe2gene, sva = T, svobj = svaobj,
                               out.file = file.path(out.path, enrdGeneList.fn))
  
  # bar plot number of genes in each cell type.
  ngenes <- sapply(cty.genes, length)
  mkAGI <- tar$Marker.AGI[match(names(ngenes), tar$Marker)]
  selfin <- sapply(1:length(ngenes), function(x) mkAGI[x] %in% cty.genes[[x]])
  inUniv <- sapply(mkAGI, function(x) x %in% gene_universe)
  
  bcol <- ifelse(selfin, "lightgreen", colRamp("red4", 1.8))
  bcol <- ifelse(tolower(mkAGI) == "unknown"|!inUniv, "gray70", bcol)
  ymax <- max(ngenes) + 100
  pdf(file.path(out.path, sprintf("Number_of_Genes_%s.pdf",dataset)), width = 8, height = 6.5)
  bx <- barplot(ngenes, las = 2, space = 0.25, col = bcol, 
                ylim = c(0, ymax), xaxt = "n",
                main = sprintf("Summary of Cell-type Enriched Genes (%s)", dataset),
                ylab = "Number of Genes")
  axis(side = 1, at = bx, las = 2, labels = names(ngenes), cex = 0.75)
  text(x = bx, y = ngenes, pos = 3, label = ngenes, cex = 0.75)
  abline(h = 0)
  legend("topright", legend = c("Self In", "Self not In", "Gene ID unknown\nor No Probe"),
         fill = c("lightgreen", colRamp("red4", 1.8), "gray70"), border = "black", pt.cex = 1.5)
  dev.off()
  
  ## compare cell type enriched genes to genes in Brady & Zhang's paper
  comp2Brady(cty.genes, enrd.brady)
 
##  ## 
## compare Cty Enrd genes to Syncytium up genes
## import Syncytium up genes.

wbc.up <- list(EIS = EIS.wbc,
               LCM = LCM.wbc) %>% 
  lapply(., function(x) filter(x, padj < 0.05, log2FoldChange > 0)) %>% 
  lapply(., function(x) select(x, Gene_id) %>% unlist %>% unname())
wbc.up$EnL <- intersect(wbc.up$EIS, wbc.up$LCM)  

vbc.up <- list(EIS = EIS.vbc,
               LCM = LCM.vbc) %>% 
  lapply(., function(x) filter(x, padj < 0.05, log2FoldChange > 0)) %>% 
  lapply(., function(x) select(x, Gene_id) %>% unlist %>% unname())
vbc.up$EnL <- intersect(vbc.up$EIS, vbc.up$LCM)  

for(gt in c("WT", "clv")) {
    ## celltype enrd genes vs Syncytium
    out.sub <- file.path(out.path, sprintf("Cty2%ssyn%s",gt, suffix))
    if(!file.exists(out.sub)) dir.create(out.sub)
    
    datlabel <- ifelse(gt == "WT", "wbc", "vbc") #wbc.up else vbc.up 
    syn.up <- get(sprintf("%s.up", datlabel))
    cty2syn <- comp2Syncytium(cty.genes, syn.up, gene_universe, outpath = out.sub)
    
    ## plot fisher test result
    f.plot <- fisherPlot(cty2syn[["fisher"]])
    ggsave(file.path(out.sub, sprintf("fisherPlot_celltype.vs.%ssyn %s.pdf", gt, suffix)), 
           device = "pdf", width = 8, height = 11, units = "in", dpi = 300)
    
    
    ## go enrichment
    go.path <- file.path(out.sub, "GO")
    if(!file.exists(go.path)) dir.create(go.path)
    
    fisher <- cty2syn[["fisher"]]
    for(r in 1:nrow(fisher)) {
      if(any(fisher[r, ] < fisher.p.thd)) {
        mk <- row.names(fisher)[r]
        genes <- cty2syn[["genes"]][[mk]]    
        
        x <- rep(1, length(gene_universe))
        names(x) <- gene_universe
        x[which(names(x) %in% genes)] <- 0.001
        
        go <- goTair(x, "BP", out.path = go.path, fn.prefix = mk)
        go <- goTair(x, "MF", out.path = go.path, fn.prefix = mk)
        go <- goTair(x, "CC", out.path = go.path, fn.prefix = mk)
        
        go.res <- readGO(go.path, patterns = mk)
        names(go.res) <- unname(sapply(names(go.res), 
                                       function(x) unlist(str_split(x, "\\."))[2]))
        
        for (nm in names(go.res)) {
          goVisBar(go.res[[nm]], img.path = go.path, pthreshold = fisher.p.thd, 
                   fn.prefix = sprintf("%s-%ssyn.%s", mk, gt, nm))
        }
        
        gene_anno <- geneInfo[geneInfo$Gene_id %in% genes, ]
        EISexp <- get(sprintf("EIS.%s", datlabel))
        gene.exp.EIS <- EISexp[match(gene_anno$Gene_id, EISexp$Gene_id),]
        
        LCMexp <- get(sprintf("LCM.%s", datlabel))
        gene.exp.LCM <- LCMexp[match(gene_anno$Gene_id, LCMexp$Gene_id),]

        
        gene_anno$EWR.lfc <- gene.exp.EIS$log2FoldChange
        gene_anno$EWR.p <- gene.exp.EIS$padj
        gene_anno$LCM.lfc <- gene.exp.LCM$log2FoldChange
        gene_anno$LCM.p <- gene.exp.LCM$padj
        
        fn <- sprintf("[Gene_Details].%s-%ssyn_Overlapped.csv", mk, gt)
        write.csv(gene_anno, file.path(out.sub, fn), quote = F)
      }
    } # end of loop for fisher score plot
    
} # end of loop for WT/clv genotypes


## Brady result vs Syncytium-up genes
b.WTout <- file.path(out.path, "Brady2WTsyn")
if(!file.exists(b.WTout)) dir.create(b.WTout)
b.clvout <- file.path(out.path, "Brady2clvsyn")
if(!file.exists(b.clvout)) dir.create(b.clvout)

brady2WTsyn <- comp2Syncytium(enrd.brady, wbc.up, gene_universe, outpath = b.WTout)
brady2clvsyn <- comp2Syncytium(enrd.brady, vbc.up, gene_universe, outpath = b.clvout)

p.brady.WT <- fisherPlot(brady2WTsyn[["fisher"]])
ggsave(file.path(b.WTout, "fisherPlot_Brady.vs.WTsyn_.pdf"),
       device = "pdf", width = 8, height = 11, units = "in", dpi = 300)
p.brady.clv <- fisherPlot(brady2clvsyn[["fisher"]])
ggsave(file.path(b.clvout, "fisherPlot_Brady.vs.clvsyn_.pdf"),
       device = "pdf", width = 8, height = 11, units = "in", dpi = 300)



######### Wrapup and save settings ####################################
## write settings and sessionInfo
thisfile <- basename(rstudioapi::getSourceEditorContext()$path)
saveSettings()
writeSessionInfo(readme)



## ATHB ovrelapped genes
ecol <- brewer.pal(8, "Set2")
dataset <- "Mock"
WTfn <- file.path(out.path, "Cty2WTsyn[Mock_after_sva]", "[Gene_Details].ATHB15-WTsyn_Overlapped.csv")
clvfn <- file.path(out.path, "Cty2clvsyn[Mock_after_sva]", "[Gene_Details].ATHB15-clvsyn_Overlapped.csv")
cty.Enrd.fn <- file.path(out.path, "[list].Cell-type_ENRD_Genes.[gcrma_fc1.2p1e02]Mock.csv")

athb15WT <- read.csv(WTfn, header = T, stringsAsFactors = F)
athb15clv <- read.csv(clvfn, header = T, stringsAsFactors = F)
cty.Enrd <- read.csv(cty.Enrd.fn, header = T, stringsAsFactors = F)

ATHB15 <- cty.Enrd$ATHB15
ATHB15 <- ATHB15[!ATHB15 == ""]

expr.list <- list(WT.EIS  = EIS.wbc,
                  clv.EIS = EIS.vbc,
                  WT.LCM = LCM.wbc,
                  clv.LCM = LCM.vbc)
synup <- lapply(expr.list, function(x) {
  df.up <- x[abs(x$log2FoldChange) > log2(fc.thd) & x$padj < p.thd, ]
  return(df.up$Gene_id)
})

elist <- c(synup, list(ATHB15 = ATHB15))
names(elist) <- c("WT-ATHB15", "clv-ATHB15")
eul <- euler(elist)

plot(eul, 
     fill = colRamp(ecol[1:2], 1.5), edges = ecol[1:2], 
     labels = T, quantities = T)


## compare to SD2009PJ result
sd.path <- "./syncytia2009"
sd.res <- read.csv(file.path(sd.path, "tableS1.csv"), header = T, stringsAsFactors = F)
fn.prefix <- sprintf("[fc%.1fp%s]", fc.thd, formatC(p.thd, format = "e", digits = 0))
names(sd.res)[1] <- "Gene_id"
sd.up <- sd.res[sd.res$M > log2(fc.thd) & sd.res$adj.q < p.thd, ]
sd.dn <- sd.res[sd.res$M < log2(fc.thd) & sd.res$adj.q < p.thd, ]

gene_universe.sd <- intersect(sd.res$Gene_id, intersect(EIS.wbc$Gene_id, LCM.wbc$Gene_id))
EIS.wbc.up <- EIS.wbc[EIS.wbc$log2FoldChange > log2(fc.thd) & EIS.wbc$padj < p.thd, ]
EIS.wbc.dn <- EIS.wbc[EIS.wbc$log2FoldChange < log2(fc.thd) & EIS.wbc$padj < p.thd, ]

LCM.wbc.up <- LCM.wbc[LCM.wbc$log2FoldChange > log2(fc.thd) & LCM.wbc$padj < p.thd, ]
LCM.wbc.dn <- LCM.wbc[LCM.wbc$log2FoldChange < log2(fc.thd) & LCM.wbc$padj < p.thd, ]

elist.sdup <- list(EIS = intersect(EIS.wbc.up$Gene_id, gene_universe.sd),
                   LCM = intersect(LCM.wbc.up$Gene_id, gene_universe.sd),
                   SD2009 = intersect(sd.up$Gene_id, gene_universe.sd))

elist.sddn <- list(EIS = intersect(EIS.wbc.dn$Gene_id, gene_universe.sd),
                   LCM = intersect(LCM.wbc.dn$Gene_id, gene_universe.sd),
                   SD2009 = intersect(sd.dn$Gene_id, gene_universe.sd))


eul.sdup <- euler(elist.sdup)
eul.sddn <- euler(elist.sddn)
pdf(file.path(sd.path, sprintf("SD2009vsEISvsLCM_up_euler.%s.pdf", fn.prefix)))
plot(eul.sdup, fill = colRamp(ecol[1:3],1.5),edge = ecol[1:3], quantities = T)
dev.off()

pdf(file.path(sd.path, sprintf("SD2009vsEISvsLCM_dn_euler.%s.pdf", fn.prefix)))
plot(eul.sddn, fill = colRamp(ecol[1:3],1.5),edge = ecol[1:3], quantities = T)
dev.off()

enrd.genes <- read.csv(file.path(out.path, "[list].Cell-type_ENRD_Genes.[gcrma_fc1.2p1e02]Mock.csv"), 
                       header = T, stringsAsFactors = F)
enrd.genes <- as.list(enrd.genes)
sapply(enrd.genes, length)
head(enrd.genes[[1]])
enrd.genes <- lapply(enrd.genes, function(x) x[!x==""])

cty2sd2009 <- comp2Syncytium(enrd.genes, list(sd = sd.up$Gene_id), 
                              gene_universe = intersect(gene_universe, gene_universe.sd), 
                      outpath = file.path(out.path, "Cty2SD2009[Mock_after_sva]"))

p.sd <- fisherPlot(cty2sd2009[["fisher"]])
ggsave(file.path(out.path, sprintf("fisherPlot_SD2009.vs.CtyEnrd.%s.pdf", fn.prefix )),
       device = "pdf", width = 8, height = 11, units = "in", dpi = 300)

athb15sd <- cty2sd2009$genes$ATHB15
athb15sd.tb <- sd.res[sd.res$Gene_id %in% athb15sd, ]
write.csv(athb15sd.tb, file.path(out.path, "Cty2SD2009[Mock_after_sva]", "[Gene_Detail].SD2009vsATHB15_overlapped.csv"),
            quote = F)

athb15sd.tb2 <- EIS.wbc[EIS.wbc$Gene_id %in% athb15sd, ]
write.csv(athb15sd.tb2, file.path(out.path, "Cty2SD2009[Mock_after_sva]", "[Gene_Detail].SD2009vsATHB15_overlapped_in_wbc.csv"),
          quote = F)

athb15WT <- read.csv(file.path(out.path, "Cty2WTsyn[Mock_after_sva]", "[Gene_Details].ATHB15-WTsyn_Overlapped.csv"), header = T, stringsAsFactors = F)
wt15sd <- sd.res[sd.res$Gene_id %in% athb15WT$Gene_id, ]

list4 <- list(EIS = EIS.wbc.up$Gene_id,
              LCM = LCM.wbc.up$Gene_id,
              SD = sd.up$Gene_id,
              ATHB15 = enrd.genes$ATHB15)
venn.diagram(list4, file.path(out.path, "venn_EISvsLCMvsSDvsATHB15.png"),fill = colRamp(ecol[1:4], 1.5))
o4 <- calculate.overlap(list4)
sapply(o4, length)
ATHB15WT_SDsyn <- o4$a6

df <- EIS.wbc[EIS.wbc$Gene_id %in% ATHB15WT_SDsyn, ]
df2 <- LCM.wbc[LCM.wbc$Gene_id %in% ATHB15WT_SDsyn, ]
df3 <- sd.res[sd.res$Gene_id %in% ATHB15WT_SDsyn, ]

write.csv(df, file.path(out.path, "[Gene_Details].ATHB15&WTsyn&SDsyn.overlapped[EIS].csv"), quote = F)
write.csv(df2, file.path(out.path, "[Gene_Details].ATHB15&WTsyn&SDsyn.overlapped[LCM].csv"), quote = F)
write.csv(df3, file.path(out.path, "[Gene_Details].ATHB15&WTsyn&SDsyn.overlapped[SD].csv"), quote = F)

x <- rep(1, length(gene_universe.sd))
names(x) <- gene_universe.sd
x[which(names(x) %in% ATHB15WT_SDsyn)] <- 0.001
go.path <- file.path(out.path, "Cty2SD2009[Mock_after_sva]", "GO")
go <- goTair(x, "BP", out.path = go.path, fn.prefix = "ATHB15&WTsyn&SDsyn")
go <- goTair(x, "MF", out.path = go.path, fn.prefix = "ATHB15&WTsyn&SDsyn")
go <- goTair(x, "CC", out.path = go.path, fn.prefix = "ATHB15&WTsyn&SDsyn")

go.res <- readGO(go.path, patterns = "ATHB15&WTsyn&SDsyn")
names(go.res) <- unname(sapply(names(go.res), 
                               function(x) unlist(str_split(x, "\\."))[2]))
for (nm in names(go.res)) {
  goVisBar(go.res[[nm]], img.path = go.path, pthreshold = fisher.p.thd, 
           fn.prefix = sprintf("ATHB15&WTsyn&SDsyn%s", nm))
}

cty15exp.in.sd <- sd.res[sd.res$Gene_id %in% enrd.genes$ATHB15, c("Sync", "Root", "M", "adj.q")]
cty15exp.in.sd$color <- ifelse(cty15exp.in.sd$adj.q < 0.05, "red", "gray70")
cty15exp.in.sd$baseMean <- (2**cty15exp.in.sd$Sync + 2**cty15exp.in.sd$Root)/2
m <- cty15exp.in.sd$M
p <- cty15exp.in.sd$adj.q
nup <- m > log2(fc.thd) & p < p.thd
ndn <- m < log2(fc.thd) & p < p.thd
col <- ifelse(p < p.thd, "red", "gray50")
col <- ifelse(ndn, "blue", col)

pdf(file.path(out.path, "Expression_of_ATHB15_Enriched_Genes_SD2009.pdf"), width = 6, height = 6)
plot(cty15exp.in.sd$baseMean, cty15exp.in.sd$M, pch = 16, col = col,
     ylab = "log2FoldChange", xlab = "Mean Expression", ylim = c(-6,6),
     main = "Expression of ATHB15 Enriched Genes in SD2009 Dataset")
abline(h=0)
text(x = 1500, y = 6, labels = sprintf("up-regulated: %d", sum(nup)), pos = 2)
text(x = 1500, y = -6, labels = sprintf("down-regulated: %d", sum(ndn)), pos = 2)
dev.off()


## expression of ATHB15&WTSyn overlapped genes in Micro-aspiration dataset
cty15WT.in.sd <- sd.res[sd.res$Gene_id %in% athb15WT$Gene_id, c("Sync", "Root", "M", "adj.q")]
cty15WT.in.sd$color <- ifelse(cty15WT.in.sd$adj.q < 0.05, "red", "gray70")
exp.s <- cty15WT.in.sd$Sync
exp.r <- cty15WT.in.sd$Root
bm <- (2**exp.s + 2**exp.r)/2
m <- cty15WT.in.sd$M
p <- cty15WT.in.sd$adj.q
nup <- m > log2(fc.thd) & p < p.thd
ndn <- m < log2(fc.thd) & p < p.thd
col <- ifelse(p < p.thd, "red", "gray50")
col <- ifelse(ndn, "blue", col)

pdf(file.path(out.path, "Expression_of_ATHB15&WTSyn_Genes_SD2009.pdf"), width = 6, height = 6)
plot(bm, m, pch = 16, col = col,
     ylab = "log2FoldChange", xlab = "Mean Expression", ylim = c(-6,6))
     #main = "Expression of ATHB15 Enriched Genes in SD2009 Dataset")
abline(h=0)
text(x = 1500, y = 6, labels = sprintf("up-regulated: %d", sum(nup)), pos = 2)
text(x = 1500, y = -6, labels = sprintf("down-regulated: %d", sum(ndn)), pos = 2)
dev.off()





