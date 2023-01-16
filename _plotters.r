## plot functions

qcPlots <- function(dataset, imgpath, filename = NULL, type, subset = "all", 
                    ref_levels = NA, cols="", col = "rev", 
                    use=c("display", "ppt", "publish", "poster"), title="") {
  
  # to generate qc plots.
  # including count density, cedf, sample distance heatmap, PCA, gene clusering heatmap

  ## check plot type and input object class
  if (!type %in% c("density", "ecdf", "heatmap", "PCA", "genecluster", "dispersion", "MA")) {
    stop("plot type unknown")
  }
  
  if (!is.na(ref_levels[1])) {
    dataset <- ddsFactorRelevel(dataset, ref_levels = ref_levels)
  }
  
  if (type %in% c("density", "ecdf")) {
    if (!class(dataset) == "DESeqDataSet") {
      stop("Wrong data type")
    }
    keep <- rowSums(counts(dataset)) >= 10
  }
  
  if (type %in% c("heatmap", "PCA", "genecluster")) {
    if (class(dataset) == "DESeqDataSet") {
      cat("DESeqDataSet detected; performaing vst log2 transformation on dataset\n")
      cat("It may take a few minutes.\n")
      vst.dataset <- vst(dataset, blind = T)
    } else if (class(dataset) == "DESeqTransform") {
      cat("DESeqTransform datatype detected.\n")
    } else {
      stop("Wrong data type\n")
    }
  }
  
  #image export settings
  if (use == "ppt") {
    img <- list(type = "pdf", w = 10, h = 8, unit = "in", dpi = 150, px = 20, lsz = 1.5, psz = 6)
  } else if (use == "publish") {
    img <- list(type = "pdf", w = 5, h = 4, unit = "in", dpi = 300, px = 12)
  } else if (use == "poster" ) {
    img <- list(type = "pdf", w = 5, h = 4, unit = "in", dpi = 150)
  } else if (use == "display" ){
    img <- list(type = NA, w = 5, h = 4, unit = "in", dpi = 300, px = 12)
  } else {
    stop("unknown usage")
  }
  
  if (!use == "display") {
    if(!is.null(filename)) {
      imgname <- filename
    } else {
      imgname <- sprintf("qc_%s_%s.%dx%din.%ddpi.pdf", type, use, img$w, img$h, img$dpi)
    }
    
    cat(imgpath)
    cat("\n")
    cat(imgname) 
  } else {
    imgpath <- NA
    imgname <- NA
  }

  fcts <- ddsGetFactors(dataset)
  dataset <- ddsCombineFactors(dataset, combined_name = "Group")
  
  fctlvls <- list()
  ncolors <- 1
  for (f in fcts) {
    fctlvls[[f]] <- levels(colData(dataset)[[f]])
    ncolors <- ncolors * length(fctlvls[[f]])
  }
  
  if(length(fctlvls[[1]])==2) {
    samplecol <- brewer.pal(10, "Paired")[1:min(ncolors,10)]
  } else if (length(fctlvls[[2]])==2) {
    samplecol <- c(brewer.pal(10, "Paired")[c(TRUE, FALSE)], brewer.pal(10, "Paired")[c(FALSE,TRUE)])[1:min(ncolors,10)]
  } else {
    samplecol <- brewer.pal(8, "Set2")[1:min(ncolors,8)]
  }
  if(col=="rev") {
    samplecol <- rev(samplecol)
  }
  labelcol <- colorMatch(colData(dataset)$Group, samplecol)
  
  ## set colors for plots
  hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(255)          ## color for heatmap
  gccol <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))  ## color for genecluster
  
  ## ploting
  if (type == "density" ) {
    #tiff(paste0(imgpath, imgname), width = img$w, height = img$h, units = img$unit, res = img$dpi, pointsize = img$px, compression = "lzw")
    multidensity(counts(dataset, normalized = T)[keep,], xlab="Normalized Counts", 
                 xlim=c(0,1000), col = labelcol, main = title)
    #dev.off()
    
  } else if (type == "ecdf" ) {
    multiecdf(counts(dataset, normalized = T)[keep,], xlab="Normalized Counts", 
              xlim=c(0,1000), col = labelcol, main = title)
    
  } else if (type == "heatmap") {
    dataset <- ddsCombineFactors(dataset, combined_name = "Group")
    distsRL <- dist(t(assay(vst.dataset)))
    mat<- as.matrix(distsRL)
    colnames(mat) <- colData(dataset)$Group
    rownames(mat) <- rownames(colData(dataset))
    
    pheatmap(mat, trace="none", 
             fontsize = img$px, angle_col = 90, 
             col=rev(hmcol), border_color = NA, 
             filename = ifelse(is.na(imgname), NA, file.path(imgpath, imgname)), 
             width = img$w, height = img$w * 0.9,
             main=title)
    
  } else if (type == "PCA") {
    
    if (subset == "all") {
      ntop=nrow(dataset)
      tisufix <- "all genes"
    } else {
      ntop=subset
      tisufix <- sprintf("top variable genes (%s)", ntop)
    }
    
    ti <- sprintf("PC1 vs PC2, %s", tisufix)
    
    Pvars<- rowVars(assay(vst.dataset))
    select <- order(Pvars, decreasing=T)[seq_len(min(ntop, length(Pvars)))]
    PCA <- prcomp(t(assay(vst.dataset)[select,]), scale = F)
    percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
    
    dataGG <- data.frame(PC1=PCA$x[,1], PC2=PCA$x[,2], PC3=PCA$x[,3], PC4=PCA$x[,4],
                         sampleNo = rownames(colData(dataset)),
                         sample = colData(dataset)$Group)
    
    PCAplot <- ggplot(dataGG, aes(PC1, PC2, label = colnames(dataset))) +
              ggtitle(title) + 
              geom_point(aes(col=factor(sample)),size=5*img$px/14) +
              scale_colour_manual(values = samplecol) + 
              geom_text_repel(size = 3*img$px/14,
                              color=colRamp(labelcol,1.5),
                              hjust="inward", vjust = "inward", 
                              show.legend = F) +
              labs(title = ti,
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

    if (!use=="display") {
        imgname <- sub(type, paste0(type, "(", ntop, ")"), imgname)
        ggsave(file.path(imgpath, imgname), device = img$type, scale = 1, width = img$w, height = img$w * 0.8, units = img$unit, dpi = img$dpi)   
      }
    
    return(PCAplot)
    
  } else if( type == "genecluster" ) {
        
        if (subset == "all") {
          ntop = nrow(dataset)
        } else {
          ntop = subset
        }
        
        topVar <- head(order(rowVars(counts(dataset, normalized=T)), decreasing = T),ntop)
        topVargenes <- as.matrix(counts(dataset, normalized=T)[topVar,])
        labels <- sub("(EWR|LCM)_", "", colnames(dataset))
        colnames(topVargenes) <- paste0(colData(dataset)$sample, "(", labels, ")") 
        title=sprintf("cluster of top %s variable genes", ntop)
        imgname <- sub(type, paste0(type, "(", ntop, ")"), imgname)
        #heatmap.2(topVargenes,scale="row", trace="none",dendrogram="column", col = gccol)
        pheatmap(topVargenes, 
                 scale="row", 
                 trace="none", 
                 col = gccol, 
                 border_color = NA,
                 angle_col = 90, 
                 fontsize = img$px,  
                 show_rownames = F,  
                 filename = ifelse(is.na(imgname), NA, file.path(imgpath, imgname)), 
                 width = img$w, height = img$w,
                 main=title)
      }
  
}

PCAplot <- function(data, facts, ntop = "all", dims = c(1,2), cols = NA, 
                    circle = F, showLab = T, labsize = 10,
                    title_suffix = "", theme = theme_bw()) {
    mat <- as.matrix(data)
    if(!is.numeric(mat)) stop("input has to be a numeric matrix or data frame")
    dims <- dims[1:2]
    ntop <- ifelse(ntop == "all", nrow(mat), ntop)
    facts <- as.factor(facts)
    if(is.na(cols)) cols <- colorGenerator(nlevels(facts))
    if(length(cols) < nlevels(facts)) {
      stop("not enough colors for factor levels")
    } else {
      cols <- cols[1:nlevels(facts)]
    }
    if(!ntop == "all") title_suffix <- sprintf("top%d most variable genes",ntop)
    title <- sprintf("PC%d vs PC%d (%s)", dims[1], dims[2], title_suffix)
    
    Pvars<- rowVars(mat)
    keep <- order(Pvars, decreasing=T)[seq_len(min(ntop, length(Pvars)))]
    PCA <- prcomp(t(mat[keep,]), scale = F)
    percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2), 1)
    
    dataGG <- data.frame(PC1=PCA$x[,dims[1]], PC2=PCA$x[,dims[2]], 
                         sampleNo = colnames(mat),
                         group = facts)
    
    p <- ggplot(dataGG, aes(PC1, PC2, label = colnames(mat), group = group, 
                            col = group, fill = group)) 
    ti <-  ggtitle(title) 
    pts <- geom_point(size=4*labsize/30) 
    color <- scale_colour_manual(values = cols) 
    
    if(showLab){
      text <- geom_text_repel(size = 4*labsize/30,
                            #color = "gray65",
                            hjust="inward", vjust = "inward", 
                            show.legend = F,
                            segment.size = 0.2,
                            segment.alpha = 0.5,
                            alpha = 0.5) 
    } else {
      text <- NULL
    }
    lab <- labs(title = ti,
               x = sprintf("PC%d (%.1f%% var explained)", 
                           dims[1], percentVar[dims[1]]),
               y = sprintf("PC%d (%.1f%% var explained)", 
                           dims[2], percentVar[dims[2]])) 
    if(circle) {
      frame <- geom_encircle(s_shape = 1, 
                             expand = 0, 
                             alpha = 0.25, 
                             show.legend = F)
      fil <- scale_fill_manual(values = colRamp(cols, 1.5))
    } else {
      frame <- NULL
      fil <- NULL
          }    
    
    plt <- p + frame + text + pts + lab + color + fil + ti + theme
    return(plt)
}

geneCluster <- function(dds, genelist, sortMethod = "groupVar", cluster_rows = T, ntop = NA, main = "", img.path = NA, img.name = NA) {
  
  if (is.na(ntop)) { 
    ntop <- length(genelist)
  } else {
    ntop <- min(length(genelist), ntop)
  }
  gccol <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))  ## color for genecluster
  
  dds <- ddsCombineFactors(dds)
  dds <- dds[match(genelist, row.names(dds)), ]
  row.names(dds) <- ifelse(mcols(dds)$Gene_name=="", 
                           row.names(dds), 
                           paste0(row.names(dds), " (", mcols(dds)$Gene_name, ")"))
  
  ## calculate group variances
  d <- counts(dds, normalized = T)
  novars <- which(rowVars(d)==0)
  d[novars, ] <- d[novars,] + abs(rnorm(ncol(d) * length(novars), 0, 0.00001))
  d.melt <- melt(d)
  d.melt <- merge(d.melt, as.data.frame(colData(dds)), by.x = "Group", by.y = 0, all.x=T)
  d.melt <- select(d.melt, Var1, Group, value)
  d.melt$Group <- as.factor(d.melt$Group)
  m <- tapply(d.melt$value, list(d.melt$Var1, d.melt$Group), mean)
  #m <- m[order(row.names(m)), drop = F]
  
  if (sortMethod == "var") {
     keep <- head(order(rowVars(d), decreasing = T),ntop)
  } else if (sortMethod == "sum") {
     keep <- head(order(rowSums(d), decreasing = T), ntop)
  } else if (sortMethod == "mean") {
     keep <- head(order(rowMeans(d), decreasing = T), ntop)
  } else if (sortMethod == "groupVar") {
     keep <- head(order(rowVars(m), decreasing = T), ntop)
  } else {
     keep <- NA
  }
  
  if (!is.na(keep[1])) {
     mat <- as.matrix(d)[keep, , drop=F]
  } else {
     mat <- as.matrix(d)
  }

  if(nrow(mat) < 2 ) {
    message("Nothing to cluster: Less then 2 genes")
    return()
  }

  #colnames(mat) <- dds$Group
  title <- main
  pheatmap(mat, 
           scale="row", 
           trace="none", 
           cluster_rows = cluster_rows,
           cellwidth = 60,
           cellheight = 12,
           col = gccol, 
           border_color = NA,
           angle_col = 90, 
           fontsize = 12,  
           show_rownames = T,  
           filename = ifelse(is.na(img.name), NA, file.path(img.path, img.name)), 
           #width = img$w, height = img$w,
           main = main)
  
}

plotCounts2 <- function(dds, genelist, intgroup="Genotype", pty = "point", 
                        countthreshold = 10, use.log = T, showLab = F, cols = NULL,
                        ncol = NA, showTrend = T, showMean = T, showErrBar = F, 
                        main="", xlab="", ylab="Normalized Counts") {
  
  if (length(genelist)==0) stop("no gene for ploting")
  
  mincount <- ncol(dds) * countthreshold
  fts <- ddsGetFactors(dds) ## extract factors from designs
  
  dds <- dds[match(genelist, row.names(dds)), ]
  dds <- dds[rowSums(counts(dds, normalized = T)) >= mincount, ]
  if ("Gene_name" %in% colnames(mcols(dds))) {
      row.names(dds) <- ifelse(!mcols(dds)$Gene_name == "", paste0(row.names(dds), " (", mcols(dds)$Gene_name, ")"), row.names(dds))
  }
  
  cnts <- counts2(dds, normalize = T)
  if(use.log) {
     cnts$Counts <- log2(cnts$Counts)
     ylab <- "log2(Expr)"
  }
  colnames(cnts)[1] <- "Gene_name"
  
  if (is.na(ncol)) {
    ncol <- ncolOptimize(nrow(cnts))
  }
  
  ## set x axis factor
  if (length(fts) == 1) {
    msg <- sprintf("only one factor found. use %s as X axis\n", fts)
    message(msg)
    x <- fts
    intgroup <- ifelse(intgroup == fts, NULL, intgroup)
  } else if(length(fts) > 1 & is.null(intgroup)) {
    intgroup <- fts[1]
    x <- fts[2]
  } else if (!intgroup %in% fts) {
    intgroup <- fts[1]
    x <- fts[2]
  } else {
    x <- fts[which(!fts %in% intgroup)]
  }
  
  grid.col <- "."
  grid.row <- ifelse(intgroup == "Gene_name", ".", "Gene_name")
  if (pty == "bar") {
    grid.row <- "."
    grid.col <- "."
    x <- "Gene_name"
    intgroup <- "Treatment"
    showTrend <- F
    showMean <- F
    showErrBar <- T
  }
  
  plt <- gg_Plot(cnts, y="Counts", x = x, group_by = intgroup, 
                 plt_type = pty, showTrend = showTrend, showMean = showMean, 
                 showErrBar = showErrBar, showLab = showLab, cols = cols,
                 grid.row = grid.row, grid.col = grid.col, ftcol = ncol, 
                 graph.title = main, xlab = xlab, ylab = ylab)
  
  return(plt)
}

gg_Plot <- function(df, y, x, group_by = NULL, plt_type = "point", 
                    showTrend = T, showMean = T, showErrBar = F, showLab = F,
                    grid.row = ".", grid.col = ".", ftcol=1, cols = NULL,
                    graph.title = "", xlab="", ylab="") {
  
  n <- ifelse(is.null(group_by), 0.1, nlevels(df[[group_by]]))
  p <- ggplot(df, aes_string(x = x, y = y, group = group_by, col = group_by),na.rm = TRUE) 
  if(plt_type %in% c("p", "point")) {
    geom_type <- geom_point(aes_string(col=group_by),
                            position = position_dodge(0.25 * (n-1)),
                            size=2)
  } else if (plt_type %in% c("b", "bar")) {
    geom_type <- geom_bar(aes_string(fill=group_by, col=group_by),
                          position = position_dodge(0.9),
                          stat= "summary",
                          fun.y = "mean")
  }
  
  if(!is.null(cols)) {
    col.fill <- scale_fill_manual(values = cols)
    col.line <- scale_color_manual(values = cols)
  } else {
    col.fill <- NULL
    col.line <- NULL
  }
  
  if (showErrBar) {
    err <- geom_errorbar(aes_string(group = group_by),
                         position = position_dodge(0.9), # make sure the number is the same as the one in geom_bar()
                         stat = "summary",
                         #fun.y = "mean_se",
                         width = 0.5)
  } else {
    err <- NULL
  }
  
  if(showLab) {
    #txt <- geom_text(aes_string(label = "Label"))              
    txt <- geom_text_repel(aes_string(label = "Label"), 
                    size = 3*12/14,
                    #color=colRamp(labelcol,1.5),
                    hjust="inward", vjust = "inward", 
                    show.legend = F)
  } else {
    txt <- NULL
  }
  if (showTrend) {
    trend <- geom_line(aes_string(group = group_by),
                   col = "gray85",
                   lwd = 0.75,
                   linetype = 1,
                   position = position_dodge(0.25* (n-1)),
                   stat = "summary",
                   fun.y = "mean")
  } else {
    trend <- NULL
  }
  
  if (showMean) {
    m <- geom_point(aes_string(x = x, y= y, group = group_by),
                    #col = "red",
                    pch = 3,
                    #lwd = 0.75,
                    size = 3,
                    position = position_dodge(0.25 * (n-1)),
                    stat = "summary",
                    fun.y = "mean")
  } else {
    m <- NULL
  }
  
  ti <- ggtitle(graph.title)
  
  if(grid.row =="." & grid.col == ".") {
    ft <- NULL
  } else {
    ft <- facet_wrap(reformulate(grid.col, grid.row), ncol=ftcol, scales = 'free_y')
  }
  yx <- scale_y_continuous(name = ylab,  
                           # limits = c(0,30),  #set limit in scale_y will throw away outbound data
                           expand = c(0,0))
  yxlim <- expand_limits(y=0)
  xx <- scale_x_discrete(name = xlab, expand = c(0,0.25)) 
  th <- theme_bw()
  
  myplot <- p + ti + trend + geom_type + m + err + txt + ft + xx + yx + yxlim + th + col.fill + col.line
  return(myplot)
}

ggRadar <- function(df, y, x, group = NULL, col = NA, threshold = 0.01, angle.start = 0, showLabel = F, title = "") {
  
  #if(showThreshold) df$thd <- -log10(threshold)
  
  p <- ggplot(df, aes_string(x=x, y=y, group=group, color=group, fill=group))
  b <- geom_bar(width = 1, stat = "identity", position = position_dodge()) 
  if(!is.na(col)) {
   b.fill <- scale_fill_manual(values = col)
   b.col <- scale_color_manual(values = col)
  } else {
   b.fill <- NULL
   b.col <- NULL
  }
  
  if(is.na(threshold)) {
   pg <- NULL
  } else {
   pg <- geom_polygon(aes_string(x=x, y=rep(-log10(threshold),nrow(df))), 
                      color = colRamp("red",1.2), 
                      lty = 1, 
                      fill = NA, 
                      size = 1.5) 
  }
  
  cp <- coord_polar(start=angle.start) 
  yx <- scale_y_continuous(expand = c(0, 0, 0, -1.2))
  thbw <-  theme_bw()
  th <-   theme(#axis.text = element_text(vjust = 4),
      axis.line.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),         
      axis.ticks = element_blank(),
      axis.line.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      panel.grid.major.y = element_line(colour = c(rep("gray90",6),"gray50"), size=c(rep(1,6),1.5)),
      panel.grid.major.x = element_line(color = c(rep("gray90", 21), "black"), size=0.5 * 2),
      panel.border = element_blank()
    )
  
  if(showLabel){
    # thx <- theme(axis.text.x = element_text(angle = seq(angle.start + 270,angle.start - 90,length.out = 22), size = 12))
    # thx <- theme(axis.text.x = element_text(size = 12))
    th <- NULL
  } else {
    # thx <- theme(axis.text.x = element_blank())
  }
  
  return(p + b + b.fill + b.col + pg + cp + yx + thbw + th)
}

colRamp <- function(col, amp){

      newcol <- col2rgb(col)
      if(amp <=1) {
         newcol <- newcol * amp # you can use of course other values than 2. Higher values the darker the output.
      } else if (amp > 1) {
         newcol <- newcol/amp + (amp-1)*255/amp
      }
      newcol <- rgb(t(newcol), maxColorValue=255)
      return(newcol)

}

ncolOptimize <- function(n) {
  div <- seq_len(n)
  fts <- div[n %% div == 0]
  f <- fts[1:ceiling(length(fts)/2)]
  
  nrow <- f[length(f)]
  ncol <- n/nrow
  
  if(ncol > 6 & ncol/nrow > 2) {
    n <- n+1
    ncol <- ncolOptimize(n)
  }
  
  return(ncol)
}

colorMatch <- function(samples, cols) {
  nrep <- ceiling(length(samples)/length(cols))
  if(nrep > 1) cols <- rep(cols, nrep)
  
  if(!class(samples) == "factor") samples <- factor(samples, levels=unique(samples))
  if(!class(cols) == "factor") cols <- factor(cols, levels=unique(cols))
  
  c.df <- data.frame(col=cols, lvls = as.numeric(cols))
  s.df <- data.frame(sample = samples, lvls = as.numeric(samples))
  
  s.df$col <- levels(cols)[s.df$lvls]

  return(s.df$col)
  
}

paretoh <- function(x, label = NULL, col = NA, rampCol = F,  addline.v = NA, 
                    addline.h = NA, title=NULL, xlabel = NULL, ylabel=NULL, 
                    img.path = NA, img.name = NA) {
  
  # take a named numeric vector and return a barplot in vertical version
  if(any(is.na(names(x)))) names(x)[is.na(names(x))] <- ""
  textlength <- max(nchar(names(x)))
  #op <- par()
  #if(min(op$pin) <=0) par(pin = c(3,3))
  col <- ifelse(is.na(col), colRamp("red3",1.4), col)
  col.step <- x/max(x)
  if(rampCol)  col <- sapply(col.step, function(cs) colRamp(col, cs))
  
  if(!is.na(img.path)) {
      pdf(file.path(img.path, img.name), 
          width = textlength * 0.2, 
          height = length(x) * 0.3 + 1)
  }
  par(mar=c(4,textlength * 0.5,1,3))
  bp <- barplot(x, #names.arg = name, cex.names = 0.7,
                horiz = T,
                width = 0.2, space = 0.2, border = NA, axes = F,
                col = col,
                xlab = xlabel,
                cex.lab = 1.4,
                ylab = "",
                yaxt = "n",
                main = title)
  tt <- text(x = x, y = bp, label = label, pos = 4, 
             cex = 0.75, col = "gray15", xpd = T)## add labels to top of bar
  ax2 <- axis(side = 2, at = bp, labels = names(x), tick = F, las= 2, 
              xpd = F, col.axis = "gray15", col = "gray15", 
              mgp=c(3,0,0), cex.axis=1.2)
  ax1 <- axis(side = 1, labels= T, 
              col.axis = "gray15", col = "gray15", 
              cex.axis = 1.2, las=1, 
              mgp=c(1,1,0))
  if(!is.na(addline.h)) lh <- abline(h = addline.h, lty = 2, col = "gray90", lwd = 1.5)
  if(!is.na(addline.v)) lv <- abline(v = addline.v, lty = 2, col = "gray90", lwd = 1.5)
  yx <- abline(v=0, col = "gray15", lwd = 1.1)
  if (!is.na(img.path)) {
    graphics.off()
  } else {
    return(bp + tt + ax2 + ax1 + lh + lv + yx)
  }
  
  
  #suppressWarnings(par(op))
  
}

colorGenerator <- function(n, type = c("qual", "div", "seq"), theme = NULL) {
  ## to generate a series of colors
  ## mainly using RColorBrewer package
  ## theme only for diverging or continous colors (blue,red etc.)
  
  if(length(type) > 1) type <- type[1]
  if(!type %in% c("qual", "div", "seq")) stop("wrong type")
  if(type == "qual") theme <- NULL
  n <- floor(abs(n))
  
  divColMsg <- 'For \"div\" types, themes should be one of these:
     1. one of "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"; 
     2. a vector of two or three color names to specify the left, middle, and right colors. 
        When only two colors are specified, a gray color will be inserted in the middle. 
        Color names should be,
        1). in hex (#FF9900) format, or
        2). one of the names from grDevices::colors()'
  
  pal.info <- RColorBrewer::brewer.pal.info
  pals <- pal.info[pal.info$category == type, ]
  
  if(type == "qual") {
      if(n < 75) {
          pals<- pals[c(1,2,3,6,7,8,4,5), ]
          cols <- mapply(RColorBrewer::brewer.pal, pals$maxcolors, rownames(pals))
          cols <- unlist(cols)[1:n]
      } else {
          cols = grDevices::colors()[grep('gray|grey|white|yellow', 
                                           grDevices::colors(), 
                                           invert = T)]
          cols <- sample(cols, n)
      }
    yelo <- grep("#FFFF", cols)
    mixee <- sample(cols[-yelo], length(yelo))
    mixed <- colorMixer(cols[yelo], mixee, 0.5)
    cols[yelo] <- mixed
    cols <- colorUniquerize(cols)
    return(cols)
    
  } else if( type == "div") {
      if(is.null(theme)) theme <- "RdYlBu"
      if(length(theme) == 1) {
          if(!theme %in% row.names(pals)) stop(divColMsg)
          col_name <- row.names(pals)[which(pals$category) == theme]
          col_n <- pals$maxcolors[which(pals$category) == theme]
          cols <- colorRampPalette(RColorBrewer::brewer.pal(col_n, col_name))(n)
      } else {
          if(length(theme) == 2) theme <- c(theme[1], "gray90", theme[2])
          cols <- colorRampPalette(theme)(n)
      }
      return(cols)
    
  } else {
      if(is.null(theme)) theme <- "Blues"
      if(theme %in% row.names(pals)) {
        col_n <- pals$maxcolors[which(row.names(pals) == theme)]
        cols <- colorRampPalette(brewer.pal(col_n, theme))(n)
      } else {
        cols <- colorRampPalette(c("gray90", theme))(n)
      }
      return(cols)
  }
  
}

generateColors <- function(n, theme = "bw") {
  
  theme <- tolower(theme)
  if(theme %in% c("r", "red", "rd")) {
    theme <- "r"
  } else if(theme %in% c("blue", "bl", "b")) {
    theme <- "b"
  } else if(theme %in% c("green", "g")) {
    theme <- "g"
  } else if(theme %in% c("yellow", "y", "yl")) {
    theme <- "gr"
  } else if(theme %in% c("cyan", "c")) {
    theme <- "bg"
  } else if(theme %in% c("magenta", "m")) {
    theme <- "br"
  }
  rgb <- lapply(1:3, function(i) floor(seq(10,  225, length.out = n)))
  rgb <- lapply(rgb, as.hexmode)
  rgb <- lapply(rgb, as.character)
  rgb.df <- as.data.frame(rgb)
  colnames(rgb.df) <- c("r", "g", "b")   
  if(theme == "bw") {
    colors <- sprintf("#%s", apply(rgb.df, 1, str_c, collapse = ""))
  } else if (theme %in% c("r", "g", "b")) {
    rgb.df[, !colnames(rgb.df) == theme] <- "00"
    colors <- sprintf("#%s", apply(rgb.df, 1, str_c, collapse = ""))
  } else if (theme %in% c("gr", "br", "bg")) {
    FFcols <- sapply(colnames(rgb.df), function(x) grepl(x, theme))
    rgb.df[, FFcols] <- "FF"
    colors <- sprintf("#%s", apply(rgb.df, 1, str_c, collapse = ""))
  }
  
  return(colors)
  
}

colorUniquerize <- function(colors) {
  if(is(colors, "matrix")) {
    if(is.rgb(colors)) {
      colors.x <- rgbFormating(colors)
      colors.x <- rgb(t(colors.x), maxColorValue = max(colors))

    } else {
      stop("Something don't look like color to me.")
    }
  } else {
    colors.x <- colors
  }
  
  ds <- duplicated(colors.x)
  if(sum(ds) == 0) {
    if(is(colors, "matrix")) return(col2rgb(colors.x))
    return(colors.x)
  }
  
  dds <- duplicated(colors.x[length(colors.x):1])[length(colors.x):1]

  mixee <- colors.x[ds]
  mixer.id <- sample(which(!(ds|dds)), sum(ds))
  mixer <- colors.x[mixer.id]
  mixed <- colorMixer(mixer, mixee, alpha = 0.5)
  colors.x[ds] <- mixed
  colorUniquerize(colors.x)
  
}
## to mix two colors,
## mixer, mixee:
##      vectors of colors (in hex) to mix
##      could also be data.frames, with rgb format, with each column represtates a color
## value,
##      a vector of colors in hex format
colorMixer <- function(mixer, mixee, alpha = 0.5, limits = c(0, 255)) {
  if(is(mixer, "matrix")) {
    if(!is.rgb(mixer, limits)) stop("color must be in RGB format when supplied as matrix")
    mixer.x <- rgbFormating(mixer)
  } else {
    mixer.x <- col2rgb(mixer)
    if(any(is.na(mixer.x))) stop("Something don't look like color to me.")
  }
  
  if(is(mixee, "matrix")){
    if(!is.rgb(mixer, limits)) stop("color must be in RGB format when supplied as matrix")
    mixee.x <- rgbFormating(mixee)
  } else {
    mixee.x <- col2rgb(mixee)
    if(any(is.na(mixee.x))) stop("Something don't look like color to me")
  }
  
  mixed <- mixer.x * 0.5 + mixee.x * (1-0.5)
  if(!is.matrix(mixer)) mixed <- rgb(t(mixed), maxColorValue = limits[2])
  
  return(mixed)
  
}

is.rgb <- function(x, limits = c(0, 255)) {
  if(limits[1] < 0 | limits[2] > 255) stop("limits out of range")
  if(is(x, "data.frame")) x <- as.matrix(x)
  if(is(x, "vector") & length(x) %% 3 == 0) x <- matrix(x, nrow = 3)
  if(is(x, "list")) {
    ns <- sapply(x, length)
    if(all(ns == 3)) x <- matrix(unlist(x), nrow = 3)
    else return(0)
  }
  if(!is(x, "matrix")) return(0)
  x <- as.numeric(x)
  if(any(is.na(x))) return(0)
  if(!any(dim(x)==3)) return(0)
  if(max(x) > limits[2]) return(0)
  if(min(x) < limits[1]) return(0)
  return(1)
  
}

## to turn rgb colors to a format that each row represent a color
rgbFormating <- function(x, limits = c(0, 255)) {
  if(!is.rgb(x, limits = limits)) stop("Don't look like rgb")
  if(is(x, "data.frame")) x <- as.matrix(x)
  if(is(x, "vector")) x <- matrix(x, nrow = 3)
  if(is(x, "list")) x <- matrix(unlist(x), nrow = 3)
  
  x <- as.numeric(x)
  if(ncol(x) == 3) x <- t(x)
  row.names(x) <- c("r", "g", "b")
  return(x)
}

## to display colors
## Arguments:
##  color: a vector of color names or hex color code, or
##         a matrix or data.frame of RGB values
##  byrow: each row represents a color, with 3 values (red, green, blue)
##         if false, each color represents a color. 
colorDisplay <- function(color, byrow = T, showValue = c("hex", "rgb", "none")){
  if(is(color, "data.frame")) {
    color <- as.matrix(color)
  }
  if(is(color, "matrix")) {
    if(max(color) > 255 | min(color) < 0) {
      stop("Don't look like color to me (out of range)")
    }
    if(byrow & !ncol(color) == 3) {
      stop("Don't look like color to me (make sure it has 3 columns)")
    }
    if(!byrow & !nrow(color) == 3) {
      stop("Don't look like color to me")
    }
    if(!byrow) {
      color <- t(color)
    }
    
    color <- apply(color, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    color <- unlist(color)
  }
  
  if(length(showValue) > 1) showValue <- showValue[1]
  
  if(showValue == "rgb") {
    rgb.values <- col2rgb(color)
    values <- apply(rgb.values, 2, function(x) sprintf("(%s)", paste0(x, collapse = ",")))
    values <- unlist(values)
  } else (
    values <- color
  )
  
  n <- length(color)
  image(1:n,1, matrix(1:n), col = color, yaxt = "n", xlab = "", ylab = "")
  if(!showValue == "none") {
    text(1:n, y = rep(0.62, n), labels = values, 
         adj = c(0, 0.5), srt = 90, cex = 0.75)
  }
}

## to display sequential cols from RColorBrewer
colorDisplaySeq <- function(x, n = 255) {
  ncolor <- length(x)
  colmax <- brewer.pal.info[row.names(brewer.pal.info) %in% x, "maxcolors"]
  colmax <- unlist(colmax)
  col_list <- mapply(brewer.pal, colmax, x)
  if(is(col_list, "matrix")){
    col_list <- split(col_list, rep(1:ncolor, each = nrow(col_list)))
  } 
  cols <- lapply(col_list, function(x) colorRampPalette(x)(n))
  ymax <- 0.05 * (ncolor + 1)
  
  xs <- round(((1:n) - mean(1:n))/max(abs((1:n) - mean(1:n))), 2)
  ys <- seq(0.05, ymax, length.out = ncolor + 1)
 
  plot(xs, rep(ys[1], n), ylim = c(0, ymax), col = "white", bty = "n",
       yaxt = "n", xlab = "", ylab = "")
  pts <- sapply(1:ncolor, function(i) {
           points(xs, rep(ys[i], n), pch = 15, cex = 2, col = cols[[i]])})
 # pts
  
}
