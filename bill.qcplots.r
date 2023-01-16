
source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")

#input files
in.path <- "./output/Bill.Analysis/"
qc.path <- "./output/Bill.Analysis/qcplots"


gff3_file <- "./genomes/gffa3/Arabidopsis_thaliana.TAIR10.59.gff3"
genedes_file <- "./genomes/gffa3/Arabidopsis_thaliana.TAIR10.59.des.gaf"
maplots_output <- paste0(qc.path, "pairwiseMAs.pdf")

cols <- brewer.pal(8, "Paired")
cols_EIS <- rep(c(cols[c(2,8,6)], cols[c(1,7,5)]), 3)
cols_EISd75 <- cols_EIS[-5]
cols_LCM <- cols_EIS[-(seq(1:6)*3-1)]
cols_CLE <- cols_EIS[-(seq(1:6)*3)]
## output files
NormmalizedCounts <- "counts_EIS_Normalized.txt"

## read gene description file
geneDes <- extractGeneDesFromGFF3(gff3_file)

## build dds
design = ~ Genotype + Treatment + Genotype:Treatment
dds.EIS <- ddsFromFeatureCounts(file.path(bill.counts$path, bill.counts$EIScnts),
                                file.path(bill.counts$path, bill.counts$EISconditions),
                                source = "tximport", gffa3=gff3_file, design=design,
                                normalize = T)
dds.LCM <- ddsFromFeatureCounts(file.path(bill.counts$path, bill.counts$LCMcnts), 
                                file.path(bill.counts$path, bill.counts$LCMconditions), 
                                source = "tximport", gffa3=gff3_file, design=design,
                                normalize = T)

keep <- rowSums(counts(dds.EIS)) >= 10
dds.EIS <- dds.EIS[keep, ]
dds.EIS$Group <- paste(as.character(dds.EIS$Genotype), 
                      as.character(dds.EIS$Treatment),
                      sep = "_")
dds.EIS$Group <- factor(dds.EIS$Group, levels = unique(dds.EIS$Group))

keep <- rowSums(counts(dds.LCM)) >= 10
dds.LCM <- dds.LCM[keep, ]
dds.LCM$Group <- paste(as.character(dds.LCM$Genotype), 
                       as.character(dds.LCM$Treatment),
                       sep = "_")
dds.LCM$Group <- factor(dds.LCM$Group, levels = unique(dds.LCM$Group))


dds.EISd75 <- ddsRemoveSample(dds.EIS, "clv.CLE.1")
dds.CLE <- ddsSlice(dds.EIS, factorName = "Treatment", factorLevel = c("Ctrl", "HsCLE2"))

# rlogTransformation
rld.EIS <- rlogTransformation(dds.EIS, blind = T)
rld.EISd75 <- rlogTransformation(dds.EISd75, blind = T)
rld.LCM <- rlogTransformation(dds.LCM, blind = T)
rld.CLE <- rlogTransformation(dds.CLE, blind = T)

# DESeq, default
dds2.EIS <-DESeq(dds.EIS)
dds2.LCM <- DESeq(dds.LCM)
dds2.EISd75 <- DESeq(dds.EISd75)

dsets <- c("EIS", "EISd75", "LCM")

## boxplot of cooks distance
for(ds in dsets) {
  fn <- file.path(qc.path, sprintf("%s_cooks_dist.pdf", ds))
  dat <- get(sprintf("dds2.%s", ds))
  pdf(fn, width = 8, height = 4.5)
  boxplot(log10(assays(dat)[["cooks"]]), las=2)
  dev.off()
}

## plot of dispersion estimation
for(ds in dsets) {
  fn <- file.path(qc.path, sprintf("%s_DispEsts.pdf", ds))
  dat <- get(sprintf("dds2.%s", ds))
  pdf(fn, width = 6, height = 6)
  plotDispEsts(dat)
  dev.off()
}

## qc plots
### density and ecdf
for (ds in dsets) {
    dds <- get(sprintf("dds.%s", ds))
    scols <- get(sprintf("cols_%s", ds))
    if(ds == "LCM") {
      lcols <- cols_LCM[1:4]
    } else {
      lcols <- cols_EIS[1:6]
    }
    colData(dds)$Group <- paste(as.character(dds$Genotype), 
                                as.character(dds$Treatment),
                                sep = "_")
    legd <- rev(unique(rev(dds$Group)))
    
    #boxplot
    ti <- sprintf("%s_boxplot", ds)
    fn <- sprintf("%s_boxplot.pdf", ds)
    pdf(file.path(qc.path, fn), width = 8, height = 6.5)
    boxplot(assay(get(sprintf("rld.%s",ds))), 
                 col = scols, 
                 main = ti,
                 las = 2)
    dev.off()   
    
    # density plot
    ti <- sprintf("%s_count_density", ds)
    fn <- sprintf("%s_density.pdf", ds)
    pdf(file.path(qc.path, fn), width = 8, height = 6.5)
    multidensity(counts(dds, normalized = T), 
                 xlim = c(0,1000), 
                 col = scols, 
                 main = ti,
                 xlab = "Normalized Counts",
                 legend = list(
                   x = "topright",
                   legend = legd,
                   pch = 15,
                   col = lcols)
                 )
    dev.off()   

    # ecdf plot
    ti <- sprintf("%s_ecdf", ds)
    fn <- sprintf("%s_ecdf2k.pdf", ds)
    pdf(file.path(qc.path, fn), width = 8, height = 6.5)
    multiecdf(counts(dds, normalized = T), xlim = c(0, 2000), col = scols, 
                 main = ti,
                 xlab = "Normalized Counts",
                 legend = list(
                   x = "bottomright",
                   legend = legd,
                   pch = 15,
                   col = lcols)
                 )
    dev.off()
    }

## PCA plots & hierarchical clustering
for(ds in c(dsets, "CLE")) {
  dds <- get(sprintf("rld.%s", ds))
  
  ## source of variances
  
  for(ntop in c(500, 1000, 5000, nrow(dds))){
    n <- ifelse(ntop == nrow(dds), "all", sprintf("top%d", ntop))
    scols <- get(sprintf("cols_%s", ds))
    
    ## source of variances
    Pvar <- rowVars(assay(dds))
    keep <- order(Pvar, decreasing = T)[seq_len(ntop)]
    pca <- prcomp(t(assay(dds)[keep, ]), scale = F)
    perVar <- round(100*pca$sdev^2/sum(pca$sdev^2), 1)
    names(perVar) <- sprintf("PC%d", seq_len(length(perVar)))
    fn <- sprintf("%s_percent_variance_explained(%s).pdf", ds, n)
    pdf(file.path(qc.path, fn), width = 4, height = 3.5)
    plot(perVar, pch = 16, xaxt = "n", 
         xlab = "PCs", 
         ylab = "% Variance Explained",
         ylim = c(0,100))
    title(main = "% Variance Explained by Each PC",
          line = 0.5,
          adj = 0)
    axis(1, at = seq_len(length(perVar)), 
         labels = sprintf("PC%d", seq_len(length(perVar))), 
         las = 2)
    text(x = seq_len(length(perVar)), y = perVar + 2.5, cex = 0.7,
         labels = sprintf("%.1f%%", perVar), adj = c(0.5,0))
    dev.off()
    next
    
    ## PCA plots
    plt <- PCAplot(assay(dds), dds$Group, ntop = ntop, cols = scols, labsize = 20)
    fn <- sprintf("%s_PCA(%s)_var_genes_print.pdf",ds, n)
    ggsave(file.path(qc.path, fn), width = 5, height = 3.5)
    
    # hierarchical clustering
    dat <- assay(dds)
    topvar <- order(-rowVars(dat))[1:ntop]
    d <- dist(t(dat[topvar, ]))
    hc <- hclust(d)
    ti <- sprintf("Hierarchical Clustering (%s)", n)
    fn <- sprintf("%s_hclu(%s).pdf", ds, n)
    pdf(file.path(qc.path, fn), width = 6, height = 4)
    myplclust(hc, main = ti)
    dev.off()
  }
  
}

#########################################################
## explor of clv CLE2 treated. rep1 
## 
out.path <- file.path(qc.path, "clvCLE_rep1_outlier")
if(!file.exists(out.path)) dir.create(out.path)
fc.thd <- 2
lfc.thd <- log2(fc.thd)
p.thd <- 0.001

EIS.res <- importResults(bill.resfile.EISd75)

## genes up/down regulated after CLE treatment in clv mutant
gene_info <- mcols(dds.EIS)

clvCLE.up.idx <- EIS.res$vec$log2FoldChange > lfc.thd & EIS.res$vec$padj < p.thd
sum(clvCLE.up.idx)
clvCLE.dn.idx <- EIS.res$vec$log2FoldChange < -lfc.thd & EIS.res$vec$padj < p.thd
sum(clvCLE.dn.idx)

clvCLE.up.df <- EIS.res$vec[clvCLE.up.idx, ]
clvCLE.up.df <- clvCLE.up.df[order(-clvCLE.up.df$log2FoldChange), ]

clvCLE.dn.df <- EIS.res$vec[clvCLE.dn.idx, ]
clvCLE.dn.df <- clvCLE.dn.df[order(clvCLE.dn.df$log2FoldChange),]


plotCounts2(dds.CLE, clvCLE.up.df$Gene_id[1:24], 
            intgroup = "Treatment", use.log = F,
            showLab = T, showTrend = F, ncol = 6)
fn <- sprintf("Expression_of_clvCLE.up_genes_in_EIS7-5(1).pdf")
ggsave(file.path(out.path, fn), width = 12, height = 10)
plt <- plotCounts2(dds.CLE, clvCLE.up.df$Gene_id[25:48], 
            intgroup = "Treatment", use.log = F,
            showLab = T, showTrend = F, ncol = 6)
fn <- sprintf("Expression_of_clvCLE.up_genes_in_EIS7-5(2.1).pdf")
ggsave(file.path(out.path, fn), width = 12, height = 10)

plotCounts2(dds.CLE, clvCLE.dn.df$Gene_id[1:24], 
            intgroup = "Treatment", use.log = F,
            showLab = T, showTrend = F, ncol = 6)
fn <- sprintf("Expression_of_clvCLE.dn_genes_in_EIS7-5(1).pdf")
ggsave(file.path(out.path, fn), width = 12, height = 10)
plotCounts2(dds.CLE, clvCLE.dn.df$Gene_id[25:48], 
            intgroup = "Treatment", use.log = F,
            showLab = T, showTrend = F, ncol = 6)
fn <- sprintf("Expression_of_clvCLE.dn_genes_in_EIS7-5(2).pdf")
ggsave(file.path(out.path, fn), width = 12, height = 10)
plotCounts2(dds.CLE, clvCLE.dn.df$Gene_id[49:nrow(clvCLE.dn.df)], 
            intgroup = "Treatment", use.log = F,
            showLab = T, showTrend = F, ncol = 6)
fn <- sprintf("Expression_of_clvCLE.dn_genes_in_EIS7-5(3).pdf")
ggsave(file.path(out.path, fn), width = 12, height = 10)

### marker genes
mkgenes <- sapply(c("CLV1", "CLV2", "RPK2", "WOX4"), #"ATHB-8", "^ANT$"), 
                  function(x) row.names(gene_info)[grep(x, gene_info$Gene_name)])
mkgenes2 <- sapply(c("CLV1", "CLV2", "RPK2", "WOX4", "ATHB-8", "^ANT$"), 
                  function(x) row.names(gene_info)[grep(x, gene_info$Gene_name)])

##plots
plotCounts2(dds.CLE, mkgenes, intgroup = "Treatment", use.log = F,
            showLab = T, showTrend = F)
fn <- sprintf("Expression_of_markder_genes_in_EIS7-5(1).pdf")
ggsave(file.path(out.path, fn), width = 10, height = 2.5)

plotCounts2(dds.CLE, mkgenes2, intgroup = "Treatment", use.log = F,
            showLab = T, showTrend = F, ncol = 4)
fn <- sprintf("Expression_of_markder_genes_in_EIS7-5(2).pdf")
ggsave(file.path(out.path, fn), width = 10, height = 8)


#########################################################
### protein coding genes only - EIS & LCM
keep <- which(mcols(dds.EIS)$Gene_type == "protein_coding")
dds.EISp <- dds.EIS[keep,]
dds.EISd75p <- dds.EISd75[keep,]
dds.CLEp <- dds.CLE[keep, ]

keep <- which(mcols(dds.LCM)$Gene_type == "protein_coding")
dds.LCMp <- dds.LCM[keep, ]

rld.EISp <- rlog(dds.EISp, blind = T)
rld.EISd75p <- rlog(dds.EISd75p, blind = T)
rld.CLEp <- rlog(dds.CLEp, blind = T)
rld.LCMp <- rlog(dds.LCMp, blind = T)

dds2.EISp <- DESeq(dds.EISp)
dds2.EISd75p <- DESeq(dds.EISd75p)
dds2.CLEp <- DESeq(dds.CLEp)
dds2.LCMp <- DESeq(dds.LCMp)

######################
qc.path.p <- file.path(qc.path, "protein_coding")
if(!exists(qc.path.p)) dir.create(qc.path.p)

## qc plots
### density and ecdf
for (ds in dsets) {
    dds <- get(sprintf("dds.%sp", ds))
    scols <- get(sprintf("cols_%s", ds))
    if(ds == "LCM") {
      lcols <- cols_LCM[1:4]
    } else {
      lcols <- cols_EIS[1:6]
    }
    # colData(dds)$Group <- paste(as.character(dds$Genotype), 
    #                             as.character(dds$Treatment),
    #                             sep = "_")
    legd <- levels(dds$Group)
    
    #boxplot
    ti <- sprintf("%s_boxplot", ds)
    fn <- sprintf("%s_boxplot.pdf", ds)
    pdf(file.path(qc.path.p, fn), width = 8, height = 6.5)
    boxplot(assay(get(sprintf("rld.%sp",ds))), 
                 col = scols, 
                 main = ti,
                 las = 2)
    dev.off()   
    
    # density plot
    ti <- sprintf("%s_count_density", ds)
    fn <- sprintf("%s_density.pdf", ds)
    pdf(file.path(qc.path.p, fn), width = 8, height = 6.5)
    multidensity(counts(dds, normalized = T), 
                 xlim = c(0,1000), 
                 col = scols, 
                 main = ti,
                 xlab = "Normalized Counts",
                 legend = list(
                   x = "topright",
                   legend = legd,
                   pch = 15,
                   col = lcols)
                 )
    dev.off()   

    # ecdf plot
    ti <- sprintf("%s_ecdf", ds)
    fn <- sprintf("%s_ecdf2k.pdf", ds)
    pdf(file.path(qc.path.p, fn), width = 8, height = 6.5)
    multiecdf(counts(dds, normalized = T), xlim = c(0, 2000), col = scols, 
                 main = ti,
                 xlab = "Normalized Counts",
                 legend = list(
                   x = "bottomright",
                   legend = legd,
                   pch = 15,
                   col = lcols)
                 )
    dev.off()
    }

## PCA plots & hierarchical clustering
for(ds in c(dsets, "CLE")) {
  dds <- get(sprintf("rld.%sp", ds))
  
  for(ntop in c(500, 1000, 5000, nrow(dds))){
    ## PCA plots
    scols <- get(sprintf("cols_%s", ds))
    plt <- PCAplot(assay(dds), dds$Group, ntop = ntop, cols = scols, labsize = 20)
    n <- ifelse(ntop == nrow(dds), "all", sprintf("top%d", ntop))
    fn <- sprintf("%s_PCA(%s)_var_genes.pdf",ds, n)
    ggsave(file.path(qc.path.p, fn), width = 8, height = 6)
    
    # hierarchical clustering
    dat <- assay(dds)
    topvar <- order(-rowVars(dat))[1:ntop]
    d <- dist(t(dat[topvar, ]))
    hc <- hclust(d)
    ti <- sprintf("Hierarchical Clustering (%s)", n)
    fn <- sprintf("%s_hclu(%s).pdf", ds, n)
    pdf(file.path(qc.path.p, fn), width = 6, height = 4)
    myplclust(hc, main = ti)
    dev.off()
  }
  
}





