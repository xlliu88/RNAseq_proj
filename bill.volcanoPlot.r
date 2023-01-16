source("./scripts/_loadPackages.r")
#source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")
library(RColorBrewer)

p.thd <- 0.05
fc.thd <- 1.0
lfc.thd <- log2(fc.thd)

fn.prefix <- sprintf("[fc%.1fp%s]", fc.thd, formatC(p.thd, format = "g", digits = 0))
### import results
EIS.res <- importResults(bill.resfile.EISd75)
LCM.res <- importResults(bill.resfile.LCM)


### DEG analysis
### MA plot of DEGs #######################################################################################
vplot.path <- "./output/Bill.Analysis/volcanoPlot2"
if(!file.exists(vplot.path)) dir.create(vplot.path)

## Volcano plot of DEGs in LCM samples
for (n in names(LCM.res)) {
  
  name <- compare.names[[n]]
  fn <- sprintf("Volcano_LCM_%s.%s.pdf", name, fn.prefix)
  df <- LCM.res[[n]]
  col <- ifelse(df$padj < p.thd & abs(df$log2FoldChange) > lfc.thd, "red", "gray")
  col <- ifelse(df$log2FoldChange < 0 & col == "red", "blue", col)
  xmax <- max(abs(df$log2FoldChange))
  
  nup <- sum(col == "red")
  ndn <- sum(col == "blue")
  cat(sprintf("%s:\tup:%d\tdown:%d\t\n", name, nup, ndn))
  pdf(file.path(vplot.path, fn), width = 4, height = 4.5)
  plot(df$log2FoldChange, -log10(df$padj), col = col, xlim = c(-xmax, xmax),
       pch = 20, cex = 0.5, cex.lab = 1, cex.axis = 1,
       main = name,
       xlab = "log2FoldChange",
       ylab = "-log10(Pvalue)")
  text(x = -xmax, y = max(-log10(df$padj) * 0.95), labels = sprintf("down-regulated\n%d", ndn),
       pos = 4, cex = 0.7)
  text(x = xmax, y = max(-log10(df$padj) * 0.95), labels = sprintf("up-regulated\n%d", nup),
       pos = 2, cex = 0.7)
  dev.off()
  
  
}

## volcano plot of DEGs in EIS samples

for (n in names(EIS.res)) {
  name <- compare.names[[n]]
  fn <- sprintf("Volcano_EIS_%s.%s.pdf", name, fn.prefix)
  df <- EIS.res[[n]]
  col <- ifelse(df$padj < p.thd & abs(df$log2FoldChange) > lfc.thd, "red", "gray")
  col <- ifelse(df$log2FoldChange < 0 & col == "red", "blue", col)
  xmax <- max(abs(df$log2FoldChange))
  nup <- sum(col == "red")
  ndn <- sum(col == "blue")
  cat(sprintf("%s:\tup:%d\tdown:%d\t\n", name, nup, ndn))
  # 
  pdf(file.path(vplot.path, fn), width = 4, height = 4.5)
  plot(df$log2FoldChange, -log10(df$padj), col = col, xlim = c(-xmax, xmax),
       pch = 20, cex = 0.5,
       main = name,
       xlab = "log2FoldChange",
       ylab = "-log10(Pvalue)")
  text(x = -xmax, y = max(-log10(df$padj)*0.95), labels = sprintf("down-regulated:\n%d", ndn),
       pos = 4, cex = 0.7)
  text(x = xmax, y = max(-log10(df$padj)*0.95), labels = sprintf("up-regulated:\n%d", nup),
       pos = 2, cex = 0.7)
  dev.off()

}

## MA plot of DEGs overlapped in LCM and EWR samples
for (n in names(LCM.res.bill)) {
  name <- compare.names[[n]]
  df <- LCM.res.bill[[n]]
  df2 <- EWR.res.bill[[n]]
  
  sel <- intersect(df$id, df2$id)
  df <- df[df$id %in% sel, ]
  df2 <- df2[df2$id %in% sel, ]
  
  df <- df[order(row.names(df)),]
  df2 <- df2[order(row.names(df2)),]
  if (!all(row.names(df) == row.names(df2))) {
    message(sprintf("Unmatch Row names for: %s", name))
    next
  }
  
  df$Sig <- df$padj < pthreshold & df2$padj < pthreshold
  
  df <- select(df, baseMean, log2FoldChange, Sig)
  nup <- sum(with(df, log2FoldChange > lfc.thd & Sig))
  ndn <- sum(with(df, log2FoldChange < lfc.thd & Sig))
  
  png(file.path(maplot.path, paste0("LCM&EWR.", name, ".p.e1-03.png")), width = 8, height = 8, unit = "in", res = 300)
  geneplotter::plotMA(df, ylim = c(-5,5), cex=0.75, 
                      colSig = colRamp("red3", 1.4), 
                      xlab = "Mean Expression", 
                      ylab = "log2 Fold Change",
                      cex.lab = 1.5,
                      cex.axis = 1.3)
  text(max(df$baseMean), 5, pos = 2, labels = paste0("Up-regulated: ", nup), cex = 1.3)
  text(max(df$baseMean), -5, pos = 2, labels = paste0("Down-regulated: ", ndn), cex = 1.3)
  title(name, line = 0.5, adj = 0, cex.main = 2)
  dev.off()
}
