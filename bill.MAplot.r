source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")
library(RColorBrewer)

pthreshold <- 0.05
fc.thd <- 1.0
lfc.thd <- log2(fc.thd)

### import results
EIS.res.bill <- importResults(bill.resfile.EISd75, "wbc")
LCM.res.bill <- importResults(bill.resfile.LCM, "wbc")
LCM.BCN <- consolidateData(data = LCM.res.bill, "BCN")
EIS.BCN <- consolidateData(EIS.res.bill, "BCN")
EIS.CLE <- consolidateData(EIS.res.bill, "CLE")

### DEG analysis
### MA plot of DEGs #######################################################################################
maplot.path <- "./output/Bill.Analysis/maplots"

## MA plot of DEGs in LCM samples
for (n in names(LCM.res.bill)) {
  
  name <- compare.names[[n]]
  df <- LCM.res.bill[[n]]
  df$Sig <- df$padj < pthreshold
  df <- select(df, baseMean, log2FoldChange, Sig)
  nup <- sum(with(df, log2FoldChange > lfc.thd & Sig))
  ndn <- sum(with(df, log2FoldChange < lfc.thd & Sig))
  
  png(file.path(maplot.path, paste0("LCM.", name, ".p.05.png")), width = 8, height = 8, unit = "in", res = 300)
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

## MA plot of DEGs in EWR samples
for (n in names(EIS.res.bill)) {
  
  name <- compare.names[[n]]
  df <- EIS.res.bill[[n]]
  df$Sig <- df$padj < pthreshold
  df <- select(df, baseMean, log2FoldChange, Sig)
  nup <- sum(with(df, log2FoldChange > lfc.thd & Sig))
  ndn <- sum(with(df, log2FoldChange < lfc.thd & Sig))
  
  png(file.path(maplot.path, paste0("EIS.", name, ".p.05.pdf")), 
      width = 8, height = 8, unit = "in", res = 300)
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

## MA plot of DEGs overlapped in LCM and EWR samples
for (n in names(LCM.res.bill)) {
  name <- compare.names[[n]]
  df <- LCM.res.bill[[n]]
  df2 <- EIS.res.bill[[n]]
  
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
  
  pdf(file.path(maplot.path, paste0("LCM&EIS.", name, ".p.05.pdf")), 
      width = 8, height = 8, unit = "in", res = 300)
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

## venn plot EIS vs LCM ##########
euler.path <- file.path("./output/Bill.Analysis/EIS.vs.LCM/euler_plot")
eis.up <- EIS.res.bill$wbc %>% filter(log2FoldChange > 0, padj < pthreshold) %>% 
  select(Gene_id) %>% unlist() %>% unname()
eis.dn <- EIS.res.bill$wbc %>% filter(log2FoldChange < 0, padj < pthreshold) %>% 
  select(Gene_id) %>% unlist() %>% unname()

lcm.up <- LCM.res.bill$wbc %>% filter(log2FoldChange > 0, padj < pthreshold) %>% 
  select(Gene_id) %>% unlist() %>% unname()
lcm.dn <- LCM.res.bill$wbc %>% filter(log2FoldChange < 0, padj < pthreshold) %>% 
  select(Gene_id) %>% unlist() %>% unname()

list.up <- list(EIS = eis.up,
                LCM = lcm.up)
list.dn <- list(EIS = eis.dn,
                LCM = lcm.dn)


fn <- file.path(euler.path, "WT-")
pdf(fn, width = 3, height = 3, unit = "in")
plot(euler(vlist), 
     fill = colRamp(brewer.pal(3, "Set2"),1.5),
     edges = colRamp(brewer.pal(3, "Set2"),1),
     quantities = T, 
     labels = T )#identical(legend, T))
# dev.off()
graphics.off()





### log2FoldChange difference in WT and clv triple mutant ##################################################
### LCM WT vs clv triple upon BCN treatment ================================================================
df <- LCM.BCN[,1:10]
df$color <- ifelse(df$p.wbc < pthreshold, "lightblue", colRamp("gray60",1)) # sigs in wbc
df$color <- ifelse(df$p.vbc < pthreshold, colRamp("palegreen",1.3), df$color)
df$color <- ifelse(df$p.wbc < pthreshold & df$p.vbc < pthreshold, colRamp("orange",2), df$color)
df$color <- ifelse(df$p.clvBCN < pthreshold, "red", df$color)
df <- df[order(df$color),]

lim <- max(abs(df$lfc.wbc), abs(df$lfc.vbc))
axs.lim <- c(-10,10) #c(-lim, lim)
png(file.path(maplot.path, paste0("LCM.lfc.wbc.vs.vbc.lim10.png")), width = 8, height = 8, unit = "in", res = 300)
plot(df$lfc.wbc, df$lfc.vbc, type = "p", 
     xlab = "log2 Fold Change in WT", 
     ylab = "log2 Fold Change in mutant", 
     main = "log2 Fold Change in WT and Mutant Upon BCN Infection (LCM)",
     pch = 19, 
     cex = 0.7,
     col = df$color,
     xlim = axs.lim, ylim = axs.lim)
abline(0,1, col = "red", lwd = 2)
abline(h=0, col = "palegreen", lty = 2, lwd = 2)
abline(v=0, col = colRamp("blue",1.6), lty = 2,lwd = 2)
dev.off()

LCM.bcn.extreme <- filter(LCM.BCN, abs(lfc.wbc) > 10 | abs(lfc.vbc) > 10)
LCM.bcn.extreme.nonsig <- filter(LCM.bcn.extreme, p.clvBCN >= 0.05)
plotCounts2(dds.LCM, LCM.bcn.extreme$id)
ggsave(file.path(maplot.path, "LCM.BCN.extremelfc.genes.png"), width = 10, height = 3, unit = "in", dpi = 300)


### EWR WT vs clv triple upon BCN treatment ================================================================
df <- EWR.BCN[,1:10]
df$color <- ifelse(df$p.wbc < pthreshold, "lightblue", colRamp("gray60",1)) # sigs in wbc
df$color <- ifelse(df$p.vbc < pthreshold, colRamp("palegreen",1.3), df$color)
df$color <- ifelse(df$p.wbc < pthreshold & df$p.vbc < pthreshold, colRamp("orange",2), df$color)
df$color <- ifelse(df$p.clvBCN < pthreshold, "red", df$color)
df <- df[order(df$color),]

lim <- max(abs(df$lfc.wbc), abs(df$lfc.vbc))
axs.lim <- c(-10,10) #c(-lim, lim)
png(file.path(maplot.path, paste0("EWR.lfc.wbc.vs.vbc.lim10.png")), width = 8, height = 8, unit = "in", res = 300)
plot(df$lfc.wbc, df$lfc.vbc, type = "p", 
     xlab = "log2 Fold Change in WT", 
     ylab = "log2 Fold Change in mutant", 
     main = "log2 Fold Change in WT and Mutant Upon BCN Infection (EWR)",
     xlim = axs.lim, 
     ylim = axs.lim,
     pch = 19, 
     cex = 0.7,
     col = df$color)
abline(0,1, col = "red", lwd = 2)
abline(h=0, col = "palegreen", lty = 2, lwd = 2)
abline(v=0, col = colRamp("blue",1.6), lty = 2,lwd = 2)
dev.off()

EWR.bcn.extreme <- filter(EWR.BCN, abs(lfc.wbc) > 10 | abs(lfc.vbc) > 10)
EWR.bcn.extreme.nonsig <- filter(EWR.bcn.extreme, p.clvBCN >= 0.05)
plotCounts2(dds.BCN, EWR.bcn.extreme$id)
ggsave(file.path(maplot.path, "EWR.BCN.extremelfc.genes.png"), width = 10, height = 6, unit = "in", dpi = 300)


### EWR WT vs clv triple upon CLE treatment ================================================================
df <- EWR.CLE[,1:10]
df$color <- ifelse(df$p.wec < 0.05, "lightblue", colRamp("gray60",1)) # sigs in wbc
df$color <- ifelse(df$p.vec < 0.05, colRamp("palegreen",1.3), df$color)
df$color <- ifelse(df$p.wec <0.05 & df$p.vec < 0.05, "orange", df$color)
df$color <- ifelse(df$p.clvCLE < 0.05, "red", df$color)
df <- df[order(df$color),]

lim <- max(abs(df$lfc.wec), abs(df$lfc.vec))
axs.lim <- c(-lim, lim)
axs.lim <- c(0,2)
png(file.path(maplot.path, paste0("EWR.lfc.wec.vs.vec.6x6.wox4.png")), width = 6, height = 6, unit = "in", res = 300)
plot(df$lfc.wec, df$lfc.vec, type = "p", 
     xlab = "log2 Fold Change in WT", 
     ylab = "log2 Fold Change in mutant", 
     main = "log2 Fold Change in WT and Mutant \nUpon HsCLE2 Treatment (EWR)",
     pch = 19, 
     cex = 1,
     col = df$color,
     xlim = axs.lim, ylim = axs.lim)
abline(0,1, col = "red", lwd = 2)
abline(h=0, col = "palegreen", lty = 3, lwd = 1.5)
abline(v=0, col = "blue", lty = 3,lwd = 1.5)
abline(h=0.4454, col = "gray75")
abline(v = 1.6388, col = "gray75")
text(x=1.6388, y=0.4454, labels = "WOX4:p=0.4452", col = "black", cex = 0.8, pos = 4)
dev.off()

CLE.extreme <- filter(EWR.CLE, abs(lfc.wec) > 10 | abs(lfc.vec) > 10)
plotCounts2(dds.CLE, CLE.extreme$id)
ggsave(file.path(maplot.path, "EWR.CLE.extremelfc.genes.png"), width = 10, height = 3, unit = "in", dpi = 300)


### EWR vs LCM for wbc samples ############
wbc.lcm <- LCM.BCN[, grepl("wbc", colnames(LCM.BCN))]
wbc.ewr <- EWR.BCN[,grepl("wbc", colnames(EWR.BCN))]
colnames(wbc.lcm) <- paste0("LCM.", colnames(wbc.lcm))
colnames(wbc.ewr) <- paste0("EWR.", colnames(wbc.ewr))
wbc <- merge(wbc.lcm, wbc.ewr, by=0)
row.names(wbc) <- wbc$Row.names
wbc <- select(wbc, -Row.names)

wbc$color <- ifelse(wbc$LCM.p.wbc < pthreshold, "lightblue3", colRamp("gray60",1)) # sigs in wbc
wbc$color <- ifelse(wbc$EWR.p.wbc < pthreshold, colRamp("orange",1.3), wbc$color)
wbc$color <- ifelse(wbc$LCM.p.wbc < pthreshold & wbc$EWR.p.wbc < pthreshold, "red", wbc$color)
#wbc$color <- ifelse(wbc$p.clvBCN < pthreshold, "red", wbc$color)
wbc <- wbc[order(wbc$color),]
#wbc$color <- as.factor(wbc$color)

lim <- max(abs(wbc$LCM.lfc.wbc), abs(wbc$EWR.lfc.wbc))
axs.lim <- c(-lim, lim)
#axs.lim <- c(-10,10) 
png(file.path(maplot.path, paste0("LCM.lfc.wbc.vs.vbc.png")), width = 8, height = 8, unit = "in", res = 300)
plot(wbc$LCM.lfc.wbc, wbc$EWR.lfc.wbc, type = "p", 
     xlab = "log2 Fold Change in LCM Samples", 
     ylab = "log2 Fold Change in EWR Samples", 
     main = "log2 Fold Change in LCM and EWR Samples Upon BCN Infection",
     pch = 19, 
     cex = 0.7,
     col = wbc$color,
     xlim = axs.lim, ylim = axs.lim)
abline(0,1, col = "red", lwd = 2)
abline(h=0, col = "orange", lty = 2, lwd = 2)
abline(v=0, col = colRamp("blue",1.6), lty = 2,lwd = 2)
dev.off()

wbc$id <- row.names(wbc)
LupEdn <- filter(wbc, LCM.p.wbc < 0.05, EWR.p.wbc < 0.05, LCM.lfc.wbc > 0, EWR.lfc.wbc < 0)


### union of clvBCN in LCM and EWR ========================================================================

clvBCN.EWR.sig <- row.names(EWR.BCN)[EWR.BCN$p.clvBCN < 0.05]
clvBCN.LCM.sig <- row.names(LCM.BCN)[LCM.BCN$p.clvBCN < 0.05]

u <- union(clvBCN.EWR.sig, clvBCN.LCM.sig)
its <- intersect(clvBCN.EWR.sig, clvBCN.LCM.sig)

df1 <- EWR.BCN[row.names(EWR.BCN) %in% u,]
df2 <- LCM.BCN[row.names(LCM.BCN) %in% u,]

df1 <- df1[, c("lfc.clvBCN", "p.clvBCN")]
df2 <- df2[, c("lfc.clvBCN", "p.clvBCN")]
colnames(df1) <- paste0("EWR.", c("lfc", "p"))
colnames(df2) <- paste0("LCM.", c("lfc", "p"))
df <- simpleMerge(df1,df2,by=0, all=T)
df$EWR.p[is.na(df$EWR.p)] <- 1
df$LCM.p[is.na(df$LCM.p)] <- 1
df$EWR.lfc[is.na(df$EWR.lfc)] <- 0
df$LCM.lfc[is.na(df$LCM.lfc)] <- 0


### plot counts of interested genes ===============================================
cntplt.path <- "./output/Bill.Analysis/countplots"
df.LnEsig <- df[df$LCM.p < 0.05 & df$EWR.p < 0.05,]
df.pos <- df.LnEsig[with(df.LnEsig, LCM.lfc > 0 & EWR.lfc > 0),]
df.neg <- df.LnEsig[with(df.LnEsig, LCM.lfc < 0 & EWR.lfc < 0),]

plotCounts2(dds.EWRd75, row.names(df.pos))
ggsave(file.path(cntplt.path, "LnE.clvBCN.pos.EWRexp.png"), width = 8, height = 8, dpi = 300, units = "in")
plotCounts2(dds.LCM, row.names(df.pos))
ggsave(file.path(cntplt.path, "LnE.clvBCN.pos.LCMexp.png"), width = 8, height = 8, dpi = 300, units = "in")

plotCounts2(dds.EWRd75, row.names(df.neg))
ggsave(file.path(cntplt.path, "LnE.clvBCN.neg.EWRexp.png"), width = 8, height = 8, dpi = 300, units = "in")
plotCounts2(dds.LCM, row.names(df.neg))
ggsave(file.path(cntplt.path, "LnE.clvBCN.neg.LCMexp.png"), width = 8, height = 8, dpi = 300, units = "in")

### list of These genes
clvBCN.pos.genes <- geneInfo[row.names(geneInfo) %in% row.names(df.pos),]
clvBCN.neg.genes <- geneInfo[row.names(geneInfo) %in% row.names(df.neg),]

write.table(clvBCN.pos.genes, file.path(cntplt.path, "GeneList.clvBCN.pos.csv"), row.names = F, col.names=T, sep = ",", quote = F)
write.table(clvBCN.neg.genes, file.path(cntplt.path, "GeneList.clvBCN.neg.csv"), row.names = F, col.names=T, sep = ",", quote = F)










############################
#Example Plot to Explain the Comparisons
n <- 5
se <- 150
a <- rnorm(n, 1000, se)
b <- rnorm(n, 5000, se)
c <- rnorm(n, 200, se)
d <- rnorm(n, 2000, se)

df <- data.frame(a,b,c,d)
df <- melt(df, value.name = "Counts")
df$gt <- factor(rep(c("WT","clv"), each = 10), levels = c("WT", "clv"))
df$trt <- factor(rep(c("Ctrl","BCN"), each = 5), levels = c("Ctrl", "BCN"))

ggplot(df, aes(trt, Counts)) +
  #geom_line(aes(group = gt), col = "red", lwd = 1, linetype = 1, stat = "summary", fun.y = "mean") +
  geom_jitter(width = 0.1, size = 2, col = "gray25") +
  facet_grid(. ~ gt) + 
  theme(axis.title.x = element_blank(),   # hide x axis title.
            axis.ticks.x = element_line(size = 0.5, color = "black"), #hide x axis ticks.
            axis.title.y = element_text(face = "plain", color = "black", size = 16),
            axis.ticks.y = element_line(size = 0.5, color = "black"),
            axis.text.x = element_text(vjust = 1, color = "black", size = 16),
            axis.text.y = element_text(face = "plain", color = "black", size = 12),
            axis.line = element_line(color = "black", size = 0.5, linetype = 1, lineend = "square"),
            strip.text.x = element_text(color = "black", size = 16),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            plot.background = element_rect(fill = "white"),
            plot.title = element_text(hjust = 0.5)) 
ggsave(file.path(maplot.path, "Example3.png"), device = "png", width = 6, height = 4, unit = "in", dpi = 300)

png("./legend.png", width = 2, height = 5, units = "in", res = 300)
plot(x = c(1,1,1,1,1,1), y = c(1,2,3,4,5,10), type = "p", pch = 19, cex = 2, 
     col = c("gray60", "lightblue", colRamp("palegreen", 1.3), colRamp("orange",2), "red","white"))
dev.off()


