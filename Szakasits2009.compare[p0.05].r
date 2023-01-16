# SD2009PJ vs RNAseq
# 
source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")

fc.thd <- 1
p.thd <- 0.05
sufix <- sprintf("[fc%s.p%s]", fc.thd, formatC(p.thd, format = "e", digits = 0))
proj.path <- "syncytia2009"
rebuild.path <- file.path(proj.path, "rebuild")
out.path <- file.path(proj.path, "Comparisons[p0.05]")
if(!file.exists(out.path)) dir.create(out.path)
affy.array.file <- file.path("rootcells2007", "affy_ATH1_array_elements-2010-12-20.txt")
probe2gene <- read.delim(affy.array.file, header = T, stringsAsFactors = F, sep = "\t", quote = "")

expr.sd <- read.table(file.path(proj.path, "gcrma.qq.PLM.expression.dat"), 
                      header = T, sep = "\t", stringsAsFactors = F, fill = T)
res.sd <- read.csv(file.path(proj.path, "tableS1.csv"), row.names = 1,header = T, stringsAsFactors = F)

expr.rebuild <- read.table(file.path(rebuild.path, "expression.matrix(gcrma).txt"), header = T, stringsAsFactors = F)
res.rebuild <- read.csv(file.path(rebuild.path, "result.Syn-Ctrl.csv"), row.names = 1, header = T, stringsAsFactors = F)

probset.sd <- probe2gene[match(row.names(res.sd), probe2gene$locus),]

# match geneID to probeset in rebuilded result
res.rebuild$GeneID <- probe2gene$locus[match(row.names(res.rebuild), probe2gene$array_element_name)]
res.rebuild <- res.rebuild[!grepl(";", res.rebuild$GeneID),]
# ## make sure GeneID correctly matched
# px <- probe2gene[probe2gene$array_element_name %in% row.names(res.rebuild),]
# res.rebuild <- res.rebuild[order(row.names(res.rebuild)),]
# px <- px[order(px$array_element_name),]
# all(px$array_element_name == row.names(res.rebuild))
# all(px$locus == res.rebuild$GeneID)

## removed duplicated GeneID, keep the one with larger adj.p
res.rebuild <- res.rebuild[order(res.rebuild$adj.P.Val, decreasing = T), ]
res.rebuild <- res.rebuild[!duplicated(res.rebuild$GeneID), ]
res.rebuild.clean <- res.rebuild

## overlapped genes in SD2009PJ result and rebuild result
# vlist <- list(sd = row.names(res.sd), rebuild = res.rebuild$GeneID)
# venn.diagram(vlist, 
#              filename = file.path(out.path,"GeneSet.SD.vs.XL.png"), 
#              fill = brewer.pal(3, "Set2")[1:2])
# 
# commgenes <- intersect(row.names(res.sd), res.rebuild$GeneID)
# res.sd <- res.sd[row.names(res.sd) %in% commgene, ]
# res.sd <- res.sd[order(row.names(res.sd)),]
# res.rebd <- res.rebuild[res.rebuild$GeneID %in% commgene, ]
# res.rebd <- res.rebd[order(res.rebd$GeneID), ]
# all(row.names(res.sd) == res.rebd$GeneID)
# 
# ## compare of logFC and q.value
# png(file.path(rebuild.path, "logFC.SD.vs.XL.png"), width = 1024, height = 1024)
# plot(res.sd$M, res.rebd$logFC, xlab = "SD2009PJ", ylab = "XL.rebuild", pch = 19, cex = 0.75)
# dev.off()
# png(file.path(rebuild.path, "logFC.XL.vs.SD.png"), width = 1024, height = 1024)
# plot(res.rebd$logFC, res.sd$M , ylab = "SD2009PJ", xlab = "XL.rebuild", pch = 19, cex = 0.75)
# dev.off()
# 
# png(file.path(rebuild.path, "q.value.SD.vs.XL.png"), width = 800, height = 800)
# plot(-log10(res.sd$adj.q), -log10(res.rebd$adj.P.Val), xlab = "SD2009PJ", ylab = "XL.rebuild", pch = 19, cex = 0.75)
# dev.off()

## up and down regulated genes
sd.up <- res.sd[res.sd$M > log2(fc.thd) & res.sd$adj.q < p.thd, ]
# rebd.up <- res.rebd[res.rebd$logFC > log2(fc.thd) & res.rebd$adj.P.Val < p.thd, ]

sd.dn <- res.sd[res.sd$M < -log2(fc.thd) & res.sd$adj.q < p.thd, ]
# rebd.dn <- res.rebd[res.rebd$logFC < -log2(fc.thd) & res.rebd$adj.P.Val < p.thd, ]

# direction <- "up"
# fn <- file.path(rebuild.path,sprintf("DEG-%s%s.SD.vs.XL.png", direction, sufix))
# vlist.up <- list(SD = row.names(sd.up),
#                  XL = rebd.up$GeneID)
# venn.diagram(vlist.up, 
#              filename = fn, 
#              fill = brewer.pal(3, "Set2")[1:2])
# 
# direction <- "dn"
# fn <- file.path(rebuild.path,sprintf("DEG-%s%s.SD.vs.XL.png", direction, sufix))
# vlist.dn <- list(SD = row.names(sd.dn),
#                  XL = rebd.dn$GeneID)
# venn.diagram(vlist.dn, 
#              filename = fn, 
#              fill = brewer.pal(3, "Set2")[1:2])
# 
# direction <- ""
# fn <- file.path(rebuild.path,sprintf("DEG-%s%s.SD.vs.Rebuild-Euler.png", direction, sufix))
# vlist.all <- list(SD.up = row.names(sd.up),
#                   SD.dn = row.names(sd.dn),
#                   Rebuild.up = rebd.up$GeneID, 
#                   Rebuild.dn = rebd.dn$GeneID)
# # venn.diagram(vlist.all, 
# #              filename = fn,
# #              fill = brewer.pal(4, "Set2"))
# png(fn, width = 800, height = 800)
# plot(euler(vlist.all), 
#      fill = colRamp(brewer.pal(4, "Set2"),1.5),
#      edges = colRamp(brewer.pal(4, "Set2"),1),
#      quantities = T, 
#      labels = T )
# dev.off()

## compare to RNAseq data
EIS.wbc <- importResults(bill.resfile.EISd75, "wbc")[["wbc"]]
LCM.wbc <- importResults(bill.resfile.LCM, "wbc")[["wbc"]]

## SD2009PJ result ###########################
res.sd <- read.csv(file.path(proj.path, "tableS1.csv"), row.names = 1,header = T, stringsAsFactors = F)
commgenes <- intersect(res.sd$GeneID, intersect(EIS.wbc$Gene_id, LCM.wbc$Gene_id))
EIS.wbc <- EIS.wbc %>% filter(Gene_id %in% commgenes) %>% arrange(Gene_id)
LCM.wbc <- LCM.wbc %>% filter(Gene_id %in% commgenes) %>% arrange(Gene_id)
res.sd <- res.sd %>% filter(GeneID %in% commgenes) %>% arrange(GeneID) 

## plot p value comparisons
sd.p <- -log10(res.sd$adj.q)
eis.p <- -log10(EIS.wbc$padj)
lcm.p <- -log10(LCM.wbc$padj)
col <- ifelse(sd.p > 2 & eis.p > 2, "red", "gray")
col <- ifelse(sd.p > 2 & eis.p <= 2, colRamp("orange", 1.5), col)
col <- ifelse(sd.p <=2 & eis.p > 2, "orange", col)
size <- abs(LCM.wbc$log2FoldChange/res.sd$M)
#size <- 1/size
size <- ifelse(size > 2, 2, size)
size <- ifelse(size < 0.5, 0.5, size)
pdf(file.path(out.path, "SD2009PJvsEIS.pvalues3.pdf"), width = 6, height = 6)
par(mar=c(5.1,5.1,2.1,2.1))
plot(sd.p, eis.p, xlab = "-log10(q.value) - Micro-aspiration",
     ylab = "-log10(q.value) - EIS", 
     cex.lab = 1, cex.axis = 1, 
     col = col, pch = 19, cex = size)
dev.off()

## 2009PJ vs EIS, logFoldChange
fe <- lm(EIS.wbc$log2FoldChange ~ res.sd$M)
# los <- loess(EIS.wbc$log2FoldChange ~ res.sd$M)
# d <- order(res.sd$M)
se <- summary(fe)
pdf(file.path(out.path, "SD2009PJvsEIS.logFC.fit2.pdf"), width = 6, height = 6)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot( EIS.wbc$log2FoldChange ~ res.sd$M,
     xlab = "log2FC - Micro-aspiration", ylab = "log2FC - EIS",
     xlim = c(-12, 12), ylim = c(-12, 12), 
     col = col, cex = 1/3, pch = 19,
     cex.lab = 1, cex.axis = 1)
#abline(fe, col = "gray20")
abline(se$coefficients[1,1], se$coefficients[2,1], col = "gray25")
# lines(res.sd$M[od], los$fitted[od], col = "blue", lty = 2)
# abline(0,1, col = "green", lty = 2)
abline(h = 0, col = "gray50", lty = 2)
abline(v = 0, col = "gray50", lty = 2)
text(x = 12,  y = se$coefficients[2,1]*12 + se$coefficients[1,1] + 0.5,
     labels = sprintf("R-squared\n%1.4f", se$r.squared), 
     adj = c(1,1), cex = 1)
dev.off()

## LCM vs SD2009, p-values
col <- ifelse(sd.p > 2 & lcm.p > 2, "red", "gray")
col <- ifelse(sd.p > 2 & lcm.p <= 2, colRamp("orange", 1.5), col)
col <- ifelse(sd.p <= 2 & lcm.p > 2, "orange", col)
size <- abs(LCM.wbc$log2FoldChange/res.sd$M)
size <- ifelse(size > 2, 2, size)
size <- ifelse(size < 0.5, 0.5, size)
pdf(file.path(out.path, "SD2009PJvsLCM.pvalues2.pdf"), width = 6, height = 6)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(sd.p, lcm.p, xlab = "-log10(q.value) - Micro-aspiration", ylab = "-log10(q.value) - LCM", 
     cex.lab = 1, cex.axis = 1, 
     col = col, pch=19, cex=size)
dev.off()

fl <- lm(LCM.wbc$log2FoldChange ~ res.sd$M)
sl <- summary(fl)
pdf(file.path(out.path, "SD2009PJvsLCM.logFC.fit.pdf"), width = 6, height = 6)
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(LCM.wbc$log2FoldChange ~ res.sd$M, 
     xlab = "log2FC - Micro-aspiration", ylab = "log2FC - LCM",
     xlim = c(-12, 12), ylim = c(-12, 12), 
     col = col, cex = 1/3, pch = 19,
     cex.lab = 1, cex.axis = 1)
abline(fl, col = "gray20")
abline(h = 0, col = "gray50", lty = 2)
abline(v = 0, col = "gray50", lty = 2)
text(x = 12,  y = sl$coefficients[2,1]*12 + sl$coefficients[1,1] + 0.5,
     labels = sprintf("R-squared\n%1.4f", sl$r.squared), 
     adj = c(1,1), cex = 1)
dev.off()


## venn plot
EIS.up <- EIS.wbc %>% filter(log2FoldChange > log2(fc.thd), padj < p.thd)
LCM.up <- LCM.wbc %>% filter(log2FoldChange > log2(fc.thd), padj < p.thd)
EIS.dn <- EIS.wbc %>% filter(log2FoldChange < -log2(fc.thd), padj < p.thd)
LCM.dn <- LCM.wbc %>% filter(log2FoldChange < -log2(fc.thd), padj < p.thd)

SD.up <- res.sd[res.sd$M > log2(fc.thd) & res.sd$adj.q < p.thd,]
SD.dn <- res.sd[res.sd$M < -log2(fc.thd) & res.sd$adj.q < p.thd,]

direction <- "up"
fn <- file.path(out.path, sprintf("DEG-%s%s.SD2009PJ.vs.RNAseq_4x4.pdf", direction, sufix))
#fn2 <- file.path(out.path, sprintf("DEG-%s%s.SD2009PJ.vs.RNAseq_labeless.pdf", direction, sufix))
vlist.up <- list(EIS =EIS.up$Gene_id,
                 LCM = LCM.up$Gene_id,
                 `Micro-aspiration` = SD.up$GeneID)
e.up <- euler(vlist.up)

pdf(fn, width = 4, height = 4)
plot(e.up, 
     fill = colRamp(brewer.pal(3, "Set2"),1.5),
     edges = colRamp(brewer.pal(3, "Set2"),1),
     quantities = T, 
     labels = T )#identical(legend, T))
dev.off()


direction <- "dn"
fn <- file.path(out.path, sprintf("DEG-%s%s.SD2009PJ.vs.RNAseq_4x4.pdf", direction, sufix))
vlist.dn <- list(EIS = EIS.dn$Gene_id, #row.names(EIS.dn),
                 LCM = LCM.dn$Gene_id, #row.names(LCM.dn),
                 `Micro-aspiration` = SD.dn$GeneID) #row.names(SD.dn))
e.dn <- euler(vlist.dn)
pdf(fn, width = 4, height = 4)
plot(euler(vlist.dn), 
     fill = colRamp(brewer.pal(3, "Set2"),1.5),
     edges = colRamp(brewer.pal(3, "Set2"),1),
     quantities = T, 
     labels = T )#identical(legend, T))
dev.off()


## write overlapped genes to dataframe
ovlp.up <- calculate.overlap(vlist.up)
ovlp.dn <- calculate.overlap(vlist.dn)
sd.df <- res.sd[res.sd$GeneID %in% c(ovlp.up$a5,ovlp.dn$a5), ]
eis.df <- EIS.wbc[EIS.wbc$Gene_id %in% c(ovlp.up$a5, ovlp.dn$a5), ]
lcm.df <- LCM.wbc[LCM.wbc$Gene_id %in% c(ovlp.up$a5, ovlp.dn$a5), ]
write.csv(sd.df, file.path(out.path, "[SD2009PJvsEISvsLCM].Co-regulated-Genes_in_SD2009PJ.csv"), row.names = T, quote = F)
write.csv(eis.df, file.path(out.path, "[SD2009PJvsEISvsLCM].Co-regulated-Genes_in_EIS.csv"), row.names = T, quote = F)
write.csv(lcm.df, file.path(out.path, "[SD2009PJvsEISvsLCM].Co-regulated-Genes_in_LCM.csv"), row.names = T, quote = F)

## compare to RNAseq data - Rebuild data
EIS.wbc <- importResults(bill.resfile.EISd75, "wbc")[["wbc"]]
LCM.wbc <- importResults(bill.resfile.LCM, "wbc")[["wbc"]]
res.rebuild <- read.csv(file.path(rebuild.path, "result.Syn-Ctrl.csv"), row.names = 1, header = T, stringsAsFactors = F)
res.rebuild <- res.rebuild[!grepl(";", res.rebuild$locus), ]
res.rebuild <- res.rebuild[order(res.rebuild$adj.P.Val, decreasing = T),]
res.rebuild <- res.rebuild[!duplicated(res.rebuild$locus),]
row.names(res.rebuild) <- res.rebuild$locus

## rebuild result ###########################
commgenes <- intersect(row.names(res.rebuild), intersect(row.names(EWR.wbc), row.names(LCM.wbc)))
EWR.wbc <- EWR.wbc[row.names(EWR.wbc) %in% commgenes, ]
LCM.wbc <- LCM.wbc[row.names(LCM.wbc) %in% commgenes, ]
res.rebuild <- res.rebuild[row.names(res.rebuild) %in% commgenes, ]
EWR.wbc <- EWR.wbc[order(row.names(EWR.wbc)),]
LCM.wbc <- LCM.wbc[order(row.names(LCM.wbc)),]
res.rebuild <- res.rebuild[order(row.names(res.rebuild)),]
## scatter plot of corresponding pvalue
rebd.p <- -log10(res.rebuild$adj.P.Val)
eis.p <- -log10(EWR.wbc$padj)
lcm.p <- -log10(LCM.wbc$padj)
col <- ifelse(rebd.p > 2 & eis.p > 2, "red", "gray")
col <- ifelse(rebd.p > 2 & eis.p <= 2, colRamp("orange", 1.5), col)
col <- ifelse(rebd.p <=2 & eis.p > 2, "orange", col)
png(file.path(out.path, "RebuildvsEIS.pvalues.png"), width = 800, height = 800)
plot(rebd.p, eis.p, xlab = "Rebuild -log10(q.value)", ylab = "EIS -log10(q.value)", col = col, pch=19, cex=0.7)
dev.off()

png(file.path(out.path, "RebuildvsEIS.logFC.png"), width = 800, height = 800)
plot(res.rebuild$logFC, EWR.wbc$log2FoldChange, 
     xlab = "log2FC(Rebuild)", ylab = "log2FC(EIS)",
     col = col, cex = 0.7, pch = 19)
dev.off()

col <- ifelse(rebd.p > 2 & lcm.p > 2, "red", "gray")
col <- ifelse(rebd.p > 2 & lcm.p <= 2, colRamp("orange", 1.5), col)
col <- ifelse(rebd.p <= 2 & lcm.p > 2, "orange", col)
size <- abs(LCM.wbc$log2FoldChange/res.rebuild$logFC)
size <- ifelse(size > 10, 10, size)
png(file.path(out.path, "RebuildvsLCM.pvalues.png"), width = 800, height = 800)
plot(rebd.p, lcm.p, xlab = "Rebuild -log10(q.value)", ylab = "LCM -log10(q.value)", col = col, pch=19, cex=0.7)
dev.off()

png(file.path(out.path, "RebuildvsLCM.logFC.png"), width = 800, height = 800)
plot(res.rebuild$logFC, LCM.wbc$log2FoldChange, 
     xlab = "log2FC(Rebuild)", ylab = "log2FC(LCM)",
     col = col, cex = 0.7, pch = 19)
dev.off()



## DEGs
EIS.up <- EWR.wbc[EWR.wbc$log2FoldChange > log2(fc.thd) & EWR.wbc$padj < p.thd, ]
LCM.up <- LCM.wbc[LCM.wbc$log2FoldChange > log2(fc.thd) & LCM.wbc$padj < p.thd, ]
EIS.dn <- EWR.wbc[EWR.wbc$log2FoldChange < -log2(fc.thd) & EWR.wbc$padj < p.thd, ]
LCM.dn <- LCM.wbc[LCM.wbc$log2FoldChange < -log2(fc.thd) & LCM.wbc$padj < p.thd, ]

rbd.up <- res.rebuild[res.rebuild$logFC > log2(fc.thd) & res.rebuild$adj.P.Val < p.thd,]
rbd.dn <- res.rebuild[res.rebuild$logFC < -log2(fc.thd) & res.rebuild$adj.P.Val < p.thd,]



direction <- "up"
fn <- file.path(out.path, sprintf("DEG-%s%s.rebuild.vs.RNAseq.jpg", direction, sufix))
fn2 <- file.path(out.path, sprintf("DEG-%s%s.rebuild.vs.RNAseq_labeless.jpg", direction, sufix))
vlist.up <- list(EIS = row.names(EIS.up),
                 LCM = row.names(LCM.up),
                 rebuild = row.names(rbd.up))
# venn.diagram(vlist.up, 
#              filename = fn,
#              fill = brewer.pal(3, "Set2"))
png(fn, width = 800, height = 800)
plot(euler(vlist.up), 
     fill = colRamp(brewer.pal(3, "Set2"),1.5),
     edges = colRamp(brewer.pal(3, "Set2"),1),
     quantities = T, 
     labels = T )#identical(legend, T))
dev.off()

png(fn2, width = 800, height = 800)
plot(euler(vlist.up), 
     fill = colRamp(brewer.pal(3, "Set2"),1.5),
     edges = colRamp(brewer.pal(3, "Set2"),1),
     quantities = T, 
     labels = F )#identical(legend, T))
dev.off()

direction <- "dn"
fn <- file.path(out.path, sprintf("DEG-%s%s.rebuild.vs.RNAseq.jpg", direction, sufix))
fn2 <- file.path(out.path, sprintf("DEG-%s%s.rebuild.vs.RNAseq_labeless.jpg", direction, sufix))
vlist.dn <- list(EIS = row.names(EIS.dn),
                 LCM = row.names(LCM.dn),
                 Rebuild = row.names(rbd.dn))
# venn.diagram(vlist.dn, 
#              filename = fn,
#              fill = brewer.pal(3, "Set2"))

png(fn, width = 800, height = 800)
plot(euler(vlist.dn), 
     fill = colRamp(brewer.pal(3, "Set2"),1.5),
     edges = colRamp(brewer.pal(3, "Set2"),1),
     quantities = T, 
     labels = T )#identical(legend, T))
dev.off()
png(fn2, width = 800, height = 800)
plot(euler(vlist.dn), 
     fill = colRamp(brewer.pal(3, "Set2"),1.5),
     edges = colRamp(brewer.pal(3, "Set2"),1),
     quantities = T, 
     labels = F )#identical(legend, T))
dev.off()

## write overlapped gene list
ovlp.up <- calculate.overlap(vlist.up)
ovlp.dn <- calculate.overlap(vlist.dn)
rebd.df <- res.rebuild[res.rebuild$locus %in% c(ovlp.up$a5, ovlp.dn$a5), ]
eis.df <- EWR.wbc[row.names(EWR.wbc) %in% c(ovlp.up$a5, ovlp.dn$a5), ]
lcm.df <- LCM.wbc[row.names(LCM.wbc) %in% c(ovlp.up$a5, ovlp.dn$a5), ]
write.csv(rebd.df, file.path(out.path, "[RebuildvsEISvsLCM].Co-regulated-Genes_in_Rebuild.csv"), row.names = T, quote = F)
write.csv(eis.df, file.path(out.path, "[RebuildvsEISvsLCM].Co-regulated-Genes_in_EIS.csv"), row.names = T, quote = F)
write.csv(lcm.df, file.path(out.path, "[RebuildvsEISvsLCM].Co-regulated-Genes_in_LCM.csv"), row.names = T, quote = F)

### comparison of gene expression (log2(TPM), vs. log2(Expr)) #######
EIS.wbc <- importResults(bill.resfile.EISd75, "wbc")[["wbc"]]
LCM.wbc <- importResults(bill.resfile.LCM, "wbc")[["wbc"]]
TPM.EIS <- read.csv(file.path("./rootcells2007", "TPM_EIS.csv"), row.names = 1, header = T, stringsAsFactors = F)
TPM.LCM <- read.csv(file.path("./rootcells2007", "TPM_LCM.csv"), row.names = 1, header = T, stringsAsFactors = F)

EIS.ct <- read.table(file.path(bill.counts$path, bill.counts$EISd75conditions), 
                     head=T, row.names = 1, stringsAsFactors = T)
LCM.ct <- read.table(file.path(bill.counts$path, bill.counts$LCMconditions), 
                               head=T, row.names = 1, stringsAsFactors = T)
EIS.ct <- EIS.ct[match(colnames(TPM.EIS), row.names(EIS.ct)), ]
LCM.ct <- LCM.ct[match(colnames(TPM.LCM), row.names(LCM.ct)), ]
EIS.ctrl.mat <- TPM.EIS[, which(EIS.ct$Genotype == "WT" & EIS.ct$Treatment == "Ctrl")]
EIS.BCN.mat <- TPM.EIS[, which(EIS.ct$Genotype == "WT" & EIS.ct$Treatment == "Infected")]
LCM.ctrl.mat <- TPM.LCM[, which(LCM.ct$Genotype == "WT" & LCM.ct$Treatment == "Ctrl")]
LCM.BCN.mat <- TPM.LCM[, which(LCM.ct$Genotype == "WT" & LCM.ct$Treatment == "Infected")]

EIS.expr <- data.frame(Ctrl = log2(rowMeans(EIS.ctrl.mat)+1),
                          BCN = log2(rowMeans(EIS.BCN.mat)+1))
LCM.expr <- data.frame(Ctrl = log2(rowMeans(LCM.ctrl.mat)+1),
                          BCN = log2(rowMeans(LCM.BCN.mat)+1))

commgenes <- intersect(row.names(res.sd), intersect(row.names(EIS.expr), row.names(LCM.expr)))
EIS.expr <- EIS.expr[match(commgenes, row.names(EIS.expr)), ]
LCM.expr <- LCM.expr[match(commgenes, row.names(LCM.expr)), ]
res.sd <- res.sd[match(commgenes, row.names(res.sd)), ]

all(row.names(EIS.expr) == row.names(res.sd))
all(row.names(LCM.expr) == row.names(res.sd))
expr <- cbind(EIS.expr, LCM.expr)
#expr[["EWR.p"]] <- EWR.wbc$padj[match(commgenes, row.names(EWR.wbc))]
#expr[["LCM.p"]] <- LCM.wbc$padj[match(commgenes, row.names(LCM.wbc))]

colnames(expr) <- c("EIS-Ctrl", "EIS-Infected", "LCM-Ctrl", "LCM-Infected")
colnames(res.sd)[colnames(res.sd) == "Root"] <- "Micro-asp.-Ctrl"
colnames(res.sd)[colnames(res.sd) == "Sync"] <- "Micro-asp.-Infected"
for (s in colnames(expr)){
  sample <- unlist(str_split(s, "-"))
  s.sd <- sprintf("Micro-asp.-%s", sample[2])
  fexp <- lm(expr[[s]] ~ res.sd[[s.sd]]) 
  sexp <- summary(fexp)
  pdf(file.path(out.path, sprintf("Expression.SD2009PJvs%s_loess_4x4.pdf", s)), width = 4, height = 4)
  par(mar = c(5.1,5.1,2.1,2.1))
  plot(res.sd[[s.sd]], expr[[s]],
       main = "log2(MeanExpr)",cex.main = 1.5,
       xlab = s.sd,
       ylab = s,
       xlim = c(0,14), ylim = c(0,14),
       col = "gray25", 
       cex.lab = 2, cex.axis = 2, cex = 0.7, pch=19)
  lf <- loess.smooth(res.sd[[s.sd]], expr[[s]])
  lines(lf, col = colRamp("red",1.5), lwd = 1.5, lty = 2)
  abline(fexp, col = "red", lwd = 1)
  x <- max(res.sd[[s.sd]])
  text(x = 14,  y = sexp$coefficients[2,1]*14 + sexp$coefficients[1,1] + 1,
      labels = sprintf("R-squared\n%1.4f", sexp$r.squared),
      adj = c(1,1), cex = 1.5,
      col = "red")
  dev.off()
}

for (s in colnames(expr)){
  sample <- unlist(str_split(s, "-"))
  s.sd <- sprintf("Micro-asp.-%s", sample[2])
  fexp <- lm(expr[[s]] ~ res.sd[[s.sd]]) 
  sexp <- summary(fexp)
  pdf(file.path(out.path, sprintf("Expression.SD2009PJvs%s_smooth_4x4.pdf", s)), width = 4, height = 4)
  par(mar = c(5.1,5.1,2.1,2.1))
  smoothScatter(res.sd[[s.sd]], expr[[s]],
       main = "log2(MeanExpr)",cex.main = 1,
       xlab = s.sd,
       ylab = s,
       xlim = c(0,14), ylim = c(0,14),
       cex.lab = 1, cex.axis = 1, cex = 0.5, pch=19)
  lf <- loess.smooth(res.sd[[s.sd]], expr[[s]])
  lines(lf, col = colRamp("red", 2), lwd = 1.5)
  dev.off()
}


expr2 <- cbind(expr, res.sd[,1:2])
png(file.path(out.path, "Expression.Pairs.png"), width = 1024, height = 1024)
pairs(expr2, pch = ".")
dev.off()



lfc.dat <- list(EIS = EWR.wbc$log2FoldChange,
            LCM = LCM.wbc$log2FoldChange,
            MA = res.sd$M)

p.dat <- list(EIS = EWR.wbc$padj,
              LCM = LCM.wbc$padj,
              MA = res.sd$adj.q)

p2.dat <- lapply(p.dat, function(x) -log10(x))

mypar(1,3)
lapply(1:3, function(x) hist(lfc.dat[[x]], breaks = 100,xlim = c(-20, 20), 
                             xlab = "log2FoldChange", 
                             main = sprintf("Histogram for %s", names(lfc.dat)[x])))

lapply(1:3, function(x) hist(p.dat[[x]], breaks = 100, #xlim = c(-20, 20), 
                             xlab = "adjusted p-vaue", 
                             main = names(lfc.dat)[x]))
lapply(1:3, function(x) hist(p2.dat[[x]], breaks = 100, #xlim = c(-20, 20), 
                             xlab = "-log10(p-value)", 
                             main = names(lfc.dat)[x]))
