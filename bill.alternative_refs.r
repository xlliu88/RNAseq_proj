##qPCR ref genes## controls
source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_celltype.fun.r")
source("./scripts/_importResult.bill.r")
library(RColorBrewer)
library(writexl)


alt.path <- "./output/Bill.Analysis/Bill_Alt"
EIS.res <- importResults(bill.resfile.EISd75)
LCM.res <- importResults(bill.resfile.LCM)
SD2009_path <- "./syncytia2009/rebuild"
SD2009_datafiles <- c("result.Syn_5-Ctrl.csv", 
                      "result.Syn_15-Ctrl.csv", 
                      "result.Syn_15-Syn_5.csv", 
                      "result.Syn-Ctrl.csv")
ref.lab <- c("AT3G04120","AT4G05320","AT5G25760","AT2G28390","AT4G26410","AT4G33380","AT4G34270")

ref.ids0 <- lapply(EIS.res, function(df) df$Gene_id[df$padj > 0.1 & df$baseMean > 1000])
ref.ids1 <- lapply(LCM.res, function(df) df$Gene_id[df$padj > 0.1 & df$baseMean > 1000])
ref.ids0 <- calculate.overlap2(ref.ids0)
ref.ids1 <- calculate.overlap2(ref.ids1)
ref.ids <- intersect(ref.ids0, ref.ids1)
ref.ids <- unique(c(ref.ids, ref.lab))
length(ref.ids)

#cat(ref.ids)
ids.idx <- EIS.res[[1]]$Gene_id %in% ref.ids
cols_sel <- c("Gene_id", "Gene_name", "baseMean","log2FoldChange", "padj", "Gene_description")
ids.geneInfo <- lapply(EIS.res, function(df) df[ids.idx, cols_sel])
#ids.geneInfo$Gene_description <- sub(",", ";", ids.geneInfo$Gene_description)
ids.geneInfo <- lapply(ids.geneInfo, function(df) {
   df$Gene_description <- sub("[Source:.*]$", "", df$Gene_description)
   return(df)})
write_xlsx(ids.geneInfo, file.path(alt.path, "ref_gene_candidate_expression_p0.1bM1000.xlsx"))


ids.info.LCM <- lapply(LCM.res, function(df) df[df$Gene_id %in% ref.ids, cols_sel])
write_xlsx(ids.info.LCM, file.path(alt.path, "ref_gene_candidate_expression_LCMp0.1bM1000.xlsx"))
## import 2009 Syncytium microarray data
SD2009 <- lapply(1:4, function(x) read.csv(file.path(SD2009_path, SD2009_datafiles[x]), header = T, stringsAsFactors = F))
SD2009 <- lapply(SD2009, function(df) df[order(df$X), ])


goi <- ref.ids #unique(c(goi1, ref.lab))
length(goi)

idx.SD <- sapply(goi, function(x) grep(x,SD2009[[1]]$locus))
idx.SD <- unlist(idx.SD)
data.cols <- c("X", "AveExpr", "logFC", "adj.P.Val", "locus")
SD.ref_expr <- lapply(SD2009, function(df) df[idx.SD, ])
names(SD.ref_expr) <- c("Syn5.vs.Ctrl", "Syn15.vs.Ctrl", "Syn15.vs.Syn5", "Syn.vs.Ctrl")
write_xlsx(SD.ref_expr, file.path(alt.path, "ref_gene_Expression_in_SD2009_p0.1bm1000.xlsx"))


## variaces
cnts <- counts(dds.EISd75, normalized = T)
ref_cnts <- cnts[row.names(cnts) %in% ref.ids, ]
ref_cnts2 <- log2(ref_cnts)

df <- data.frame(bM = rowMeans(ref_cnts), Var = rowVars(ref_cnts), SD = rowSds(ref_cnts))
df2 <- data.frame(bM = rowMeans(ref_cnts2), Var = rowVars(ref_cnts2), SD = rowSds(ref_cnts2))
fit <- lm(SD ~ bM + 0, data = df)
fit2 <- lm(SD ~ bM, data = df2)

pdf(file.path(alt.path, "qPCR", "ref_gene_sd.vs.baseMean(log2).pdf"), width = 12, height = 8)
plot(df2$bM-mean(df2$bM), df2$SD, pch = 19,ylim = c(0,1), xlab = "centrized baseMean", ylab = "sd", main = "BaseMean .vs. SD (log2 transformed)")
text(x = df2$bM - mean(df2$bM), y = df2$SD * 1.1, label = row.names(df), cex = 0.75)
abline(fit2$coefficients[1] + mean(df2$bM)*fit2$coefficients[2], fit2$coefficients[2], col = "gray")
dev.off()

df <- df[order(df$bM), ]
sel.ids <- c("AT4G33380", "AT2G28390", "AT4G29040", "AT1G76810", "AT2G05710", "AT1G65930", "AT3G04120")
isSelected <- sapply(row.names(df), function(x) x %in% sel.ids)
cols <- ifelse(isSelected, "red", "gray25")
pdf(file.path(alt.path, "qPCR", "ref_gene_sd.vs.baseMean.pdf"), width = 12, height = 8)
plot(df$bM, df$SD, pch = 19,
     xlim = c(0,max(df$bM)), 
     ylim = c(0,max(df$SD)), 
     col = cols,
     xlab = "centrized baseMean", 
     ylab = "sd", 
     main = "BaseMean .vs. SD")
text(x = df$bM,, y = df$SD * 1.1, label = row.names(df), cex = 0.5, srt = 0)
abline(0, fit$coefficients[1], col = "gray")
dev.off()
