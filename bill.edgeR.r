## edgeR
## 
### bill's analysis results
library(edgeR)
source("./scripts/_edgeR.vars.r")

out.path <- "./output/edgeR"
res.path <- file.path(out.path, "results")
if(!file.exists(out.path)) dir.create(out.path)
if(!file.exists(res.path)) dir.create(res.path)

geneInfo <- read.table(file.path(anno$path, anno$annotation), 
                       header = T, sep = "\t", quote = "", stringsAsFactors = F)

cond <- read.table(file.path(counts$path,counts$EWRconditions), header = T)
cond$Genotype <- factor(cond$Genotype, levels = c("WT", "MU"))
cond$Treatment <- factor(cond$Treatment, levels = c("Con", "CLE", "Wrm"))
group <- paste(cond$Genotype, cond$Treatment, sep = ".")
group <- factor(group, levels = unique(group))
design <- model.matrix( ~0+group)
colnames(design) <- levels(group)

count.table <- as.matrix(read.table(file.path(counts$path, counts$EWRcnts)))


y <- DGEList(count.table, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = F]
y <- calcNormFactors(y)
y <- estimateDisp(y, design = design)
mds <- plotMDS(y, labels = colnames(y), top = 500, pch = 19)


fit <- glmQLFit(y, design)
contrasts.mat <- makeContrasts(MU.vs.WT = MU.Con - WT.Con,
                               WT_CLE.vs.WT_Con = WT.CLE - WT.Con,
                               WT_Wrm.vs.WT_Con = WT.Wrm - WT.Con,
                               MU_CLE.vs.MU_Con = MU.CLE - MU.Con,
                               MU_Wrm.vs.MU_Con = MU.Wrm - MU.Con,
                               interactionCLE = (MU.CLE - MU.Con) - (WT.CLE - WT.Con),
                               interactionBCN = (MU.Wrm - MU.Con) - (WT.Wrm - WT.Con),
                               levels = design)

for(comp in colnames(contrasts.mat)) {
  file.name <- sprintf("%s.csv", comp)
  qlf <- glmQLFTest(fit, contrast = contrasts.mat[, comp])
  res <- topTags(qlf, n = sum(keep))
  restable <- res$table
  gene_anno <- geneInfo[match(row.names(restable), geneInfo$Gene_id), ]
  restable <- cbind(restable, gene_anno)
  
  write.csv(restable, file.path(res.path, file.name), row.names = T, quote = F)
}
