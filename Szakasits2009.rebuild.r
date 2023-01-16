## to rebuild the Szakasits data
## and compare to RNAseq dataset
## 
source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")
source("./scripts/_celltype.vars.r")
source("./scripts/_celltype.fun.r")

## variables
proj.path <- "syncytia2009"
out.path <- file.path(proj.path, "rebuild")
image.path <- file.path(proj.path, "images")
if(!file.exists(out.path)) dir.create(out.path)
if(!file.exists(image.path)) dir.create(image.path)
expr.sd <- read.table(file.path(proj.path, "gcrma.qq.PLM.expression.dat"), 
                      header = T, sep = "\t", stringsAsFactors = F, fill = T)
deg.sd <- read.csv(file.path(proj.path, "tableS1.csv"), header = T, stringsAsFactors = F)

p.thd <- 0.001
fc.thd <- 1.2

cel.path <- file.path(proj.path, "CEL")
tar.file <- file.path(proj.path, "samples.szakasits2009.csv")
tar <- readTargets(tar.file, sep = ",")
row.names(tar) <- paste(tar$Description, tar$Rep, sep = ".")
tar$Treatment <- unlist(sapply(tar$Description, function(x) unlist(str_split(x, "_"))[1]))

## read affy data to affy batch data
files <- file.path(cel.path, tar$File.name) 
cat(sprintf("reading microarray data: %d files\n", length(files)))
ab.raw <- ReadAffy(filenames = files)
# aff.info <- compute.affinities.local(ab.raw)
ab.bg <- bg.adjust.gcrma(ab.raw, affinity.source = "local") # background correction
ab.bgqq <- normalize.AffyBatch.quantiles(ab.bg) # noromalization

# save images of raw data, bgcorrected, and normalized data
pdf(file.path(image.path, "Images.raw.pdf"))
for (i in 1:length(sampleNames(ab.raw))) {
  image(ab.raw[,i])
}
dev.off()
pdf(file.path(image.path, "Hist.raw.pdf"))
hist(ab.raw)
dev.off()

pdf(file.path(image.path, "Images.bgAdjusted.pdf"))
for (i in 1:length(sampleNames(ab.bg))) {
  image(ab.bg[,i])
}
dev.off()
pdf(file.path(image.path, "Hist.bgAdjusted.pdf"))
hist(ab.bg)
dev.off()

pdf(file.path(image.path, "Images.qq.pdf"))
for (i in 1:length(sampleNames(ab.bgpp))) {
  image(ab.bgpp[,i])
}
dev.off()
pdf(file.path(image.path, "Hist.qq.pdf"))
hist(ab.bgpp)
dev.off()

plm <- fitPLM(ab.bgqq, background = F, normalize = F)
eset.plm <- PLMset2exprSet(plm)

mt <- match(sampleNames(eset.plm), tar$File.name)
row.names(probe2gene) <- probe2gene$array_element_name
fmt <- match(row.names(eset.plm), row.names(probe2gene))
# fData(ab.raw) <- probe2gene[fmt,]

phenoData(eset.plm) <- new("AnnotatedDataFrame", data = tar[mt, ])
protocolData(eset.plm) <- new("AnnotatedDataFrame", data = tar[mt, ])
fData(eset.plm) <- probe2gene[fmt,]

## MA plot
jpeg(file.path(image.path, "MA.pairs.raw.jpg"), width = 2048, height = 2048, pointsize = 10)
mva.pairs(exprs(ab.raw), log.it = T)
dev.off()

jpeg(file.path(image.path, "MA.pairs.gcrma.qq.plm.jpg"), width = 2048, height = 2048, pointsize = 10)
mva.pairs(exprs(eset.plm), log.it = F)
dev.off()

## normalization
# eset.raw <- call.exprs(ab.raw, algorithm = "gcrma")
# write.exprs(eset, file.path(out.path,"expression.matrix(gcrma).txt"))
expr.mat <- exprs(eset.plm)

design <- model.matrix(~ 0 + as.factor(tar$Description))
colnames(design) <- levels(as.factor(tar$Description))
fit <- lmFit(eset.plm, design)

qlist <- c("Syn_5-Ctrl", "Syn_15-Ctrl", "Syn_15-Syn_5")
contrast.matrix <- makeContrasts(contrasts = qlist, levels = fit$design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
for(q in qlist) {
  res <- topTable(fit2, coef = q, number = nrow(fit2))
  res <- select(res, AveExpr, logFC, t, P.Value, adj.P.Val,B, locus, description)
  res <- res[order(res$logFC, decreasing = T),]
  
  pair <- unlist(str_split(q, "-"))
  for(g in pair) {
    g.exp <- expr.mat[, which(phenoData(eset.plm)$Description == g)]
    g.exp <- g.exp[match(row.names(res), row.names(g.exp)),]
    res[[sprintf("%s.mean",g)]] <- apply(g.exp, 1, mean)
  }
  
  res$description <- gsub(",", ";", res$description)
  res$description <- gsub("\'", "", res$description)
  res$description <- gsub("\"", "", res$description)
  
  write.csv(res, file.path(out.path, sprintf("result.%s.csv", q)), row.names = T, quote = F)
}


design.a <- model.matrix(~ 0 + factor(tar$Treatment))
colnames(design.a) <- levels(as.factor(tar$Treatment))
fit.a <- lmFit(eset.plm, design.a)
contrast.matrix.a <- makeContrasts(Syn-Ctrl, levels = fit.a$design)
fit2.a <- contrasts.fit(fit.a, contrast.matrix.a)
fit2.a <- eBayes(fit2.a)
res.a <- topTable(fit2.a, coef = 1, number = nrow(fit2.a))

res.a <- select(res.a, AveExpr, logFC, t, P.Value, adj.P.Val,B, locus, description)
res.a <- res.a[order(res.a$logFC, decreasing = T),]

res.a$description <- gsub(",", ";", res.a$description)
res.a$description <- gsub("\'", "", res.a$description)
res.a$description <- gsub("\"", "", res.a$description)

## add mean expression level of each group
pair <- levels(as.factor(tar$Treatment))
for(g in pair) {
  g.exp <- expr.mat[, which(phenoData(eset.plm)$Treatment == g)]
  g.exp <- g.exp[match(row.names(res.a), row.names(g.exp)),]
  res.a[[sprintf("%s.mean",g)]] <- apply(g.exp, 1, mean)
}

write.csv(res.a, file.path(out.path,"result.Syn-Ctrl.csv"), row.names = T, quote = F)


## compare of expression matrix from
## 1. matrix from the paper. - file "gcrma.qq.PLM.expression.dat"
## 2. rebuilt data. gcrma normalized