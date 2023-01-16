source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")

## settings
gffa3.en <- file.path(bill.anno$path, bill.anno$gff3)
design = ~ Genotype + Treatment + Genotype:Treatment

## output files
out.path <- "./output/Bill.Analysis/XL.Rebuild2"
if(!file.exists(out.path)) dir.create(out.path)

## build dds objects

#dds.BCN <- ddsSlice(dds.EWRd75, factorName = "Treatment", factorLevel = c("Con", "Wrm"))

## get results
## csv file for each comparison will be saved to out_path 

DEseq2contrast(dds.LCM,  rebuild.LCM$path, ref_levels = c("WT", "Con"), countThreshold = 0, overwrite = T)
DEseq2contrast(dds.EWRd75, rebuild.EWRd75$path, ref_levels = c("WT", "Con"), countThreshold = 0, overwrite = T)

DEseq2contrast(dds.CLE,  out.path, ref_levels = c("WT", "Con"), countThreshold = 0, overwrite = T)
DEseq2contrast(dds.BCN, out.path, ref_levels = c("WT", "Con"), countThreshold = 0, overwrite = T)

## qcplots
### plot of CLV1, CLV2, and PRK2
genes <- as.data.frame(mcols(dds.EWR))
genes$id <- rownames(genes)
goi <- filter(genes, grepl("RPK2|CLV1|CLV2", Gene_description, ignore.case = T))
plotCounts2(dds.EWR, goi$id)
plotCounts2(dds.EWRd75, goi$id)
plotCounts2(dds.LCM, goi$id)

### sample distance and PCA
### rlog2transformation
rld.EWR <- rlogTransformation(dds.EWR)
rld.EWRd75 <- rlogTransformation(dds.EWRd75)
rld.LCM <- rlogTransformation(dds.LCM)

## qcplots
qcPlots(dds.EWR, file.path(rebuild.EWR$path, "qcplots"), addname="fullset", type="density", use = "ppt")
qcPlots(dds.EWR, file.path(rebuild.EWR$path, "qcplots"), addname="fullset", type ="ecdf", use = "ppt")
qcPlots(dds.EWR, file.path(rebuild.EWR$path, "qcplots"), addname="fullset", type="heatmap", use = "ppt")
qcPlots(dds.EWR, file.path(rebuild.EWR$path, "qcplots"), addname="fullset", ref_levels = c("WT", "Con"), type="PCA", subset="all", use = "display")
qcPlots(dds.EWR, file.path(rebuild.EWR$path, "qcplots"), addname="fullset", ref_levels = c("WT", "Con"), type="genecluster", subset=5000, use = "ppt")

qcPlots(dds.EWRd75, file.path(rebuild.EWR$path, "qcplots"), type="density", use = "ppt")
qcPlots(dds.EWRd75, file.path(rebuild.EWR$path, "qcplots"), type ="ecdf", use = "ppt")
qcPlots(dds.EWRd75, file.path(rebuild.EWR$path, "qcplots"), type="heatmap", use = "ppt")
qcPlots(dds.EWRd75, file.path(rebuild.EWR$path, "qcplots"), type="PCA", subset="all", use = "display")
qcPlots(dds.EWRd75, file.path(rebuild.EWR$path, "qcplots"), type="genecluster", subset=5000, use = "ppt")

qcPlots(dds.LCM, file.path(rebuild.EWR$path, "qcplots"), type="density", use = "ppt")
qcPlots(dds.LCM, file.path(rebuild.EWR$path, "qcplots"), type ="ecdf", use = "ppt")
qcPlots(dds.LCM, file.path(rebuild.EWR$path, "qcplots"), type="heatmap", use = "display")
qcPlots(dds.LCM, file.path(rebuild.EWR$path, "qcplots"), type="PCA", subset="all", use = "ppt")
qcPlots(dds.LCM, file.path(rebuild.EWR$path, "qcplots"), type="genecluster", subset=5000, use = "ppt")

