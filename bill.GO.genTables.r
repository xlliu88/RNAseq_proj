## modified from bill.GO.genTables.r
## use weight01 and Fisher's exact test

source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")
source("./scripts/_go_kegg.r")

########################################################
#                                                    ###
#      This Script Generates GO Enrichment Tables     ###
#                                                    ###
########################################################


pthreshold <- 0.05
go.path <- "./output/Bill.Analysis/topGO3"
if(!dir.exists(go.path)) dir.create(go.path)
fisher.path <- file.path(go.path, "fisherKS.tables")
if(!dir.exists(fisher.path)) dir.create(fisher.path)


### Rebuild DDS Objects
### import results
EIS.res <- importResults(bill.resfile.EISd75)
LCM.res <- importResults(bill.resfile.LCM)
LCM.BCN <- consolidateData(data = LCM.res, "BCN")
EIS.BCN <- consolidateData(EIS.res, "BCN")
EIS.CLE <- consolidateData(EIS.res, "CLE")


## calculater GO Enrichment in EIS and LCM individully and save to files ####
## use weight01 alg, and Fisher's exact test
goterms <- c("BP", "MF", "CC")
sampleSource <- c("LCM", "EIS")

LCM.go <- list()
EIS.go <- list()
for (source in sampleSource) {
  for (compname in names(compare.names)) {
      if(!(compname == "wbc")) next
      for (whichgo in goterms) {
        for (direction in c("up", "dn")) {
           
            fn.pfx <- paste(source, compname, direction, sep = ".")
            result.name <- paste(compname, direction, whichgo, sep = ".")
            if (source == "LCM" & compname %in% names(LCM.res)) {
               df <- LCM.res[[compname]]
            } else if (source == "EIS" & compname %in% names(EIS.res)) {
               df <- EIS.res[[compname]]
            } else {
               cat("no combination found; skip to next\n")
               next
            }
            
            if (direction == "up") {
               p <- with(df, ifelse(log2FoldChange > 0 & padj < pthreshold, padj, 1))
            } else if (direction == "dn") {
               p <- with(df, ifelse(log2FoldChange < 0 & padj < pthreshold, padj, 1))
            }
            
            names(p) <- df$Gene_id
            go.res <- goTair(p, whichGO = whichgo, topNodes = 30, stat = c("fisher", "KS"), 
                             out.path = fisher.path, fn.prefix = fn.pfx, outFormat = "all")
            # go.res$fdr <- p.adjust(go.res$weight01fisher, method = "fdr")
            
            if (source == "LCM") LCM.go[[result.name]] <- go.res
            if (source == "EIS") EIS.go[[result.name]] <- go.res
       }
    }
  }
}

## calculate GO Enrichment of EIS&LCM coregulated genes and save to file ####
LvsE.go <- list()
source <- "LnE"
for (compname in names(compare.names)) {
    for (whichgo in goterms) {
      for (direction in c("up", "dn")) {
        
        fn.pfx <- paste(source, compname, direction,"top20", sep = ".")
        result.name <- paste(compname, direction, whichgo, sep = ".")
        if (source == "LnE" & compname %in% names(LCM.res)) {
          df.LCM <- LCM.res[[compname]]
          df.EIS <- EIS.res[[compname]]
         } else {
          cat("no combination found; skip to next\n")
          next
        }
        
        df <- inner_join(df.LCM, df.EIS, by = c("Gene_id", "Gene_name")) %>% 
          mutate(padj = pmax(padj.x, padj.y)) %>% 
          mutate(lfc = pmin(abs(log2FoldChange.x), abs(log2FoldChange.y))) %>% 
          mutate(lfc = ifelse(log2FoldChange.x * log2FoldChange.y < 0, 0, lfc)) %>% 
          mutate(lfc = ifelse(log2FoldChange.x < 0, -lfc, lfc))
        
        if (direction == "up") {
          p <- with(df, ifelse(lfc < 0, 1, padj))
        } else if (direction == "dn") {
          p <- with(df, ifelse(lfc > 0, 1, padj))
        }
        
        names(p) <- df$Gene_id
        
        go.res <- goTair(p, whichGO = whichgo, topNodes = 20, stat = c("Fisher", "KS"),
                         outFormat = "all",
                         out.path = fisher.path, fn.prefix = fn.pfx)
        
        LvsE.go[[result.name]] <- go.res

      }
    }
}







###############################################################################################################################################
## go enrichment of genes overlapped in EIS.BCN and EIS.CLE ####
df.bcn <- EIS.res[["wbc"]]
df.cle <- EIS.res[["wec"]]
go.res <- list()

up.p1 <- with(df.bcn, ifelse(log2FoldChange > 0 & padj < pthreshold, padj, 1))
up.p2 <- with(df.cle, ifelse(log2FoldChange > 0 & padj < pthreshold, padj, 1))
names(up.p1) <- names(up.p2) <- row.names(df.bcn)
up.p <- merge(as.data.frame(up.p1), as.data.frame(up.p2), by=0)
up.p$p <- pmax(up.p$up.p1, up.p$up.p2)
p <- up.p$p
names(p) <- up.p$Row.names

go.res$BP.up <- goTair(p, whichGO = "BP", topNodes = 150, stat = c("Fisher"), 
                 out.path = "./output/Bill.Analysis/topGO", fn.prefix = "EIS.WT-BCN.vs.WT-CLE2.up")
go.res$MF.up <- goTair(p, whichGO = "MF", topNodes = 150, stat = c("Fisher"), 
                 out.path = "./output/Bill.Analysis/topGO", fn.prefix = "EIS.WT-BCN.vs.WT-CLE2.up")
go.res$CC.up <- goTair(p, whichGO = "CC", topNodes = 150, stat = c("Fisher"), 
                 out.path = "./output/Bill.Analysis/topGO", fn.prefix = "EIS.WT-BCN.vs.WT-CLE2.up")

dn.p1 <- with(df.bcn, ifelse(log2FoldChange < 0 & padj < pthreshold, padj, 1))
dn.p2 <- with(df.cle, ifelse(log2FoldChange < 0 & padj < pthreshold, padj, 1))
names(dn.p1) <- names(dn.p2) <- row.names(df.bcn)
dn.p <- merge(as.data.frame(dn.p1), as.data.frame(dn.p2), by=0)
dn.p$p <- pmax(dn.p$dn.p1, dn.p$dn.p2)
p <- dn.p$p
names(p) <- dn.p$Row.names

go.res$BP.dn <- goTair(p, whichGO = "BP", topNodes = 150, stat = c("Fisher"), 
                 out.path = "./output/Bill.Analysis/topGO", fn.prefix = "EIS.WT-BCN.vs.WT-CLE2.dn")
go.res$MF.dn <- goTair(p, whichGO = "MF", topNodes = 150, stat = c("Fisher"), 
                 out.path = "./output/Bill.Analysis/topGO", fn.prefix = "EIS.WT-BCN.vs.WT-CLE2.dn")
go.res$CC.dn <- goTair(p, whichGO = "CC", topNodes = 150, stat = c("Fisher"), 
                 out.path = "./output/Bill.Analysis/topGO", fn.prefix = "EIS.WT-BCN.vs.WT-CLE2.dn")

for (n in names(go.res)) {
  d <- go.res[[n]]
  prefix <- paste0("EIS.WT_BCN.vs.WT_CLE2.", n)
  if(grepl("\\.dn", n)) {
    goVisBar(d, col = "skyblue3", img.path = "./output/Bill.Analysis/topGO", fn.prefix = prefix)
  } else {
    goVisBar(d, img.path = "./output/Bill.Analysis/topGO", fn.prefix = prefix)
  }
  
}



## genelist: EIS, wt-bcn & wt-CLE ####
clvBCN <- EIS.res$clvBCN
clvCLE <- EIS.res$clvCLE

pos.bcn <- filter(clvBCN, log2FoldChange > 0, padj < 0.05)
pos.cle <- filter(clvCLE, log2FoldChange > 0, padj < 0.05)
pos.ov <- calculate.overlap(list(BCN = pos.bcn$id, CLE=pos.cle$id))
pos.ov.genes <- filter(geneInfo, Gene_id %in% pos.ov$a3)
write.table(pos.ov.genes, "./output/Bill.Analysis/List.genes/WT_BCN.vs.WT_CLE.up.overlap.csv", row.names=F, col.names=T, sep=",", quote=F)

neg.bcn <- filter(clvBCN, log2FoldChange < 0, padj < 0.05)
neg.cle <- filter(clvCLE, log2FoldChange < 0, padj < 0.05)
neg.ov <- calculate.overlap(list(BCN = neg.bcn$id, CLE=neg.cle$id))
neg.ov.genes <- filter(geneInfo, Gene_id %in% neg.ov$a3)
write.table(neg.ov.genes, "./output/Bill.Analysis/List.genes/WT_BCN.vs.WT_CLE.negative.overlap.csv", row.names=F, col.names=T, sep=",", quote=F)

## genelist: EIS&LCM, clv_BCN.vs.WT_BCN
EIS.vwb <- EIS.res$vwb
lcm.vwb <- LCM.res$vwb

up.EIS.vwb <- filter(EIS.vwb, log2FoldChange > 0, padj < 0.05)
up.lcm.vwb <- filter(lcm.vwb, log2FoldChange > 0, padj < 0.05)
up.vwb <- calculate.overlap(list(EIS = up.EIS.vwb$id, LCM=up.lcm.vwb$id))
up.vwb.genes <- filter(geneInfo, Gene_id %in% up.vwb$a3)
write.table(up.vwb.genes, "vwb.up.genes.csv", sep = ",", row.names = F, quote=F)

dn.EIS.vwb <- filter(EIS.vwb, log2FoldChange < 0, padj < 0.05)
dn.lcm.vwb <- filter(lcm.vwb, log2FoldChange < 0, padj < 0.05)
dn.vwb <- calculate.overlap(list(EIS = dn.EIS.vwb$id, LCM=dn.lcm.vwb$id))
dn.vwb.genes <- filter(geneInfo, Gene_id %in% dn.vwb$a3)
write.table(dn.vwb.genes, "vwb.dn.genes.csv", sep = ",", row.names = F, quote=F)


## some other code ##########

GOobj <- new("topGOdata",
             ontology = "BP",
             allGenes = lcm.wbc.up,
             #allScore = 0.05,
             geneSel = function(x) x < 0.05 ,
             nodeSize = 10, 
             annot = annFUN.org,
             mapping = "org.At.tair.db")


method <- c("weight01", "classic", "elim")
test <- c("Fisher", "KS")
par <- expand.grid(method, test, stringsAsFactors = F)

result <- list()
GOtable <- list()
for (i in 1:nrow(par)) {
  name <- paste0(par[i,], collapse = "")
  res <- runTest(GOobj, algorithm = par[i,1], statistic = par[i,2])
  
  result[[name]] <- res
  # printGraph(GOobj, res, firstSigNodes = nrow(gtable), fn.prefix = paste0("LCM.dn.GO.MF_", name), pdfSW=T, useInfo="all")
  
}

## box-plot of pvalues from different method
ps <- lapply(result, score)
ps <- lapply(ps, function(x) x[order(names(x))])
ps.df <- as.data.frame(ps)
ps.df <- ps.df[order(row.names(ps.df)),]

# pvalue
boxplot(ps.df, las=2, ylab = "p-value")
stripchart(ps.df, vertical=T, add=T, method = "jitter", jitter = 0.3, pch=20, col = "lightgreen", cex = 1.5)
# log10(pvalue)
boxplot(log10(ps.df), las=2, ylab = "log10(p-value)")
stripchart(log10(ps.df), vertical=T, add=T, method = "jitter", jitter = 0.3, pch=20, col = "lightgreen", cex = 1.5)

## save all result to a csv file.
goInfo <- GenTable(GOobj, p = result[[1]], topNodes = length(score(res)))
goInfo <- goInfo[order(goInfo$GO.ID),1:(ncol(goInfo)-1)]
res <- cbind(goInfo, ps.df)
res$Term <- gsub(",", ";", res$Term)
res <- res[order(res$weight01Fisher),]
write.table(res, "LCM.BCN.up.MF_compare.csv", col.names = T, row.names = F, sep = ",", quote = F)

  
