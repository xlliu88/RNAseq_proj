source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")
source("./scripts/_go_kegg.r")


######################################################################################################
###                                                                                                ###
###         This Script is used to analysis GO Enrichment                                          ###
###         it use tables generated by "GO.genTables.r"                                            ###
###                                                                                                ###
######################################################################################################

pthreshold <- 0.05
go.pthreshold <- 0.01
go.path <- "./output/Bill.Analysis/topGO"

######################################################################################################
### import GO result
LCM.go <- readGO(file.path(go.path, "Fisher.tables"), patterns = "LCM")
EWR.go <- readGO(go.path, "Fisher.tables", patterns = "EWR")
LnE.go <- readGO(go.path, "Fishser.tables", patterns = "LnE")         ## genes co-regulated in LCM and EWR samples


### barplot for significantly enriched go terms ######################################################
for (n in names(LCM.go)) {
   d <- LCM.go[[n]]
   cat(n,"\n")
   if(grepl("\\.dn\\.", n)) {
      goVisBar(d, col = "skyblue3", img.path = "./output/Bill.Analysis/topGO", fn.prefix = n)
   } else {
      goVisBar(d, img.path = "./output/Bill.Analysis/topGO", fn.prefix = n)
   }
}

for (n in names(EWR.go)) {
  d <- EWR.go[[n]]
  cat(n,"\n")
  if(grepl("\\.dn\\.", n)) {
    goVisBar(d, col = "skyblue3", img.path = "./output/Bill.Analysis/topGO", fn.prefix = n)
  } else {
    goVisBar(d, img.path = "./output/Bill.Analysis/topGO", fn.prefix = n)
  }
  
}

for (n in names(LnE.go)) {
  d <- LnE.go[[n]]
  cat(n,"\n")
  if(grepl("\\.dn\\.", n)) {
    goVisBar(d, col = "skyblue3", img.path = "./output/Bill.Analysis/topGO", fn.prefix = n)
  } else {
    goVisBar(d, img.path = "./output/Bill.Analysis/topGO", fn.prefix = n)
  }
  
}

### get number of Significant GOs ##################
go.summary <- data.frame(matrix(nrow = 0, ncol = 5))
LCM.go.sig <- lapply(LCM.go, function(x) sum(x$pvalue < go.pthreshold))
LCM.go.sig <- t(as.data.frame(LCM.go.sig))
ids <- sapply(row.names(LCM.go.sig), function(x) str_split(x, "\\."))
ids <- t(as.data.frame(ids))
all(row.names(ids) == row.names(LCM.go.sig))
LCM.go.sig <- cbind(ids, LCM.go.sig)
colnames(LCM.go.sig) <- c("Source", "Compare", "Direction", "GOtype", "No.Sig.GO")
go.summary <- rbind(go.summary, LCM.go.sig)

LCM.go.sig <- lapply(EWR.go, function(x) sum(x$pvalue < go.pthreshold))
LCM.go.sig <- t(as.data.frame(LCM.go.sig))
ids <- sapply(row.names(LCM.go.sig), function(x) str_split(x, "\\."))
ids <- t(as.data.frame(ids))
all(row.names(ids) == row.names(LCM.go.sig))
LCM.go.sig <- cbind(ids, LCM.go.sig)
colnames(LCM.go.sig) <- c("Source", "Compare", "Direction", "GOtype", "No.Sig.GO")
go.summary <- rbind(go.summary, LCM.go.sig)

LCM.go.sig <- lapply(LnE.go, function(x) sum(x$pvalue < go.pthreshold))
LCM.go.sig <- t(as.data.frame(LCM.go.sig))
ids <- sapply(row.names(LCM.go.sig), function(x) str_split(x, "\\."))
ids <- t(as.data.frame(ids))
all(row.names(ids) == row.names(LCM.go.sig))
LCM.go.sig <- cbind(ids, LCM.go.sig)
colnames(LCM.go.sig) <- c("Source", "Compare", "Direction", "GOtype", "No.Sig.GO")
go.summary <- rbind(go.summary, LCM.go.sig)

write.table(go.summary, file.path("./output/Bill.Analysis/topGO", "GO_enrichmary_summary(p0.01).csv"), 
                                  row.names = F, col.names = T, sep = ",", quote = F)



### combine EIS, LCM, and LnE ##############################################################################

goplot.path <- file.path(go.path, "Fisher.plots[EIS-LCM-LnE]")
if(!file.exists(goplot.path)) dir.create(goplot.path)
go_cols <- brewer.pal(11, "RdBu")
for(goterm in c("BP", "MF", "CC")) {
  for(d in c("up", "dn")) {
    for(comp in c("wbc", "vbc", "clvBCN", "vwc", "vwb")) {
      fn <- sprintf("GO_Enrichment_EIS_LCM_LnE.%s.%s.%s.pdf", comp, d, goterm)
      go.res <- readGO(file.path(go.path, "Fisher.tables"), patterns = c(comp, d, goterm))
      go.res <- dfIntersect(go.res, by.col = "GO.ID")
      
      fun <- ifelse(d == "up", "head", "tail")
      cols <- do.call(fun, list(go_cols, length(go.res)))
      if(d == "dn") cols <- rev(cols)
      goVisBar2(go.res, showThreshold = T, select.method = "anySig",
                cols = cols, 
                cat.names = c("EIS", "LCM", "EIS & LCM"), 
                img.path = goplot.path, #"./output/Bill.Analysis/topGO/Fisher.plot",
                img.name = fn)
      
    }
  }
}


		  
## subGO Terms  ###############################################################################################
key <- list()
key$hormones <- "hormone|auxin|indoleacetic.*acid|cytokinin|ethylene|Brassinosteroid|abscisic.*acid|gibberellin|salicylic.*acid|jasmonic.*acid|Strigolactone"
key$development <- "root|shoot|petal|sepal|floral|morphogenesis|development|polarity|pattern formation|cell fate|differentiatioin"
key$cellfate <- "meristem|stem cell|differentiation|specification|pattern formation|cell fate"
key$vascular <- "xylem|phloem|cambium|meristem|embryo"
key$defense <- "bacterium|bacterial|chitin|insect|fungus|fungal|virus|defense|resistance|immune|cell death|autophagy"
key$abiotic <- "response to|stimulus"
key$chromsome <- "DNA methylation|DNA recombination|histone|chromosom|recombination|DNA conformation"
key$transport <- "transport"
key$deve2 <- paste(key$development, key$cellfate, sep = "|")
key$expreg <- "translation|transcription|DNA methylation|RNA methylation|histone.*methylation|mRNA splicing|histone acetylation"

wbc.up.BP <- readGO(go.path, patterns = c("wbc", "up", "BP"))
wbc.up.BP <- dfListTrim(wbc.up.BP, "GO.ID")
wbc.dn.BP <- readGO(go.path, patterns = c("wbc", "dn", "BP"))
wbc.dn.BP <- dfListTrim(wbc.dn.BP, "GO.ID")

for (k in names(key)) {
  cat(k, "\n")
  cat("BP up regulated: \n")
  wbc.up.BP.x <- lapply(wbc.up.BP, function(df) filter(df, grepl(key[[k]], Term)))#, !grepl(key[["defense"]], Term), !grepl(key[["hormones"]], Term)))
  goVisBar2(wbc.up.BP.x, showThreshold = T, select.method = "anySig",
            cols = head(brewer.pal(11, "RdBu"),length(wbc.up.BP)), 
            cat.names = c("EWR", "LCM", "EWR & LCM"), 
            showLegend = T,
            img.path = "./output/Bill.Analysis/topGO",
            fn.prefix = k)
  
  cat("BP down regulated: \n")
  wbc.dn.BP.x <- lapply(wbc.dn.BP, function(df) filter(df, grepl(key[[k]], Term)))#, !grepl(key[["defense"]], Term), !grepl(key[["hormones"]], Term)))
  goVisBar2(wbc.dn.BP.x, showThreshold = T, select.method = "anySig",
            cols = rev(tail(brewer.pal(11, "RdBu"),length(wbc.up.BP))), 
            cat.names = c("EWR", "LCM", "EWR & LCM"), 
            showLegend = T,
            img.path = "./output/Bill.Analysis/topGO",
            fn.prefix = k)
}


### combine LCM and EWR - vwc: clv vs WT. Ctrl##############################################################################
vwc.up.BP <- readGO(go.path, patterns = c("vwc", "up", "BP"))
vwc.up.BP <- dfListTrim(vwc.up.BP, "GO.ID")
goVisBar2(vwc.up.BP, showThreshold = T,
          cols = head(brewer.pal(11, "RdBu"),length(vwc.up.BP)), 
          cat.names = c("EWR", "LCM", "EWR & LCM"), 
          select.method = "allSig",
          img.path = "./output/Bill.Analysis/topGO")

vwc.dn.BP <- readGO(go.path, patterns = c("vwc", "dn", "BP"))
vwc.dn.BP <- dfListTrim(vwc.dn.BP, "GO.ID")
goVisBar2(vwc.dn.BP, showThreshold = T,
          cols = rev(tail(brewer.pal(11, "RdBu"),length(vwc.up.BP))), 
          cat.names = c("EWR", "LCM", "EWR & LCM"), 
          select.method = "allSig",
          img.path = "./output/Bill.Analysis/topGO")

vwc.up.MF <- readGO(go.path, patterns = c("vwc", "up", "MF"))
vwc.up.MF <- dfListTrim(vwc.up.MF, "GO.ID")
goVisBar2(vwc.up.MF, showThreshold = T,
          cols = head(brewer.pal(11, "RdBu"),length(vwc.up.BP)), 
          cat.names = c("EWR", "LCM", "EWR & LCM"), 
          select.method = "allSig",
          img.path = "./output/Bill.Analysis/topGO")

vwc.dn.MF <- readGO(go.path, patterns = c("vwc", "dn", "MF"))
vwc.dn.MF <- dfListTrim(vwc.dn.MF, "GO.ID")
goVisBar2(vwc.dn.MF, showThreshold = T,
          cols = rev(tail(brewer.pal(11, "RdBu"),length(vwc.up.BP))), 
          cat.names = c("EWR", "LCM", "EWR & LCM"), 
          select.method = "allSig",
          img.path = "./output/Bill.Analysis/topGO")


vwc.up.CC <- readGO(go.path, patterns = c("vwc", "up", "CC"))
vwc.up.CC <- dfListTrim(vwc.up.CC, "GO.ID")
goVisBar2(vwc.up.CC, showThreshold = T,
          cols = head(brewer.pal(11, "RdBu"),length(vwc.up.BP)), 
          cat.names = c("EWR", "LCM", "EWR & LCM"), 
          select.method = "allSig",
          img.path = "./output/Bill.Analysis/topGO")

vwc.dn.CC <- readGO(go.path, patterns = c("vwc", "dn", "CC"))
vwc.dn.CC <- dfListTrim(vwc.dn.CC, "GO.ID")
goVisBar2(vwc.dn.CC, showThreshold = T,
          cols = rev(tail(brewer.pal(11, "RdBu"),length(vwc.up.BP))), 
          cat.names = c("EWR", "LCM", "EWR & LCM"), 
          select.method = "allSig",
          img.path = "./output/Bill.Analysis/topGO")
