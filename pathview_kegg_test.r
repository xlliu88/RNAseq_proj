source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")
source("./scripts/_go_kegg.r")
select <- dplyr::select
anno.select <- AnnotationDbi::select

pw.path <- "./output/Bill.Analysis/pathview"
graph.path <- file.path(pw.path, "pdf")
dir.create(graph.path)
## KEGG Pathway Analysis ##################################################################
## list genes in each KEGG pathway
kegg.genes <- keggGeneSet()

EIS.res <- importResults(bill.resfile.EISd75)
LCM.res <- importResults(bill.resfile.LCM)

pathways.list <- keggList("pathway", "ath")
head(pathways.list)

pathway.codes <- sub("path:", "", names(pathways.list)) 

res <- keggTestFisher("ath00062", geneListUP, kegg.genes[["ath00062"]])

df.F <- keggTest(pathway.codes[1:6], geneListUP, kegg.genes, method = "Fisher")
#df.W <- keggTest(pathway.codes[1:6], geneListUP, method = "Wilcox")
fisher.path <- file.path(pw.path, "fisher")
wilcox.path <- file.path(pw.path, "wilcox")
if(!dir.exists(fisher.path)) dir.create(fisher.path)
if(!dir.exists(wilcox.path)) dir.create(wilcox.path)

for (source  in c("EIS", "LCM")) {
  resList <-get(sprintf("%s.res", source))
  
  for(drct in c("up", "down")) {
    for(name in names(resList)) {
      cat(sprintf("KEGG Testing:\t%s - %s - %s\n", source, name, drct))
      
      df <- resList[[name]]
      geneList <- df$padj
      if(drct == "up"){
        geneList <- ifelse(df$log2FoldChange < 0, 1, geneList)
      } else {
        geneList <- ifelse(df$log2FoldChange > 0, 1, geneList)
      }
      
      names(geneList) <- df$Gene_id
      fn.prefix <- sprintf("%s.%s.%s", source, name, drct)
      fn.Fisher <- sprintf("%s_fisher_test.csv", fn.prefix)
      fn.Wilcox <- sprintf("%s_wilcox_test.csv", fn.prefix)
      
      res.F <- keggTest(pathway.codes, geneList, gSet = kegg.genes, method = "Fisher")
      cat("writing Fisher test result...\n")
      write_csv(res.F, file.path(fisher.path, fn.Fisher))
      
      res.W <- keggTest(pathway.codes, geneList, gSet = kegg.genes, method = "Wilcox")
      cat("writing Wilcox test result...\n")
      write_csv(res.W, file.path(wilcox.path, fn.Wilcox))
    }
  }
}






df <- EIS.res$wbc
for (k in names(EIS.kegg.wbc)) {
  ids <- EIS.kegg.wbc[[k]]$KEGG.ID

  for (id in ids) {
    pathwayGraph(df = df, pathway.id = id, img.path = graph.path, out.prefix = k, 
                 same.layer = F, file.type = "png")
  }
  
}



df <- LCM.res$wbc
for (k in names(LCM.kegg.wbc)) {
  ids <- LCM.kegg.wbc[[k]]$KEGG.ID
  
  for (id in ids) {
    pathwayGraph(df = df, pathway.id = id, img.path = graph.path, out.prefix = k, 
                 same.layer = F, file.type = "png")
  }
}

png.path <- file.path(pw.path, "png")
fs <- list.files(png.path) %>% 
  gsub("wbc\\.up\\.", "", .) %>% 
  sub("EWR", "EIS", .) %>% 
  sub("\\.png$", "", .) %>% 
  sub("\\.", "_", .) %>% 
  str_split("_") %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame()
row.names(fs) <- NULL
colnames(fs) <- c("Source", "KEGG_id", "KEGG_name")

EIS.fs <- fs %>% filter(Source == "EIS") %>% 
  mutate(Page = seq_along(1:nrow(.))) %>% 
  write.csv(., file = file.path(graph.path, "EIS_pathview_page_index.csv"))

EIS.fs <- fs %>% filter(Source == "LCM") %>% 
  mutate(Page = seq_along(1:nrow(.))) %>% 
  write.csv(., file = file.path(graph.path, "LCM_pathview_page_index.csv"))

## genes in each up/down regulated pathways
pw.up.genes <- genesInKEGG(row.names(pw.up), species = "ath")

## heatmap of up/down regulated genes in a select pathway
pathwayGeneCluster(dds.LCM, row.names(pw.up)[19], expr = LCM.res$wbc, ntop = 100, direction = "down")





