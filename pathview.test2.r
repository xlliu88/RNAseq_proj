source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")
source("./scripts/_go_kegg.r")
anno.select <- AnnotationDbi::select

pw.path <- "./output/Bill.Analysis/pathview"
graph.path <- file.path(pw.path, "pdf")
dir.create(graph.path)
## KEGG Pathway Analysis ##################################################################
## list genes in each KEGG pathway
kegg.genes <- keggGeneSet()

### pathview graphs##########
EIS.kegg.wbc <- readPathwayRes(pw.path, patterns = c("EWR", "wbc"))
LCM.kegg.wbc <- readPathwayRes(pw.path, patterns = c("LCM", "wbc"))

EIS.res <- importResults(bill.resfile.EISd75)
LCM.res <- importResults(bill.resfile.LCM)

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





