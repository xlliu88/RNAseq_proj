
## import files
source_files <- c("_loadPackages.r", 
                  "_DESeq2Processor.r", 
                  "_plotters.r", 
                  "_importResult.bill.r",
                  "_bill.vars.r", 
                  "_celltype.vars.r", 
                  "_celltype.fun.r", 
                  "_go_kegg.r")

sapply(file.path("./scripts", source_files), source)


# out.path0 <- file.path(proj.path, "Enrd.SVA3.[gcrma_fc1.2p1e03]")
# assign("p.thd", 0.001, envir = .GlobalEnv)
# assign("fc.thd", 1.2, envir = .GlobalEnv)
# assign("mthd", "gcrma", envir = .GlobalEnv)
# prefix <- sprintf("Enrd.SVA3.[%s_fc%1.1f-p%s]", mthd, fc.thd, formatC(p.thd, format = "e", digits = 0))
# prefix <- gsub("-", "", prefix)
# out.path <- file.path(proj.path, prefix)
proj.path <- "rootcells2007"

p <-  0.001
fc <- c(1.0, 1.2)# 1.5, 2.0)
sample <- c("WT", "clv")
suffix <- "[Mock_after_sva]"

cutoffs <- expand.grid(fc, p) %>% 
  rename(Var1 = "fc", Var2 = "p")

## ATHB15 vs EIS vs LCM, overlap of WT and clv

for (i in 1:nrow(cutoffs)) {
    path2 <- file.path(proj.path,
                       sprintf("Enrd.SVA3.[gcrma_fc%1.1fp%s]", 
                               cutoffs$fc[i], 
                               formatC(cutoffs$p[i], format = "e", digits = 0))) %>% 
      str_replace(., "-", "")
    path.wt <- sprintf("Cty2WTsyn%s", suffix)
    path.clv <- sprintf("Cty2clvsyn%s", suffix)
    file.wt <- "[Gene_Details].ATHB15-WTsyn_Overlapped.csv"
    file.clv <- "[Gene_Details].ATHB15-clvsyn_Overlapped.csv"
    
    genelist.wt <- read_csv(file.path(path2, path.wt, file.wt))
    genes.wt <- genelist.wt %>% 
      select(Gene_id) %>% 
      unlist() %>% 
      unname()
    
    genes.clv <- read_csv(file.path(path2, path.clv, file.clv)) %>% 
      select(Gene_id) %>% 
      unlist() %>% 
      unname()
    
    vlist <- list(WT = genes.wt, clv = genes.clv)
    v.cols <- brewer.pal(3, "Set2")[1:2] #palette.colors()[2:3]
    pdf(file.path(path2,path.clv,"[clv-WT].ATHB15-EIS-LCM_overlapped.pdf"), width = 3, height = 3)
    plot(euler(vlist), fill = colRamp(v.cols, 1.5), edge = v.cols, label = T, quantities = T)
    dev.off()
    
    g <- genelist.wt %>% 
      filter(Gene_id %in% intersect(genes.wt, genes.clv)) %>% 
      select(-`...1`)
    write_csv(g, file.path(path2, path.clv, "[clv-WT][Gene_Details].ATHB-EIS-LCM_overlapped.csv"))
    
}
## Go enrichment of ATHB15-EIS-LCM overlapped genes
mk <- "ATHB15"
for (i in 1:nrow(cutoffs)){
  for (gt in c("WT", "clv")) {
    path2 <- file.path(proj.path,
                       sprintf("Enrd.SVA3.[gcrma_fc%1.1fp%s]", 
                               cutoffs$fc[i], 
                               formatC(cutoffs$p[i], format = "e", digits = 0))) %>% 
      str_replace(., "-", "")
    path.gt <- sprintf("Cty2%ssyn%s", gt, suffix)
    file.name <- sprintf("[Gene_Details].%s-%ssyn_Overlapped.csv", mk, gt)
    
    overlap_genes <- read_csv(file.path(path2, path.gt, file.name))
    genelist <- overlap_genes %>% 
      select(Gene_id) %>% 
      unlist() %>% 
      unname()
    
    x <- rep(1, length(gene_universe))
    names(x) <- gene_universe
    x[names(x) %in% genelist] <- 0.001
    
    for (goterms in c("BP", "CC", "MF")) {
      outfn.prefix <- sprintf("[GO][%s&EIS&LCM]", mk)
      go <- goTair(x, goterms, 
                   algorithm = "weight01",
                   stat = "KS",
                   out.path = file.path(path2, path.gt), 
                   fn.prefix = outfn.prefix) 
      
      goVisBar(go, img.path = file.path(path2, path.gt), pthreshold = fisher.p.thd, 
               fn.prefix = sprintf("%s-%ssyn.%s", mk, gt, goterms))
      
    }

  }
}
   
## down regulated genes
x <- EIS.wbc$padj
names(x) <- EIS.wbc$Gene_id
x[which(EIS.wbc$log2FoldChange > 0)] <- 1


outfn.prefix <- "[GO][EIS.down]"
go <- goTair(x, goterms, 
             algorithm = "weight01",
             stat = c("fisher", "KS"),
             out.path = file.path(path2, path.gt), 
             fn.prefix = outfn.prefix) 

goVisBar(go, img.path = file.path(path2, path.gt), pthreshold = fisher.p.thd, 
         fn.prefix = sprintf("%s-%ssyn.%s", mk, gt, goterms))
