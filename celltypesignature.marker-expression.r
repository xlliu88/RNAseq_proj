## expression of marker genes
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)

source("./scripts/_celltype.vars.r")

proj_path <- "rootcells2007"
xpt_path <- "Enrd.SVA3.[gcrma_fc1.2p1e03]"
expr_files <- list(raw = "Expression_matrix_Mock_raw.csv",
                 gc = "Expression_matrix_Mock_gcrma.csv",
                 sva = "Expression_matrix_Mock_gcrma_sva.csv")

expr_mats <- lapply(expr_files, 
                    function(x) {
                      read.csv(file.path(proj_path, xpt_path, x), row.names = 1) %>% 
                        as.matrix()
                      })
names(expr_mats) <- names(expr_files)


mk_df <- cel.targets %>% 
  select(Marker, Marker.AGI) %>% 
  mutate(Marker = ifelse(grepl("ARRI", Marker), "ARR15", Marker)) %>%
  filter(!Marker%in% c("None", "A8")) %>%
  filter(!Marker.AGI == "Unknown") %>% 
  unique() %>% 
  separate(Marker.AGI, into = c("AGI1", "AGI2"), sep = "/", fill = "left") %>%
  gather(key = "N", value = "AGI", 2:3, na.rm = T) %>%
  mutate(probe = probe2gene$array_element_name[match(AGI, probe2gene$locus)]) %>%
  filter(!is.na(probe)) %>% 
  #mutate(Gene = sprintf("%s (%s)", AGI, Marker)) #%>% 
  select(probe, AGI, Marker)

mk_path <- file.path(proj_path, xpt_path, "marker_expression")
if(!file.exists(mk_path)) dir.create(mk_path)
for(i in 1:3) {
  f <- file.path(proj_path, xpt_path, expr_files[[i]])
  dataset <- names(expr_mats)[i]
  mk_tbl <- read_csv(f, col_types = cols()) %>% 
    rename(probe = "X1") %>%
    filter(probe %in% mk_df$probe) %>%
    left_join(mk_df, by = "probe") %>% 
    select(colnames(mk_df), everything()) %>% 
    gather(key = "sample", value = "y", 4:ncol(.)) %>% 
    separate(sample, into = c("sample", "r"), sep = "\\.") %>% 
    mutate(sample = ifelse(grepl("ARRI", sample), "ARR15", sample)) %>% 
    select(-r) %>% 
    mutate(color = ifelse(Marker == sample, "self", "gray50")) %>% 
    mutate(Gene = sprintf("%s (%s)", AGI, Marker)) %>% 
    mutate(Gene = as.factor(Gene))
  
  for(gene in levels(mk_tbl$Gene)) {
    ti <- sprintf("expression of %s in self-enriched dataset (red) - %s",gene, dataset)
    #gene <- levels(mk_tbl$Gene)
    p <- mk_tbl %>% 
      filter(Gene == gene) %>% 
            ggplot(aes(x = sample, y = 2^y, fill = color)) +
            geom_bar(stat = "summary",
                     fun.y = "mean") +
            geom_errorbar(stat = "summary") +
            #facet_wrap(. ~ Gene, scales = "free_y")  + 
            scale_fill_manual(values = c("#A9A9A9", "#FF8888")) + 
            labs(title = ti,
                 y = "gene expression",
                 x = "samples") +
            theme_bw() + 
            theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
                  legend.position = "none")
    fn <- sprintf("%s_Expression_%s.pdf", gene, dataset)
    ggsave(file.path(mk_path, fn), plot = p, width = 4, height = 3)
  }

}


################################################################################
  tar <- cel.targets %>% 
    filter(DuplicateOf == 0,
           !Marker %in% c("A8", "None"),
           Treatment == "Mock")
   qc.path <- file.path(proj.path, "Enrd.SVA3.[gcrma_fc1.2p1e03]", "qc2")
   dir.create(out.path)
   mthd <- "gcrma"
   suffix <- "[fc1.2p1e03]"
  
  mat <- expr_mats$sva
  ss <- svd(mat - rowMeans(mat))
  var_xpln.sva <- ss$d^2/sum(ss$d^2)
  
  ## source of variances  
  var_source <- c("Marker", "Date", "GSE", "Treatment", "Age", "Substrate", "Ecotype", "Lab")
  for (vs in var_source) {
     tar[[vs]] <- as.factor(tar[[vs]])
     lvls <- nlevels(tar[[vs]])
     if(lvls == 1) var_source <- var_source[!var_source == vs]
  }
  
  pltn <- length(var_source) + 1
  pltr <- ifelse(pltn > 3, 2, 1)
  pltc <- ifelse(pltr == 1, pltn, ceiling(pltn/pltr))
  for(i in 1:4) {
    fn <- sprintf("e2.PC%i_var_source_%s %s.pdf", i, mthd, suffix)
    pdf(file.path(qc.path, fn), width = pltc * 2.5, height = pltr * 2.5, pointsize = 8)
    par(mfrow = c(pltr, pltc), mar = c(5.1, 4.1, 1.1, 1.1))
    
    for(vs in var_source) {
      boxplot(split(ss$v[, i], tar[[vs]]), las = 2, outline = F, outpch = "")
      stripchart(split(ss$v[, i], tar[[vs]]),  method = "jitter", 
                 add = T, vertical = T, cex = 0.5, pch = 19, col = "gray30")
      title(vs, line = -1)
    }
    plot(var_xpln.sva[1:20] * 100, pch = 16,
         #main = ti,
         ylab = "% variance explained",
         xlab = "Principal component")
    title(sprintf("PC%i (%.1f%%)", i, ss$d[i]^2/sum(ss$d^2) * 100), line = -1)
    points(i, var_xpln.sva[i] * 100, pch = 16, col = "red")
    dev.off()
  }
  
  pltc <- ceiling((pltn-1)/pltr)
  xlab <- sprintf("PC1 (%.1f%% Var)", ss$d[1]^2/sum(ss$d^2) * 100)
  ylab <- sprintf("PC2 (%.1f%% Var)", ss$d[2]^2/sum(ss$d^2) * 100)
  fn <- sprintf("f2a.PCA(svd)_by_factors_%s.%s.pdf", mthd, suffix)
  pdf(file.path(qc.path, fn), width = 2.5 * pltc, height = 2.5 * pltr)
  par(mfrow = c(pltr,pltc), mar = c(4.1,2.1,1.1,1.1))
  for(vs in var_source) {
    fac <- factor(tar[[vs]])
    legcol <- ifelse(nlevels(fac) > 15, 2, 1)
    plot(ss$d[1]*ss$v[, 1], ss$d[2]*ss$v[, 2], pch = 16, 
         col = cols[fac], cex = 0.5, cex.axis = 0.7, xlab = "", ylab = "")
         #xlab = xlab, ylab = ylab)
    legend("topright", levels(fac), pch = 19, xjust = 1, cex = 0.5,
           col = cols[1:nlevels(fac)], ncol = legcol)
    title(main = sprintf("colored_by_%s", vs), 
         line = 0.2, font.main = 1, cex.main = 0.8)
    title(xlab = xlab, ylab = ylab, cex.lab = 0.7, line = 2)
  }
  dev.off()
# 