source("./scripts/_loadPackages.r")
source("./scripts/_DESeq2Processor.r")
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")
source("./scripts/_go_kegg.r")
library(RColorBrewer)
library(writexl)

p.thd <- 0.05
fc.thd <- 1.0
lfc.thd <- log2(fc.thd)

fn.prefix <- sprintf("[fc%.1fp%s]", fc.thd, formatC(p.thd, format = "g", digits = 0))
### import results
EIS.res <- importResults(bill.resfile.EISd75)
LCM.res <- importResults(bill.resfile.LCM)

### DEG analysis
alt.path <- "./output/Bill.Analysis/Bill_Alt"
if(!file.exists(alt.path)) dir.create(alt.path)
go.path <- file.path(alt.path, "analysis", "topGO")
if(!file.exists(go.path)) dir.create(go.path)
# altCLE
# up-regulated in CLE2 treatment, 
# down-regulated in clv mutant
# log2FC smaller in clv mutant than WT upon CLE2 treatment


cr1 <- EIS.wec$padj < p.thd & EIS.wec$log2FoldChange > lfc.thd
cr2 <- EIS.vwc$padj < p.thd & EIS.vwc$log2FoldChange < -lfc.thd
cr3 <- near(EIS.vec$log2FoldChange, EIS.wec$log2FoldChange) | EIS.vec$log2FoldChange < EIS.wec$log2FoldChange

idx <- cr1 & cr2 & cr3  ## 336 genes

df <- EIS.res[["wbc"]][idx,c("Gene_id", "Gene_name", "Gene_description", "baseMean")]
df$Gene_name[df$Gene_name == ""] <- "NaN"
df$Gene_description <- gsub(",", ";", df$Gene_description)
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df.x <- EIS.res[[nm]][match(df$Gene_id, EIS.res[[nm]]$Gene_id), c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df <- cbind(df, df.x)
  
}

fn <- sprintf("CLE_Alt1_CLE.up_clv.dn.%s.csv", fn.prefix)
write.csv(df, file.path(alt.path, fn), quote = F, row.names = F)


## Alt2: add: clvBCN interaction in EIS
EIS.wec <- EIS.res[["wec"]]
EIS.vec <- EIS.res[["vec"]]
EIS.vwc <- EIS.res[["vwc"]]
EIS.wbc <- EIS.res[["wbc"]]
EIS.vbc <- EIS.res[["vbc"]]
EIS.cB <- EIS.res[["clvBCN"]]

cr1 <- EIS.wec$padj < p.thd & EIS.wec$log2FoldChange > lfc.thd ## up in WT-CLE2 treatment
cr2 <- EIS.vwc$padj < p.thd & EIS.vwc$log2FoldChange < -lfc.thd ## down in clv mutant control vs WT control
cr3 <- near(EIS.vec$log2FoldChange, EIS.wec$log2FoldChange) | EIS.vec$log2FoldChange < EIS.wec$log2FoldChange
cr4 <- EIS.cB$padj < p.thd # & EIS.cB$log2FoldChange > lfc.thd
#cr5 <- EIS.vbc$padj < p.thd & EIS.vbc$log2FoldChange > lfc.thd
cr3a <- EIS.vec$padj > p.thd
cr4a <- EIS.wbc$padj < p.thd & EIS.wbc$log2FoldChange > lfc.thd
# cr5 <- EIS.vbc$padj > p.thd
cr6 <- EIS.res[["clvCLE"]]$padj < p.thd
idx <- cr1 & cr2 & cr3 & cr4
idx <- cr1 & cr2 & cr3a & cr6
#idx <- cr6
df <- EIS.res[["wbc"]][idx,c("Gene_id", "Gene_name", "Gene_description", "baseMean")]
df$Gene_name[df$Gene_name == ""] <- "NaN"
df$Gene_description <- gsub(",", ";", df$Gene_description)
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df.x <- EIS.res[[nm]][match(df$Gene_id, EIS.res[[nm]]$Gene_id), c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df <- cbind(df, df.x)
  
}

fn <- sprintf("CLE_Alt3a-clvCLE_overlap.%s.csv", fn.prefix)
write.csv(df, file.path(alt.path, fn), quote = F, row.names = F)

## Alt2 - down: add: clvBCN interaction in EIS

cr1 <- EIS.wec$padj < p.thd & EIS.wec$log2FoldChange < -lfc.thd ## down in WT-CLE2 treatment
cr2 <- EIS.vwc$padj < p.thd & EIS.vwc$log2FoldChange > lfc.thd ## up in clv mutant control vs WT control
#cr3 <- near(EIS.vec$log2FoldChange, EIS.wec$log2FoldChange) | EIS.vec$log2FoldChange < EIS.wec$log2FoldChange
cr4 <- EIS.cB$padj < p.thd # & EIS.cB$log2FoldChange > lfc.thd

idx <- cr1 & cr2 & cr4

df <- EIS.res[["wbc"]][idx,c("Gene_id", "Gene_name", "Gene_description", "baseMean")]
df$Gene_name[df$Gene_name == ""] <- "NaN"
df$Gene_description <- gsub(",", ";", df$Gene_description)
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df.x <- EIS.res[[nm]][match(df$Gene_id, EIS.res[[nm]]$Gene_id), c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df <- cbind(df, df.x)
  
}

fn <- sprintf("CLE_Alt2.dn_clv.up_clvBCN.%s.csv", fn.prefix)
write.csv(df, file.path(alt.path, fn), quote = F, row.names = F)


## Alt3: add: up-regulated in both WT and clv in BCN treatment

cr1 <- EIS.wec$padj < p.thd & EIS.wec$log2FoldChange > lfc.thd ## up in WT-CLE2 treatment
cr2 <- EIS.vwc$padj < p.thd & EIS.vwc$log2FoldChange < -lfc.thd ## down in clv mutant control vs WT control
cr3 <- near(EIS.vec$log2FoldChange, EIS.wec$log2FoldChange) | EIS.vec$log2FoldChange < EIS.wec$log2FoldChange
cr4 <- EIS.wbc$padj < p.thd & EIS.wbc$log2FoldChange > lfc.thd
#cr5 <- EIS.vbc$padj < p.thd & EIS.vbc$log2FoldChange > lfc.thd

idx <- cr1 & cr2 & cr3 & cr4

df <- EIS.res[["wbc"]][idx,c("Gene_id", "Gene_name", "Gene_description", "baseMean")]
df$Gene_name[df$Gene_name == ""] <- "NaN"
df$Gene_description <- gsub(",", ";", df$Gene_description)
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df.x <- EIS.res[[nm]][match(df$Gene_id, EIS.res[[nm]]$Gene_id), c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df <- cbind(df, df.x)
  
}

fn <- sprintf("CLE_Alt3.up_clv.dn_BCN.up.%s.csv", fn.prefix)
write.csv(df, file.path(alt.path, fn), quote = F, row.names = F)

## Alt3. down
cr1 <- EIS.wec$padj < p.thd & EIS.wec$log2FoldChange < -lfc.thd ## down in WT-CLE2 treatment
cr2 <- EIS.vwc$padj < p.thd & EIS.vwc$log2FoldChange > lfc.thd ## up in clv mutant control vs WT control
cr3 <- near(EIS.vec$log2FoldChange, EIS.wec$log2FoldChange) | EIS.vec$log2FoldChange > EIS.wec$log2FoldChange
cr4 <- EIS.wbc$padj < p.thd & EIS.wbc$log2FoldChange < -lfc.thd
#cr5 <- EIS.vbc$padj < p.thd & EIS.vbc$log2FoldChange > lfc.thd

idx <- cr1 & cr2 & cr3 & cr4

df <- EIS.res[["wbc"]][idx,c("Gene_id", "Gene_name", "Gene_description", "baseMean")]
df$Gene_name[df$Gene_name == ""] <- "NaN"
df$Gene_description <- gsub(",", ";", df$Gene_description)
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df.x <- EIS.res[[nm]][match(df$Gene_id, EIS.res[[nm]]$Gene_id), c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df <- cbind(df, df.x)
  
}

fn <- sprintf("CLE_Alt3.CLE.down_clv.up_BCN.down.%s.csv", fn.prefix)
write.csv(df, file.path(alt.path, fn), quote = F, row.names = F)

## Alt4.1
## up in WT HsCLE2 treatment
## not-significant in clv HsCLE2 treatment
## down in clv mutant control vs WT control
## clvBCN significant

cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange > lfc.thd
cr2 <- EIS.res$vec$padj > p.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange < -lfc.thd
cr4 <- EIS.res$clvBCN$padj < p.thd
idx <- cr1 & cr2 & cr3 & cr4
info.cols <- c("Gene_id", "Gene_name", "Gene_description", "baseMean")

df.c <- EIS.res[[1]][idx, info.cols]
df.c <- df.c %>%
  mutate(Gene_name = ifelse(Gene_name == "", "NaN", Gene_name)) %>%
  mutate(Gene_description = gsub(",", ";", Gene_description)) 

res <- list()
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df <- EIS.res[[nm]][idx, ] 
  res[[pnm]] <- df
  
  df.x <- df[, c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df.c <- cbind(df.c, df.x)
  
}


res[["Consolidated"]] <- df.c
fn <- "Alt4.1_CLE-responsive.up_clvBCN.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))


## Alt4.2
## down in WT HsCLE2 treatment
## not-significant in clv HsCLE2 treatment
## up in clv mutant control vs WT control
## clvBCN significant
cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange < -lfc.thd
cr2 <- EIS.res$vec$padj > p.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange > lfc.thd
cr4 <- EIS.res$clvBCN$padj < p.thd
idx <- cr1 & cr2 & cr3 & cr4
info.cols <- c("Gene_id", "Gene_name", "Gene_description", "baseMean")

df.c <- EIS.res[[1]][idx, info.cols]
df.c <- df.c %>%
  mutate(Gene_name = ifelse(Gene_name == "", "NaN", Gene_name)) %>%
  mutate(Gene_description = gsub(",", ";", Gene_description)) 

res <- list()
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df <- EIS.res[[nm]][idx, ] 
  res[[pnm]] <- df
  
  df.x <- df[, c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df.c <- cbind(df.c, df.x)
  
}


res[["Consolidated"]] <- df.c
fn <- "Alt4.2_CLE-responsive.down_clvBCN.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))


############################################################################
## Alt 5.1 a was used in final analysis                                   ## 
############################################################################
##                                                                        ##
## Alt 5.1                                                                ##
## up in WT HsCLE2 treatment                                              ##
## not-significant in clv HsCLE2 treatment                                ##
## down in clv mutant control vs WT control                               ##
## up in WT BCN infection                                                 ##
##                                                                        ##
## 5.1a: clv_BCN vs WT_BCN up               --> 147 genes                 ##
## 5.1b: clv_BCN vs WT_BCN not_significant  --> 0 genes                   ##
## 5.1c: clv_BCN vs WT_BCN down             --> 41 genes                  ##
## 5.1d: clv_BCN vs WT_BCN down & clvBCN significant                      ##
##                                                                        ##
############################################################################

## when involving multiple (n) comparisons, the overall p-value would be
##  1-(1-pth.d)^n < 0.05
## set p.thd = 0.01; overall p value is: 0.029
## if use p_overall <- 1-(1-EIS.res$wec$padj) * (1-EIS.res$wbc$padj) * (1-EIS.res$vwc$padj)
##    gene list would be the same as use individual p < 0.05
p.thd2 <- 1-(1-p.thd)**(1/3) ## 0.017
p.thd3 <- 0.01
cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange > lfc.thd
cr2 <- EIS.res$vec$padj > p.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange < -lfc.thd
cr4 <- EIS.res$wbc$padj < p.thd & EIS.res$wbc$log2FoldChange > lfc.thd

cr5.1 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange < lfc.thd
cr5.2 <- EIS.res$clvBCN$padj >= p.thd 
cr5.3 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange > lfc.thd
cr6.1 <- EIS.res$vwb$padj < p.thd & EIS.res$vwb$log2FoldChange < lfc.thd
cr6.2 <- EIS.res$vwb$padj >= p.thd
cr6.3 <- EIS.res$vwb$padj < p.thd & EIS.res$vwb$log2FoldChange > lfc.thd

crx <- expand.grid(c("cr5.1", "cr5.2", "cr5.3"), 
                   c("cr6.1", "cr6.2", "cr6.3"),
                   stringsAsFactors = F)
crx_nms <- expand.grid(c("clvBCN.DOWN", "clvBCN.NS", "clvBCN.UP"), 
                       c("vwb.DOWN", "vwb.NS", "vwb.UP"))
#crx_nms <- crx_nms[, 2:1]
row.names(crx) <- apply(crx_nms, 1, paste0, collapse = "_")
#crx <- crx[, 2：1]
idx0 <- cr1 & cr2 & cr3 & cr4


tt <- 0
res <- list()
for(i in 1:nrow(crx)) {
  idx <- idx0 & get(crx[i,1]) & get(crx[i,2]) 
  cat(sprintf("%s:\t%d\n", row.names(crx)[i], sum(idx)))
  tt <- tt + sum(idx)
  
  if(sum(idx) == 0) next
  lab <- data.frame(label = rep(row.names(crx)[i], sum(idx)))
  dfs <- lapply(EIS.res, function(df) df[idx, ])
  dfs <- lapply(dfs, function(df) cbind(df, lab))
  if(length(res) == 0) { 
    res <- dfs
  } else {
    res <- lapply(1:length(res), function(x) rbind(res[[x]], dfs[[x]]))
  }
  
}
names(res) <- names(EIS.res)
res <- lapply(res, function(df) select(df, label, everything()))

info.cols <- c("label", "Gene_id", "Gene_name", "Gene_description", "baseMean")
df.c <- res[[1]][, info.cols]
for(i in 1:length(res)) {
  dfx <- res[[i]][, c("log2FoldChange", "padj")]
  colnames(dfx) <- c("log2FC", "padj")
  colnames(dfx) <- paste0(colnames(dfx), "(", names(res)[i], ")")
  df.c <- cbind(df.c, dfx)
}
df.c$Gene_name[which(df.c$Gene_name == "")] <- "NaN"
names(res) <- unlist(compare.names[names(res)])
names(res) <- paste(names(res), "(", names(EIS.res), ")", sep = "")
res[["Consolidated"]] <- df.c
fn <- "Alt5.1[p_overall0.03]_CLE_responsive_all_lower.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))

## go analysis of 5.1a list.

p <- ifelse(idx0, p.thd/10, 1)
names(p) <- EIS.res[[1]]$Gene_id
for (whichgo in c("BP", "CC", "MF")) {
  go.res <- goTair(p, 
                   whichGO = whichgo, 
                   topNodes = 20, 
                   stat = c("Fisher","KS"), 
                   out.path = go.path, 
                   fn.prefix = "alt5.1a")
  
  
}

## go visual plot
go.res <- readGO(go.path, "BP.*Fisher")
goVisBar(go.res[[1]], col = 'red3', img.path = go.path, fn.prefix = "goPlot.BP2")
goVisBar2(go.res, cols = "red3", img.path = go.path)
#########################################################################

## Alt 5.1b
## up in WT HsCLE2 treatment  (2550 genes)
## down in clv HsCLE2 treatment (0 genes)
## down in clv mutant control vs WT control
## up in WT BCN infection

### Alt 5.1c
## up in WT HsCLE2 treatment  (2550 genes)
## up in clv HsCLE2 treatment (790 genes)
## down in clv mutant control vs WT control (126 genes)
## up in WT BCN infection (55 genes) 

##########################################################################
cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange > lfc.thd
cr2 <- EIS.res$vec$padj < p.thd & EIS.res$vec$log2FoldChange > lfc.thd
cr2a <- EIS.res$vwe$padj < p.thd & EIS.res$vwe$log2FoldChange < -lfc.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange < -lfc.thd
cr4 <- EIS.res$wbc$padj < p.thd & EIS.res$wbc$log2FoldChange > lfc.thd

cr5.1 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange < lfc.thd
cr5.2 <- EIS.res$clvBCN$padj >= p.thd 
cr5.3 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange > lfc.thd

cr6.1 <- EIS.res$vwb$padj < p.thd & EIS.res$vwb$log2FoldChange < lfc.thd
cr6.2 <- EIS.res$vwb$padj >= p.thd
cr6.3 <- EIS.res$vwb$padj < p.thd & EIS.res$vwb$log2FoldChange > lfc.thd

crx <- expand.grid(c("cr5.1", "cr5.2", "cr5.3"), 
                   c("cr6.1", "cr6.2", "cr6.3"),
                   stringsAsFactors = F)
crx_nms <- expand.grid(c("clvBCN.DOWN", "clvBCN.NS", "clvBCN.UP"), 
                       c("vwb.DOWN", "vwb.NS", "vwb.UP"))
#crx_nms <- crx_nms[, 2:1]
row.names(crx) <- apply(crx_nms, 1, paste0, collapse = "_")
# crx <- crx[, 2：1]
idx0 <- cr1 & cr2  & cr2a & cr3 & cr4
cat(sprintf("Total genes: %d\n", sum(idx0)))

tt <- 0
res <- list()
for(i in 1:nrow(crx)) {
  idx <- idx0 & get(crx[i,1]) & get(crx[i,2]) 
  cat(sprintf("%s:\t%d\n", row.names(crx)[i], sum(idx)))
  tt <- tt + sum(idx)
  
  if(sum(idx) == 0) next
  lab <- data.frame(label = rep(row.names(crx)[i], sum(idx)))
  dfs <- lapply(EIS.res, function(df) df[idx, ])
  dfs <- lapply(dfs, function(df) cbind(df, lab))
  if(length(res) == 0) { 
    res <- dfs
  } else {
    res <- lapply(1:length(res), function(x) rbind(res[[x]], dfs[[x]]))
  }
  
}
names(res) <- names(EIS.res)
res <- lapply(res, function(df) select(df, label, everything()))

info.cols <- c("label", "Gene_id", "Gene_name", "Gene_description", "baseMean")
df.c <- res[[1]][, info.cols]
for(i in 1:length(res)) {
  dfx <- res[[i]][, c("log2FoldChange", "padj")]
  colnames(dfx) <- c("log2FC", "padj")
  colnames(dfx) <- paste0(colnames(dfx), "(", names(res)[i], ")")
  df.c <- cbind(df.c, dfx)
}
df.c$Gene_name[which(df.c$Gene_name == "")] <- "NaN"
names(res) <- unlist(compare.names[names(res)])
names(res) <- paste(names(res), "(", names(EIS.res), ")", sep = "")
res[["Consolidated"]] <- df.c
fn <- "Alt5.1b_CLE_responsive_all_lower.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))

## Alt 5.dn_a
## down in WT HsCLE2 treatment
## not-significant in clv HsCLE2 treatment
## up in clv mutant control vs WT control
## down in WT BCN infection

cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange < -lfc.thd
cr2 <- EIS.res$vec$padj > p.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange > lfc.thd
cr4 <- EIS.res$wbc$padj < p.thd & EIS.res$wbc$log2FoldChange < -lfc.thd

cr5.1 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange < lfc.thd
cr5.2 <- EIS.res$clvBCN$padj >= p.thd 
cr5.3 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange > lfc.thd

cr6.1 <- EIS.res$vwb$padj < p.thd & EIS.res$vwb$log2FoldChange < lfc.thd
cr6.2 <- EIS.res$vwb$padj >= p.thd
cr6.3 <- EIS.res$vwb$padj < p.thd & EIS.res$vwb$log2FoldChange > lfc.thd

crx <- expand.grid(c("cr5.1", "cr5.2", "cr5.3"), 
                   c("cr6.1", "cr6.2", "cr6.3"),
                   stringsAsFactors = F)
crx_nms <- expand.grid(c("clvBCN.DOWN", "clvBCN.NS", "clvBCN.UP"), 
                       c("vwb.DOWN", "vwb.NS", "vwb.UP"))
#crx_nms <- crx_nms[, 2:1]
row.names(crx) <- apply(crx_nms, 1, paste0, collapse = "_")
#crx <- crx[, 2：1]
idx0 <- cr1 & cr2 & cr3 & cr4


tt <- 0
res <- list()
for(i in 1:nrow(crx)) {
  idx <- idx0 & get(crx[i,1]) & get(crx[i,2]) 
  cat(sprintf("%s:\t%d\n", row.names(crx)[i], sum(idx)))
  tt <- tt + sum(idx)
  
  if(sum(idx) == 0) next
  lab <- data.frame(label = rep(row.names(crx)[i], sum(idx)))
  dfs <- lapply(EIS.res, function(df) df[idx, ])
  dfs <- lapply(dfs, function(df) cbind(df, lab))
  if(length(res) == 0) { 
    res <- dfs
  } else {
    res <- lapply(1:length(res), function(x) rbind(res[[x]], dfs[[x]]))
  }
  
}
names(res) <- names(EIS.res)
res <- lapply(res, function(df) select(df, label, everything()))

info.cols <- c("label", "Gene_id", "Gene_name", "Gene_description", "baseMean")
df.c <- res[[1]][, info.cols]
for(i in 1:length(res)) {
  dfx <- res[[i]][, c("log2FoldChange", "padj")]
  colnames(dfx) <- c("log2FC", "padj")
  colnames(dfx) <- paste0(colnames(dfx), "(", names(res)[i], ")")
  df.c <- cbind(df.c, dfx)
}
df.c$Gene_name[which(df.c$Gene_name == "")] <- "NaN"
names(res) <- unlist(compare.names[names(res)])
names(res) <- paste(names(res), "(", names(EIS.res), ")", sep = "")
res[["Consolidated"]] <- df.c
fn <- "Alt5.2a_CLE_responsive_all_higher.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))

## Alt 5.dn_b
## down in WT HsCLE2 treatment  
## up in clv HsCLE2 treatment   (4 genes)
## up in clv mutant control vs WT control (1 gene)
## down in WT BCN infection               (1 gene)

cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange < -lfc.thd
cr2 <- EIS.res$vec$padj < p.thd & EIS.res$vec$log2FoldChange > lfc.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange > lfc.thd
cr4 <- EIS.res$wbc$padj < p.thd & EIS.res$wbc$log2FoldChange < -lfc.thd

info.cols <- c("Gene_id", "Gene_name", "Gene_description", "baseMean")
idx0 <- cr1 & cr2 & cr3 & cr4
res <- lapply(EIS.res, function(df) df[idx0, info.cols])
names(res) <- unname(unlist(compare.names))
names(res) <- paste(names(res), "(", names(EIS.res), ")", sep = "")
fn <- "Alt5.2b_CLE_responsive_all_higher.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))


## Alt 5.dn_c
## down in WT HsCLE2 treatment  
## down in clv HsCLE2 treatment 
## up in clv vs WT in CLE2 treated samples
## up in clv mutant control vs WT control 
## down in WT BCN infection (145 genes)

cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange < -lfc.thd
cr2 <- EIS.res$vec$padj < p.thd & EIS.res$vec$log2FoldChange < -lfc.thd
cr2a <- EIS.res$vwe$padj < p.thd & EIS.res$vwe$log2FoldChange > lfc.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange > lfc.thd
cr4 <- EIS.res$wbc$padj < p.thd & EIS.res$wbc$log2FoldChange < -lfc.thd

cr5.1 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange < lfc.thd
cr5.2 <- EIS.res$clvBCN$padj >= p.thd 
cr5.3 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange > lfc.thd

cr6.1 <- EIS.res$vwb$padj < p.thd & EIS.res$vwb$log2FoldChange < lfc.thd
cr6.2 <- EIS.res$vwb$padj >= p.thd
cr6.3 <- EIS.res$vwb$padj < p.thd & EIS.res$vwb$log2FoldChange > lfc.thd

info.cols <- c("Gene_id", "Gene_name", "Gene_description", "baseMean")
idx0 <- cr1 & cr2 &cr2a & cr3 & cr4
res <- lapply(EIS.res, function(df) df[idx0, info.cols])
names(res) <- unname(unlist(compare.names))
names(res) <- paste(names(res), "(", names(EIS.res), ")", sep = "")


crx <- expand.grid(c("cr5.1", "cr5.2", "cr5.3"), 
                   c("cr6.1", "cr6.2", "cr6.3"),
                   stringsAsFactors = F)
crx_nms <- expand.grid(c("clvBCN.DOWN", "clvBCN.NS", "clvBCN.UP"), 
                       c("vwb.DOWN", "vwb.NS", "vwb.UP"))
row.names(crx) <- apply(crx_nms, 1, paste0, collapse = "_")
idx0 <- cr1 & cr2 & cr2a & cr3 & cr4

tt <- 0
res <- list()
for(i in 1:nrow(crx)) {
  idx <- idx0 & get(crx[i,1]) & get(crx[i,2]) 
  cat(sprintf("%s:\t%d\n", row.names(crx)[i], sum(idx)))
  tt <- tt + sum(idx)
  
  if(sum(idx) == 0) next
  lab <- data.frame(label = rep(row.names(crx)[i], sum(idx)))
  dfs <- lapply(EIS.res, function(df) df[idx, ])
  dfs <- lapply(dfs, function(df) cbind(df, lab))
  if(length(res) == 0) { 
    res <- dfs
  } else {
    res <- lapply(1:length(res), function(x) rbind(res[[x]], dfs[[x]]))
  }
  
}
names(res) <- names(EIS.res)
res <- lapply(res, function(df) select(df, label, everything()))

info.cols <- c("label", "Gene_id", "Gene_name", "Gene_description", "baseMean")
df.c <- res[[1]][, info.cols]
for(i in 1:length(res)) {
  dfx <- res[[i]][, c("log2FoldChange", "padj")]
  colnames(dfx) <- c("log2FC", "padj")
  colnames(dfx) <- paste0(colnames(dfx), "(", names(res)[i], ")")
  df.c <- cbind(df.c, dfx)
}
df.c$Gene_name[which(df.c$Gene_name == "")] <- "NaN"
names(res) <- unlist(compare.names[names(res)])
names(res) <- paste(names(res), "(", names(EIS.res), ")", sep = "")
res[["Consolidated"]] <- df.c
fn <- "Alt5.2c_CLE_responsive_all_higher.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))



## ==========================================================================================
## CLE responsive genes .vs. clvCLE interaction genes. 
## up in WT HsCLE2 treatment
## not-significant in clv HsCLE2 treatment
## down in clv mutant control vs WT control
## compare to clvCLE interaction list

cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange > lfc.thd
cr2 <- EIS.res$vec$padj > p.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange < -lfc.thd
idx.Cr <- cr1 & cr2 & cr3
idx.clvCLE <- EIS.res$clvCLE$padj < p.thd & EIS.res$clvCLE$log2FoldChange < -lfc.thd
idx <- c(which(idx.Cr), which(idx.clvCLE))
info.cols <- c("Gene_id", "Gene_name", "Gene_description", "baseMean")

Filter_grp <- c(rep("CLE_responsive", sum(idx.Cr)), rep("clvBCN_interaction", sum(idx.clvCLE)))

df.c <- EIS.res[[1]][idx, info.cols]
df.c <- df.c %>%
  mutate(Gene_name = ifelse(Gene_name == "", "NaN", Gene_name)) %>%
  mutate(Gene_description = gsub(",", ";", Gene_description)) %>%
  mutate(Filter = Filter_grp) %>%
  select(Filter, everything())

res <- list()
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df <- EIS.res[[nm]][idx, ] %>%
    mutate(Filter = Filter_grp) %>%
    select(Filter, everything())
  res[[pnm]] <- df
  
  df.x <- df[, c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df.c <- cbind(df.c, df.x)
  
}


res[["Consolidated"]] <- df.c
fn <- "Alt4_CLE-responsive_not_filter_for_BCN.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))


## the opposite regulated genes: 
## 1. WT CLE treatment down regulated
## 2. clv mutant CLE treatment not-significant in 
## 3. up-regulated in clv mutant

cr1 <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange < -lfc.thd
cr2 <- EIS.res$vec$padj > p.thd
cr3 <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange > lfc.thd
idx.Cr <- cr1 & cr2 & cr3
idx.clvCLE <- EIS.res$clvCLE$padj < p.thd & EIS.res$clvCLE$log2FoldChange > lfc.thd
idx <- c(which(idx.Cr), which(idx.clvCLE))
info.cols <- c("Gene_id", "Gene_name", "Gene_description", "baseMean")

Filter_grp <- c(rep("CLE_responsive", sum(idx.Cr)), rep("clvBCN_interaction", sum(idx.clvCLE)))

df.c <- EIS.res[[1]][idx, info.cols]
df.c <- df.c %>%
  mutate(Gene_name = ifelse(Gene_name == "", "NaN", Gene_name)) %>%
  mutate(Gene_description = gsub(",", ";", Gene_description)) %>%
  mutate(Filter = Filter_grp) %>%
  select(Filter, everything())

res <- list()
for (nm in names(EIS.res)) {
  pnm <- compare.names[[nm]]
  df <- EIS.res[[nm]][idx, ] %>%
    mutate(Filter = Filter_grp) %>%
    select(Filter, everything())
  res[[pnm]] <- df
  
  df.x <- df[, c("log2FoldChange", "padj")]
  colnames(df.x)[1] <- "log2FC"
  colnames(df.x) <- paste(nm, colnames(df.x), sep = "_")
  df.c <- cbind(df.c, df.x)
  
}


res[["Consolidated"]] <- df.c
fn <- "Alt4_CLE-responsive(dn)_not_filter_for_BCN.xlsx"
writexl::write_xlsx(res, path = file.path(alt.path, fn))


cr1d <- EIS.res$wec$padj < p.thd & EIS.res$wec$log2FoldChange < -lfc.thd
cr2d <- EIS.res$vec$padj > p.thd
cr3d <- EIS.res$vwc$padj < p.thd & EIS.res$vwc$log2FoldChange > lfc.thd
idx.d <- cr1d & cr2d & cr3d


EIS.clvBCN.up <- with(EIS.res$clvBCN, log2FoldChange > lfc.thd & padj < p.thd)
EIS.clvBCN.dn <- with(EIS.res$clvBCN, log2FoldChange < -lfc.thd & padj < p.thd)
LCM.clvBCN.up <- with(LCM.res$clvBCN, log2FoldChange > lfc.thd & padj < p.thd)
LCM.clvBCN.dn <- with(LCM.res$clvBCN, log2FoldChange < -lfc.thd & padj < p.thd)

vlist.up <- list(CLEresp = EIS.res$wec$Gene_id[idx],
                 EIS.clvBCN = EIS.res$wbc$Gene_id[EIS.clvBCN.up],
                 LCM.clvBCN = LCM.res$wbc$Gene_id[LCM.clvBCN.up])
venn.diagram(vlist.up, 
             filename = file.path(alt.path, "clvBCN.up.png"), 
             fill = venn.col[1:length(vlist.up)])

vlist.dn <- list(CLEresp = EIS.res$wbc$Gene_id[idx.d],
                 EIS.clvBCN = EIS.res$wbc$Gene_id[EIS.clvBCN.dn],
                 LCM.clvBCN = LCM.res$wbc$Gene_id[LCM.clvBCN.dn])
venn.diagram(vlist.dn, 
             filename = file.path(alt.path, "clvBCN.dn.png"), 
             fill = venn.col[1:length(vlist.dn)])




gid <- EIS.res$wbc$Gene_id
for (nm in names(EIS.res)) {
  x <- EIS.res[[nm]][["Gene_id"]]
  y <- all(x == gid)
  cat(nm, "\t", y, "\n")
}


## genes that overlapped in clvBCN and clvCLE

cr1 <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange > 0
cr2 <- EIS.res$clvCLE$padj < p.thd & EIS.res$clvCLE$log2FoldChange > 0
idx <- cr1 & cr2

genes <- EIS.res[[1]][idx, info.cols]

cr1d <- EIS.res$clvBCN$padj < p.thd & EIS.res$clvBCN$log2FoldChange < 0
cr2d <- EIS.res$clvCLE$padj < p.thd & EIS.res$clvCLE$log2FoldChange < 0
idx <- cr1d & cr2d 

genes2 <- EIS.res[[1]][idx, info.cols]

res.interaction <- list(positive = genes,
                        negative = genes2)

writexl::write_xlsx(res.interaction, file.path(alt.path, "EIS_clvBCN_clvCLE_overlapped.xlsx"))


## controls

ref.ids0 <- lapply(EIS.res, function(df) df$Gene_id[df$padj > 0.05 & df$baseMean > 5000])
ref.ids <- calculate.overlap2(ref.ids0)
length(ref.ids)
#cat(ref.ids)
ids.idx <- EIS.res[[1]]$Gene_id %in% ref.ids
ids.geneInfo <- EIS.res[[1]][ids.idx, c("Gene_id", "Gene_name", "baseMean", "Gene_description")]
ids.geneInfo$Gene_description <- sub(",", ";", ids.geneInfo$Gene_description)
ids.geneInfo$Gene_description <- sub("[Source:.*]$", "", ids.geneInfo$Gene_description)
write.csv(ids.geneInfo, file.path(alt.path, "ref_gene_candidate_expression.csv"), sep = ",", quote = F)
