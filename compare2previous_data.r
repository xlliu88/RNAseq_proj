# compare RNAseq dataset to previous published dataset:
# 
# 1. Szakasits, D., Heinen, P., Wieczorek, K., Hofmann, J., Wagner, F., Kreil, D.P., Sykacek, P., Grundler, F.M.W., and Bohlmann, H. (2009). 
#    The transcriptome of syncytia induced by the cyst nematode Heterodera schachtii in Arabidopsis roots. 
#    Plant J. 57: 771–784.
# 2. Barcala, M., García, A., Cabrera, J., Casson, S., Lindsey, K., Favery, B., García-Casado, G., Solano, R., Fenoll, C., and Escobar, C. (2010). 
#    Early transcriptomic events in microdissected Arabidopsis nematode-induced giant cells. 
#    Plant J. 61: 698–712.
# 3. Yamaguchi, Y.L. et al. (2017). 
#    Root-Knot and Cyst Nematodes Activate Procambium-Associated Genes in Arabidopsis Roots. 
#    Front. Plant Sci. 8: 1195.
#    

source('./scripts/_loadPackages.r')
source("./scripts/_plotters.r")
source("./scripts/_bill.vars.r")
source("./scripts/_importResult.bill.r")

proj.path <- "syncytia2009"
out.path <- file.path(proj.path, "SynSeq.vs.Syn2009")
p.thd <- 0.001
fc.thd <- 1.2

EWR.wbc <- importResults(bill.resfile.EWRd75, "wbc")[["wbc"]]
LCM.wbc <- importResults(bill.resfile.LCM, "wbc")[["wbc"]]

syn2009 <- read.csv(file.path(proj.path, "tableS1.csv"), header = T, stringsAsFactors = F)
head(syn2009)

## filter out genes overlapped in RNAseq and microarray datasets
ids <- intersect(EWR.wbc$id, LCM.wbc$id)
ids <- intersect(ids, syn2009$GeneID)

EWR.wbc.up <- filter(EWR.wbc, log2FoldChange >= log2(fc.thd), padj <= p.thd)
LCM.wbc.up <- filter(LCM.wbc, log2FoldChange >= log2(fc.thd), padj <= p.thd)
syn2009.up <- filter(syn2009, M >= log2(fc.thd), adj.q <= p.thd)

EWR.wbc.dn <- filter(EWR.wbc, log2FoldChange <= log2(fc.thd), padj <= p.thd)
LCM.wbc.dn <- filter(LCM.wbc, log2FoldChange <= log2(fc.thd), padj <= p.thd)
syn2009.dn <- filter(syn2009, M <= log2(fc.thd), adj.q <= p.thd)

uplist <- list(EIS = EWR.wbc.up$id,
               LCM = LCM.wbc.up$id,
               DS09 = syn2009.up$GeneID)
dnlist <- list(EIS = EWR.wbc.dn$id,
               LCM = LCM.wbc.dn$id,
               DS09 = syn2009.dn$GeneID)

venn.diagram(uplist, filename = NA)
