
## import libraries
## common library
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(fdrtool)) # false discovery rate tools
suppressMessages(library(qvalue))
suppressMessages(library(sva))
suppressMessages(library(pdftools))

## parallel computing
suppressMessages(library(parallel))
suppressMessages(library(doParallel))

## graphic
# suppressMessages(library(ggplot2))  # included in tidyverse
suppressMessages(library(RColorBrewer))
suppressMessages(library(rafalib))
suppressMessages(library(ggrepel))
suppressMessages(library(ggalt))
suppressMessages(library(pheatmap))
suppressMessages(library(gplots))
suppressMessages(library(VennDiagram))
suppressMessages(library(eulerr))
suppressMessages(library(geneplotter)) # multiecdf plotChr
suppressMessages(library(dendextend))

## high throughput data analysis
suppressMessages(library(DESeq2))
suppressMessages(library(limma))
#suppressMessages(library(simpleaffy))
suppressMessages(library(affyPLM))
suppressMessages(library(affyio))

## genome annotation
suppressMessages(library(AnnotationDbi))
suppressMessages(library(biomaRt))
suppressMessages(library(org.At.tair.db))
suppressMessages(library(topGO))
suppressMessages(library(KEGGREST))
suppressMessages(library(pathview))

select <- dplyr::select
anno.select <- AnnotationDbi::select
invisible(utils::memory.limit(16000))

# suppressMessages(library(genefilter))
# suppressMessages(library(LSD))			# a colorful plot tool
# suppressMessages(library(gage))
# suppressMessages((library(fmsb)))			# statistics datasets and functions; radarplot



# suppressMessages(library(EDASeq)) # Seq quality control
# 
# library(BiocStyle) # HTML, PDF output
# library(GOstats)