## convert gtf to saf file
library(stringr)
library(dplyr)
select <- dplyr::select
p <- "../genomes/gtf"
gtf <- "Arabidopsis_thaliana.TAIR10.59.gtf"
saf <- "Arabidopsis_thaliana.TAIR10.59.saf"

gtf.mtx <- read.table(file.path(p,gtf), skip=5, header=F, sep = "\t",stringsAsFactors = F)
names(gtf.mtx) <- c("Chr", "Source","Feature", "Start", "End", "Score", "Strand","Frame", "Attribute")
geneid<- str_extract(gtf.mtx$Attribute, "gene_id\\s\\w*\\;")
geneid <- sub("gene_id\\s", "", geneid)
geneid <- sub(";", "", geneid)
gtf.mtx$GeneID <- geneid
saf.mtx <- gtf.mtx %>%
     filter(Feature == "exon")%>%
     select(GeneID, Chr, Start, End, Strand)

write.table(saf.mtx, file = paste0(p, saf), sep="\t", quote=F, row.names = F, col.names = T)
