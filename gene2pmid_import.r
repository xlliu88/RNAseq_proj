## to read and filter publications related to arabidopsis;
## 

path <- "./genomes"
file <- "gene2pubmed.gz"
pubmed.url <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"


pubmedfile <- gzfile(file.path(path, file), open = "r")

con <- gzcon(url(pubmed.url))
pubmedfile <- textConnection(readLines(con))

gene2pm <- read.table(pubmedfile, header = T, stringsAsFactors = F, 
                          sep = "\t", quote = "", comment.char = "")
colnames(gene2pubmed) <- c("tax_id", "gene_id", "pmid")
gene2pm.ath <- gene2pm[gene2pm$tax_id == 3702, ]


