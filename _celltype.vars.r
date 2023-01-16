## set up variables for the celltype signagure project

## directories and files
proj.path <- "rootcells2007"
cel.path <- file.path(proj.path, "CEL")
out.path <- file.path(proj.path, "output")

tar.file <- file.path(proj.path,"targets.csv")
mkmap.file <- file.path(proj.path, "mark2type.csv")
query.file <- file.path(proj.path, "query.matrix.csv")
tpm.ewr.file <- file.path(proj.path, "TPM_EWR.csv")
tpm.lcm.file <- file.path(proj.path, "TPM_LCM.csv")
affy.array.file <- file.path(proj.path, "affy_ATH1_array_elements-2010-12-20.txt")
enrd.file.brady <- file.path(proj.path, "rootcell.ENRD.genes.brady2007.csv")
enrd.file.zhang <- file.path(proj.path, "rootcell.ENRD.genes.Zhang2019.csv")


damagedfiles <- c("GSM133972.CEL.gz", "GSM133990.CEL.gz", "GSM133994.CEL.gz", 
                  "GSM226522_S17_3.CEL.gz", "GSM133982.CEL.gz")
notCEL <- c("GSE6349")

## parameters
p.thd <- 0.001
fc.thd <- 1.0
fisher.p.thd <- 0.01

## import vars.
cel.targets <- limma::readTargets(tar.file, sep = ",") %>%
  filter(!GSE %in% notCEL, !File.name %in% damagedfiles) %>%
  mutate(Date = sub("T", " ", Date)) %>%
  mutate(Date = sapply(Date, function(x) unlist(str_split(x, " "))[1])) %>%
  mutate(Date = ifelse(grepl("-", Date), 
                       format(as.Date(Date, "%Y-%m-%d")),
                       format(as.Date(Date, "%m/%d/%Y")))) %>% 
  arrange(Marker)
 

query.mat <- as.matrix(read.csv(query.file, sep = ",", row.names = 1, header = T))
query.mat[is.na(query.mat)] <- 0

probe2gene <- read.delim(affy.array.file, header = T, stringsAsFactors = F, sep = "\t", quote = "")
mk2cty <- read.csv(mkmap.file, header = T, stringsAsFactors = F)

gene_univ_seq <- read_delim(file.path(bill.anno$path, bill.anno$annotation), 
                            delim = "\t", show_col_types = F) %>% 
  select(Gene_id) %>% unlist %>% unname
gene_univ_ma <- probe2gene %>% 
  filter(!locus == "no_match") %>% 
  select(locus) %>% unlist() %>% 
  str_split(., ";") %>% unlist %>% unname %>% unique

# genes only exist in micro-array probe sets are usually,
# psudogenes, transposable element, or discontinued
gene_universe <- intersect(gene_univ_ma, gene_univ_seq)


enrd.brady <- read.csv(enrd.file.brady, header = T, stringsAsFactors = F)
enrd.arr15 <- read.csv(enrd.file.zhang, header = T, stringsAsFactors = F)
enrd.brady <- merge(enrd.brady, enrd.arr15, by = 0, all = T)
enrd.brady <- enrd.brady[order(as.numeric(enrd.brady$Row.names)),]
enrd.brady <- select(enrd.brady, -Row.names)
enrd.brady$ARR15[is.na(enrd.brady$ARR15)] <- ""
enrd.brady <- lapply(enrd.brady, function(x) unlist(str_split(x, ";")))
enrd.brady <- lapply(enrd.brady, function(x) x[!x %in% c("", "no_match")])




