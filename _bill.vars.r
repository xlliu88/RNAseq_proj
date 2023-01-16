### constants
### file path and names

### count tables, condition list, and annotation file
bill.counts <- list(path="./output/Bill.Analysis/Salmoncounts",
                    EIScnts = "EIS_Counts_Per_Sample.txt",
                    EISconditions = "EIS.samples.txt",
                    EISd75cnts = "EISd75_Counts_Per_Sample.txt",
                    EISd75conditions = "EISd75.samples.txt",
                    LCMcnts = "LCM_Counts_Per_Sample.txt",
                    LCMconditions = "LCM.samples.txt")

bill.anno <- list(path="./genomes/gff3.ensemble",
                  gff3 = "Arabidopsis_thaliana.TAIR10.42.gff3",
                  annotation = "Arabidopsis_thaliana.TAIR10.42.annotation")

### bill's analysis results
compare.names <- list(wec = "WT_CLE2.vs.WT_Ctrl",
                      wbc = "WT_BCN.vs.WT_Ctrl",
                      
                      vec = "clv_CLE2.vs.clv_Ctrl",
                      vbc = "clv_BCN.vs.clv_Ctrl",
                      
                      vwc = "clv_Ctrl.vs.WT_Ctrl",
                      vwe = "clv_CLE2.vs.WT_CLE2",
                      vwb = "clv_BCN.vs.WT_BCN",
                      
                      clvCLE = "interaction_clvCLE2",
                      clvBCN = "interaction_clvBCN")

bill.resfile.EISd75 <- list(path = "./output/Bill.Analysis/EIS_results_all_genes",
                        
                        wec = "EIS_wt.cle_cont_results.txt", 
                        wbc = "EIS_wt.worm_cont_results.txt",
                                         
                        vec = "EIS_mut.cle_cont_results.txt", 
                        vbc = "EIS_mut.worm_cont_results.txt",
                         
                        vwc = "EIS_control.mutant_wt_results.txt", 
                        vwe = "EIS_cle.mutant_wt_results.txt",  
                        vwb = "EIS_worm.mutant_wt_results.txt", 
                                       
                        clvCLE = "EIS_CLE_interact_results.txt",
                        clvBCN = "EIS_Worm_interact_results.txt")

bill.resfile.LCM <- list(path="./output/Bill.Analysis/LCM_results_all_genes",
                    
						wbc = "LCM_wt.worm_cont_results.txt",
						vbc = "LCM_mut.worm_cont_results.txt",
						vwc = "LCM_cont.mut_wt_results.txt",
						vwb = "LCM_worm.mut_wt_results.txt",
						clvBCN = "LCM_interact_results.txt")

### results rebuilt from bill's count files
rebuild.EIS <- list(path = "./output/salmon/EIS_results_contrast",
                    wec = "EIS.WT_HsCLE2.vs.WT_Ctrl.csv",
          					wbc = "EIS.WT_Infected.vs.WT_Ctrl.csv",
                               
          					vec = "EIS.clv_HsCLE2.vs.clv_Ctrl.csv",
          					vbc = "EIS.clv_Infected.vs.clv_Ctrl.csv",
                               
          					vwc = "EIS.clv_Ctrl.vs.WT_Ctrl.csv",
          					vwe = "EIS.clv_HsCLE2.vs.WT_HsCLE2.csv",
          					vwb = "EIS.clv_Infected.vs.WT_Infected.csv",
                               
          					clvCLE = "EIS.interaction_clv.HsCLE2.csv",
          					clvBCN = "EIS.interaction_clv.Infected.csv")

rebuild.EISd75 <- list(path = "./output/Bill.Analysis/XL.Rebuild",
                    wec = "EISd75.WT_CLE.vs.WT_Con.csv",
                    wbc = "EISd75.WT_Wrm.vs.WT_Con.csv",
                    
                    vec = "EISd75.MU_CLE.vs.MU_Con.csv",
                    vbc = "EISd75.MU_Wrm.vs.MU_Con.csv",
                    
                    vwc = "EISd75.MU_Con.vs.WT_Con.csv",
                    vwe = "EISd75.MU_CLE.vs.WT_CLE.csv",
                    vwb = "EISd75.MU_Wrm.vs.WT_Wrm.csv",
                    
                    clvCLE = "EISd75.interaction_MU.CLE.csv",
                    clvBCN = "EISd75.interaction_MU.Wrm.csv")

rebuild.LCM <- list(path = "./output/salmon/LCm_results_contrast",
                    wbc = "LCM.WT_Infected.vs.WT_Ctrl.csv",
                    vbc = "LCM.clv_Infected.vs.clv_Ctrl.csv",
                    vwc = "LCM.clv_Ctrl.vs.WT_Ctrl.csv",
                    vwb = "LCM.clv_Infected.vs.WT_Infected.csv",
                    clvBCN = "LCM.interaction_clv.Infected.csv")
# rebuild.LCM <- list(path = "./output/Bill.Analysis/XL.Rebuild",
#                     wbc = "LCM.WT_Wrm.vs.WT_Con.csv",
#                     vbc = "LCM.MU_Wrm.vs.MU_Con.csv",
#                     vwc = "LCM.MU_Con.vs.WT_Con.csv",
#                     vwb = "LCM.MU_Wrm.vs.WT_Wrm.csv",
#                     clvBCN = "LCM.interaction_MU.Wrm.csv")


### rebuild dds
dds.EIS <- ddsFromFeatureCounts(file.path(bill.counts$path, bill.counts$EIScnts), 
                                file.path(bill.counts$path, bill.counts$EISconditions), 
                                source = "tximport", 
                                gffa3 = file.path(bill.anno$path, bill.anno$gff3), 
                                design = ~ Genotype + Treatment + Genotype:Treatment)

dds.EISd75 <- ddsFromFeatureCounts(file.path(bill.counts$path, bill.counts$EISd75cnts), 
                                   file.path(bill.counts$path, bill.counts$EISd75conditions),
                                   source = "tximport", 
                                   gffa3 = file.path(bill.anno$path, bill.anno$gff3), 
                                   design = ~ Genotype + Treatment + Genotype:Treatment)


dds.LCM <- ddsFromFeatureCounts(file.path(bill.counts$path, bill.counts$LCMcnts), 
                                file.path(bill.counts$path, bill.counts$LCMconditions),
                                source = "tximport", 
                                gffa3 = file.path(bill.anno$path, bill.anno$gff3), 
                                design = ~ Genotype + Treatment + Genotype:Treatment)


dds.BCN <- ddsSlice(dds.EIS, factorName = "Treatment", factorLevel = c("Ctrl", "Infected"))
dds.CLE <- ddsSlice(dds.EIS, factorName = "Treatment", factorLevel = c("Ctrl", "HsCLE2"))

geneInfo <- as.data.frame(mcols(dds.EISd75))
geneInfo$Gene_id <- row.names(geneInfo)
geneInfo <- select(geneInfo, Gene_id, Gene_name, Gene_description)

TPM.EIS <- read.csv(file.path("./rootcells2007", "TPM_EIS.csv"), row.names = 1, header = T, stringsAsFactors = F)
TPM.LCM <- read.csv(file.path("./rootcells2007", "TPM_LCM.csv"), row.names = 1, header = T, stringsAsFactors = F)

## Syncytium up-regulated genes
## cat("Comparing Gene Enrichment Results to RNAseq dataset...\n")
EIS.wbc <- importResults(bill.resfile.EISd75, "wbc")[["wbc"]]
EIS.vbc <- importResults(bill.resfile.EISd75, "vbc")[["vbc"]]

LCM.wbc <- importResults(bill.resfile.LCM, "wbc")[["wbc"]]
LCM.vbc <- importResults(bill.resfile.LCM, "vbc")[["vbc"]]