### count tables, condition list, and annotaiton file
counts <- list(path="./output/Bill.Analysis/Salmoncounts",
                    EWRcnts = "EWR_Counts_Per_Sample.txt",
                    EWRconditions = "EWR.samples.txt",
                    EWRd75cnts = "EWRd75_Counts_Per_Sample.txt",
                    EWRd75conditions = "EWRd75.samples.txt",
                    LCMcnts = "LCM_Counts_Per_Sample.txt",
                    LCMconditions = "LCM.samples.txt")

anno <- list(path="./genomes/gff3.ensemble",
                  gff3 = "Arabidopsis_thaliana.TAIR10.42.gff3",
                  annotation = "Arabidopsis_thaliana.TAIR10.42.annotation")

### bill's analysis results
compare.names <- list(wec="WT_CLE2.vs.WT_Ctrl",
                      wbc="WT_BCN.vs.WT_Ctrl",
                      
                      vec="clv_CLE2.vs.clv_Ctrl",
                      vbc="clv_BCN.vs.clv_Ctrl",
                      
                      vwc="clv_Ctrl.vs.WT_Ctrl",
                      vwe="clv_CLE2.vs.WT_CLE2",
                      vwb="clv_BCN.vs.WT_BCN",
                      
                      clvCLE = "interaction_clvCLE2",
                      clvBCN = "interaction_clvBCN")