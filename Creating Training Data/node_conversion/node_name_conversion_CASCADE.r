# node name conversion method to conert from from HGNC.symbols to node names in CASCADE 2.0 network 
# using file: https://github.com/druglogics/influential-nodes/blob/master/data/features/20190130_nodes_to_HGNC.txt
# specifically designed for decopleR methods (change column "source" for other methods)

#lines commended out are nodes with two alternative conversion names.

node_name_conversion <- function(df) {
                    df$source <-
                    ifelse(df$source == "STAT1" | df$source == "STAT2", "ISGF3_c", 
                    ifelse(df$source == "TCF7" | df$source == "TCF7L2", "TCF7_f",df$source))
                    return(df)
}

library(dplyr)

node_name_conversion <- function(df) {
                     df$source <-
                        ifelse(df$source == "AKT1" | df$source == "AKT2" | df$source == "AKT3",  "AKT_f", 
                        ifelse(df$source == "SMAD4"| df$source == "SMAD3" | df$source == "ATF2" | df$source == "FOS" | df$source == "JUN" , "AP1_c",
                        ifelse(df$source == "CREBBP" | df$source == "EP300", "CBPp300_c",          
                        ifelse(df$source == "CFL1" | df$source == "CFL2",  "CFL_f",                    
                        ifelse(df$source == "CSNK1A1"| df$source == "CSNK1D" | df$source == "CSNK1E" , "CK1_f",         
                        ifelse(df$source == "DKK4" | df$source == "DKK3" | df$source == "DKK2" | df$source == "DKK1", "DKK_f",
                        ifelse(df$source == "DKK4" | df$source == "DKK3" | df$source == "DKK2" | df$source == "DKK1", "DKK_g",
                       #ifelse(df$source == "DUSP1", "DUSP1",
                        ifelse(df$source == "DUSP1", "DUSP1_g",
                        ifelse(df$source == "DVL3" | df$source == "DVL2" | df$source == "DVL1", "DVL_f",
                        ifelse(df$source == "MAPK1" | df$source == "MAPK3", "ERK_f",
                        ifelse(df$source == "FOXO1" | df$source == "FOXO3" | df$source == "FOXO4", "FOXO_f",
                        ifelse(df$source == "FZD10"| df$source == "FZD9" | df$source == "FZD8"| df$source == "FZD7" | df$source == "FZD6" | df$source == "FZD5"| df$source == "FZD4" | df$source == "FZD3"| df$source == "FZD2" | df$source == "FZD1" , "FZD_f",
                        ifelse(df$source == "GAB1" | df$source == "GAB2", "GAB_f",
                        ifelse(df$source == "GSK3B" | df$source == "GSK3A", "GSK3_f",
                        ifelse(df$source == "IL6R" | df$source == "IL6ST" | df$source == "LIFR", "ILR_f",   
                        ifelse(df$source == "STAT2" | df$source == "STAT1", "ISGF3_c",
                        ifelse(df$source == "JAK2" | df$source == "JAK1", "JAK_f",
                        ifelse(df$source ==	"MAPK8" | df$source == "MAPK9", "JNK_f",
                        ifelse(df$source == "LRP6" | df$source == "LRP5", "LRP_f",
                       #ifelse(df$source == "MDM2", "MDM2",
                        ifelse(df$source == "MDM2", "MDM2_g",
                        ifelse(df$source == "MAP2K2" | df$source == "MAP2K1", "MEK_f",
                        ifelse(df$source == "MMP9" | df$source == "MMP7" | df$source == "MMP2", "MMP_f",
                        ifelse(df$source == "RPS6KA5"| df$source == "RPS6KA4", "MSK_f",
                        ifelse(df$source == "MLST8" | df$source == "RPTOR" | df$source == "MTOR", "mTORC1_c",
                        ifelse(df$source == "MLST8" | df$source == "RICTOR" | df$source == "MTOR", "mTORC2_c",
                        ifelse(df$source == "NFKB1" | df$source == "NFKB2", "NFKB_f",
                        ifelse(df$source == "NA", "PIP3",
                       #ifelse(df$source == "PTEN", "PTEN",
                        ifelse(df$source == "PTEN", "PTEN_g",
                        ifelse(df$source == "RAC3" | df$source == "RAC2" | df$source == "RAC1", "RAC_f",
                        ifelse(df$source == "RAF1" | df$source == "BRAF" | df$source == "ARAF", "RAF_f",
                        ifelse(df$source == "REL" | df$source == "RELB" | df$source == "RELA", "REL_f",
                        ifelse(df$source == "RPS6KA3" | df$source == "RPS6KA2" | df$source == "RPS6KA1", "RSK_f",
                        #ifelse(df$source == "EGFR"| df$source == "ERBB2"| df$source == "ERBB3" | df$source == "ERBB4"| df$source == "IGF1R"| df$source == "PDGFRA"| df$source == "PDGFRB"| df$source == "CSF1R"| df$source == "KIT"| df$source == "FLT3"| df$source == "FGFR1"| df$source == "FGFR2"| df$source == "FGFR3"| df$source == "FGFR4"| df$source == "FLT1"| df$source == "FLT4"| df$source == "KDR"| df$source == "MET", "RTPK_f",
                        ifelse(df$source == "EGFR"| df$source == "ERBB2"| df$source == "ERBB3" | df$source == "ERBB4"| df$source == "IGF1R"| df$source == "PDGFRA"| df$source == "PDGFRB"| df$source == "CSF1R"| df$source == "KIT"| df$source == "FLT3"| df$source == "FGFR1"| df$source == "FGFR2"| df$source == "FGFR3"| df$source == "FGFR4"| df$source == "FLT1"| df$source == "FLT4"| df$source == "KDR"| df$source == "MET", "RTPK_g",
                        ifelse(df$source == "RPS6KB2" | df$source == "RPS6KB1", "S6K_f",
                       #ifelse(df$source == "SFRP1", "SFRP1",
                        ifelse(df$source == "SFRP1", "SFRP1_g",
                       #ifelse(df$source == "SMAD6", "SMAD6",
                        ifelse(df$source == "SMAD6", "SMAD6_g",
                       #ifelse(df$source == "SMAD7", "SMAD7",
                        ifelse(df$source == "SMAD7", "SMAD7_g",
                       #ifelse(df$source == "SOCS1", "SOCS1",
                        ifelse(df$source == "SOCS1", "SOCS1_g",
                        ifelse(df$source == "LEF1", "LEF",
                        ifelse(df$source == "TAB2" | df$source == "TAB1", "TAB_f",
                        ifelse(df$source == "TCF7" | df$source == "TCF7L1"| df$source == "TCF7L2", "TCF7_f",
                        ifelse(df$source == "TSC1" | df$source == "TSC2", "TSC_f", df$source))))))))))))))))))))))))))))))))))))))))))
                        return(df)
              
}