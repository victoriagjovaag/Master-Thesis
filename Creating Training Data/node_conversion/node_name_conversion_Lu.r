# node name conversion method to conert from from HGNC.symbols to node names in the Lu network 


library(dplyr)

node_name_conversion <- function(df) {
  df$source <- ifelse(df$source == "AKT1", "AKT",
               ifelse(df$source == "MAP3K5", "ASK1",
               ifelse(df$source == "PTGS2", "COX2",
               ifelse(df$source == "DCX", "DC",
               ifelse(df$source == "PTGER2", "EP2",
               ifelse(df$source == "MAPK1", "ERK",
               ifelse(df$source == "IL6ST", "GP130",
               ifelse(df$source == "ALPI", "IAP",
               ifelse(df$source == "MAPK8", "JNK",
               ifelse(df$source == "MAP3K1", "MEKK1",
               ifelse(df$source == "CDKN1A", "P21",
               ifelse(df$source == "TP53", "P53",
               ifelse(df$source == "PIK3CA", "PI3K",
               ifelse(df$source == "ROS1", "ROS",
               ifelse(df$source == "MBTPS1", "S1P",
               ifelse(df$source == "DIABLO", "SMAC",
               ifelse(df$source == "CISH", "SOCS",
               ifelse(df$source == "TGFB1", "TGFB",
               ifelse(df$source == "NELFCD", "TH1",
               ifelse(df$source == "TNF", "TNFA", df$source))))))))))))))))))))
  return(df)
}