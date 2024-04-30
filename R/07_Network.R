#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================




#===============================================================================
## SET UP ----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

PATH      <- 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC'
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq')
setwd(PATH)

library(igraph)
source("C:/Users/E15639P/Doctorat/scRNA-seq/Usefull_functions.R")

## Read files
Candidates.UP   <- read.table(paste0(PATH, '/Saves/Bulk-RNA/Candidates_UP.txt'))[,1]
Candidates.DOWN <- read.table(paste0(PATH, '/Saves/Bulk-RNA/Candidates_DOWN.txt'))[,1]


LEAP_WT <- read.table(paste0(PATH_SAVE, '/WT/LEAP/results_indexed_filtered.txt'), header = T)
TF_list <- read.table('C:/Users/E15639P/Data/Lists/masterTFlist.txt')[,1]



#===============================================================================
## ANALYSIS --------------------------------------------------------------------
#===============================================================================


## Level 1 : subset of interaction involving most probably Nfil3 as TF regulator
# Subset Nfil3 regulators
Level1 <- LEAP_WT[LEAP_WT$gene_row %in% 'Nfil3' & LEAP_WT$gene_col %in% TF_list,]

# For each TF, check if reciprocal interaction has not a better correlation
remove <- c()
for(g in Level1$gene_col){
  cor1 <- Level1$Correlation[Level1$gene_col %in% g]
  cor2 <- LEAP_WT$Correlation[LEAP_WT$gene_row %in% g & LEAP_WT$gene_col %in% 'Nfil3'] 
  if(length(cor2) != 0 && abs(cor2) > abs(cor1)){
    remove <- c(remove, g)
  }
}
Level1       <- Level1[Level1$gene_col %!in% remove,]
Level1$Level <- 1

# Subset to Bulk candidates
Candidates <- c(Candidates.DOWN, Candidates.UP)
Level1 <- Level1[Level1$gene_col %in% Candidates,]

# Remove co-expression interactions and apply correlation cutoff
Level1 <- Level1[Level1$Lag %!in% 0,]
Level1 <- Level1[abs(Level1$Correlation) > 0.3,]

#### OLD -----------------------------------------------------------------------
## Level 2 : get all interactions involving Level1 and get their best direction
#Level2 <- data.frame()
#for(g in Level1$gene_col){
#  print(g)
#  remove <- c()
#  # Subset interactions involving g as regulator
#  subreg <- LEAP_WT[LEAP_WT$gene_row %in% g & LEAP_WT$gene_col %in% TF_list,]
#  for(d in subreg$gene_col){
#    # Check if reciprocal interaction has not a better correlation
#    cor1 <- subreg$Correlation[subreg$gene_col %in% d]
#    cor2 <- LEAP_WT$Correlation[LEAP_WT$gene_row %in% d & LEAP_WT$gene_col %in% g]
#    if(length(cor2) != 0 && abs(cor2) > abs(cor1)){
#      remove <- c(remove, g)
#    }
#  }
#  subreg <- subreg[subreg$gene_col %!in% remove,]
#  Level2 <- rbind(Level2, subreg)
#}
#Level2$Level <- 0
#### ---------------------------------------------------------------------------


Level2 <- data.frame()
for(g in Level1$gene_col){
  print(g)
  remove <- c()
  # Subset interactions involving g as regulator
  subreg <- LEAP_WT[LEAP_WT$gene_row %in% g & LEAP_WT$gene_col %in% Level1$gene_col,]
  for(d in subreg$gene_col){
    # Check if reciprocal interaction has not a better correlation
    cor1 <- subreg$Correlation[subreg$gene_col %in% d]
    cor2 <- LEAP_WT$Correlation[LEAP_WT$gene_row %in% d & LEAP_WT$gene_col %in% g]
    if(length(cor2) != 0 && abs(cor2) > abs(cor1)){
      remove <- c(remove, d)
    }
  }
  #subreg <- subreg[subreg$gene_col %!in% remove,]
  Level2 <- rbind(Level2, subreg)
}
Level2$Level <- 2

# Remove co-expression interactions
Level2 <- Level2[Level2$Lag %!in% 0,]
Level2 <- Level2[abs(Level2$Correlation) > 0.3,]



## iGraph generation for Candidates only and saving to GML file
Network    <- rbind(Level1, Level2)
Genes      <- unique(c(Network$gene_row, Network$gene_col))


nodes <- data.frame(Name   = Genes,
                    TF     = ifelse(Genes %in% TF_list, 1, 0))
edges <- data.frame(From   = Network$gene_row,
                    To     = Network$gene_col,
                    Cor    = Network$Correlation,
                    AbsCor = abs(Network$Correlation),
                    Lag    = Network$Lag,
                    Neglag = 1.5-as.numeric(Network$Lag/max(Network$Lag)),
                    Level  = Network$Level) 

graph <- graph_from_data_frame(edges, directed = T, vertices = nodes)
plot(graph, layout = layout.fruchterman.reingold, vertex.size = 8,
     edge.arrow.size = 0.5, vertex.label = NA)
write_graph(graph, paste0(PATH_SAVE, '/Network.gml'), format='gml')




LEAP_WT[LEAP_WT$gene_row %in% 'Id2' & LEAP_WT$gene_col %in% 'Tox',]
LEAP_WT[LEAP_WT$gene_row %in% 'Tox' & LEAP_WT$gene_col %in% 'Id2',]


#===============================================================================
## NFIL3, ID2, TOX & GATA3 -----------------------------------------------------
#===============================================================================

Candidates <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_All.txt')[,1]

## Initialize with NFIL3
Level1 <- LEAP_WT[LEAP_WT$gene_row %in% 'Nfil3',]
Level1$Level <- 1

## Extend to other factors
Factors <- c('Id2', 'Tox', 'Gata3')

for(f in Factors){
  sub       <- LEAP_WT[LEAP_WT$gene_row %in% f,]
  sub$Level <- 1
  Level1    <- rbind(Level1, sub)
}

# Remove co-expression interactions and subset to our candidates
Level1 <- Level1[Level1$Lag %!in% 0,]
Level1 <- Level1[Level1$gene_col %in% Candidates,]
#Level1 <- Level1[abs(Level1$Correlation) > 0.3,]

## iGraph generation for Candidates only and saving to GML file
Network    <- Level1
Genes      <- unique(c(Network$gene_row, Network$gene_col))


nodes <- data.frame(Name   = Genes,
                    TF     = ifelse(Genes %in% TF_list, 1, 0))
edges <- data.frame(From   = Network$gene_row,
                    To     = Network$gene_col,
                    Cor    = Network$Correlation,
                    AbsCor = abs(Network$Correlation),
                    DirCor = ifelse(Network$Correlation > 0, 'Up', 'Down'),
                    Lag    = Network$Lag,
                    Neglag = 1.5-as.numeric(Network$Lag/max(Network$Lag)),
                    Level  = Network$Level) 

graph <- graph_from_data_frame(edges, directed = T, vertices = nodes)
plot(graph, layout = layout.fruchterman.reingold, vertex.size = 8,
     edge.arrow.size = 0.5, vertex.label = NA)
#write_graph(graph, paste0(PATH_SAVE, '/Network_4factors_updated.gml'), format='gml')

#write.table(Network, paste0(PATH_SAVE, '/WT/LEAP/TABLE_4Factors_Network.txt'),quote = F, row.names = F, col.names = T)




#===============================================================================
## ONLY NFIL3, ID2, TOX, GATA3 TCF7, ZBTB16 & RUNX3 ----------------------------
#===============================================================================

Candidates <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_All.txt')[,1]
Factors <- c('Id2', 'Tox', 'Gata3', 'Tcf7', 'Zbtb16', 'Runx3')

## Initialize with NFIL3
Level1 <- LEAP_WT[LEAP_WT$gene_row %in% 'Nfil3' & LEAP_WT$gene_col %in% Factors,]
Level1$Level <- 1

## Extend to other factors


for(f in Factors){
  sub       <- LEAP_WT[LEAP_WT$gene_row %in% f & LEAP_WT$gene_col %in% Factors,]
  sub$Level <- 1
  Level1    <- rbind(Level1, sub)
}

# Remove co-expression interactions and subset to our candidates
Level1 <- Level1[Level1$Lag %!in% 0,]
Level1 <- Level1[Level1$gene_col %in% Candidates,]


## iGraph generation for Candidates only and saving to GML file
Network    <- Level1
Genes      <- unique(c(Network$gene_row, Network$gene_col))


nodes <- data.frame(Name   = Genes,
                    TF     = ifelse(Genes %in% TF_list, 1, 0))
edges <- data.frame(From   = Network$gene_row,
                    To     = Network$gene_col,
                    Cor    = Network$Correlation,
                    AbsCor = abs(Network$Correlation),
                    DirCor = ifelse(Network$Correlation > 0, 'Up', 'Down'),
                    Lag    = Network$Lag,
                    Neglag = 1.5-as.numeric(Network$Lag/max(Network$Lag)),
                    Level  = Network$Level) 

graph <- graph_from_data_frame(edges, directed = T, vertices = nodes)
plot(graph, layout = layout.fruchterman.reingold, vertex.size = 8,
     edge.arrow.size = 0.5, vertex.label = NA)
#write_graph(graph, paste0(PATH_SAVE, '/Network_7factors_updated.gml'), format='gml')

#write.table(Network, paste0(PATH_SAVE, '/WT/LEAP/TABLE_7Factors_Network.txt'),quote = F, row.names = F, col.names = T)



#===============================================================================
## RECIPROCAL NFIL3, ID2, TOX, GATA3 TCF7, ZBTB16 & RUNX3 ----------------------
#===============================================================================

Candidates <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_All.txt')[,1]
Factors <- c('Nfil3', 'Id2', 'Tox', 'Gata3', 'Tcf7', 'Zbtb16', 'Runx3')

Level1 <- LEAP_WT[LEAP_WT$gene_row %in% Factors & LEAP_WT$gene_col %in% Factors,]


# Remove co-expression interactions and subset to our candidates
Level1 <- Level1[Level1$Lag %!in% 0,]
Level1 <- Level1[Level1$Correlation > 0.25,]

## iGraph generation for Candidates only and saving to GML file
Network    <- Level1
Genes      <- unique(c(Network$gene_row, Network$gene_col))


nodes <- data.frame(Name   = Genes,
                    TF     = ifelse(Genes %in% TF_list, 1, 0))
edges <- data.frame(From   = Network$gene_row,
                    To     = Network$gene_col,
                    Cor    = Network$Correlation,
                    AbsCor = abs(Network$Correlation),
                    DirCor = ifelse(Network$Correlation > 0, 'Up', 'Down'),
                    Lag    = Network$Lag,
                    Neglag = 1.5-as.numeric(Network$Lag/max(Network$Lag))) 

graph <- graph_from_data_frame(edges, directed = T, vertices = nodes)
plot(graph, layout = layout.fruchterman.reingold, vertex.size = 8,
     edge.arrow.size = 0.5, vertex.label = NA)
#write_graph(graph, paste0(PATH_SAVE, '/Network_Reciprocal.gml'), format='gml')

#write.table(Network, paste0(PATH_SAVE, '/WT/LEAP/TABLE_Reciprocal_Network.txt'),quote = F, row.names = F, col.names = T)

# Save Supplementary Table 2
write.table(Level1, paste0(PATH_SAVE, '/WT/LEAP/Supplementary_Table_2.txt'), quote = F, row.names = F, col.names = T, sep = '\t')











