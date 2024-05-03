#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================




#===============================================================================
## SET UP ----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

PATH      <- 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC'
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/WT')
setwd(PATH)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

source("C:/Users/E15639P/Doctorat/scRNA-seq/Usefull_functions.R")

data    <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))

#
TF.table <- read.table('C:/Users/E15639P/Data/Lists/TFClass_mouse_encoded.txt')

#===============================================================================
## LEAP MATRIX PREPARATION -----------------------------------------------------
#===============================================================================

## Generate matrix and index for ILC trajectory
result_S1 <- LeapMatrix(data, trajectory = "Slingshot_1", suffix = "WT",  min.cells = 50, write = F, dir = paste0(PATH_SAVE, '/LEAP'), slot = "data")
index_S1  <- as.data.frame(result_S1[2])
matrix_S1 <- as.data.frame(result_S1[1])



#===============================================================================
## GENE EXPRESSION VISUALIZATION -----------------------------------------------
#===============================================================================

table <- OrderMatrix(data, "Slingshot_1", min.cells = 50, slot = 'scaledata')

# Optimization
CalculateBin(table, 11)


# Draw superposed curves - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DrawExpr(table, bin.number = 11, by.order = F, feature = c("Nfil3", "Tox", "Zbtb16", "Tcf7", "Id2", "Gata3"), scale = F, ylim = c(-0.8,1.75), col = ColorBlind[c(1,2,3,8,9,10)], lwd = 3,
         superposed = T, compare.with = F, write = F, dir = paste0(PATH, '/Figures'))

# Draw multiple but distinct curves - - - - - - - - - - - - - - - - - - - - - -
output_dir   <- "C:/Users/E15639P/Desktop/plot_florent"
gene_list    <- "Nfil3"
absent_genes <- c()

for(g in gene_list){
  if(g %in% rownames(table)){
    #png(paste(output_dir, "/", g, ".png", sep = ""))
    DrawExpr(table, bin.number = 9, scale = T, by.order = F, feature = g, std = F, col = "darkblue")
    #dev.off()
  }else{
    absent_genes <- c(absent_genes, g)
  }
}

# Draw multiple but distinct curves compared to one gene - - - - - - - - - - - -
gene_list    <- c("Tox", "Tcf7", "Id2", "Zbtb16")
absent_genes <- c()

for(g in gene_list){
  if(g %in% rownames(table)){
    #png(paste(output_dir, "/", g, ".png", sep = ""))
    DrawExpr(table, bin.number = 9, scale = T, by.order = F, feature = g, compare.with = 'Nfil3', std = F, col = ColorBlind[c(8,1)])
    #dev.off()
  }else{
    absent_genes <- c(absent_genes, g)
  }
}

################################################################################
TF_fam <- list()

TF_fam[["BZIP"]] <- c("Cebpa","Cebpb","Cebpg","Cebpd","Cebpe","Ddit3")
TF_fam[["PAR"]]  <- c("Dbp","Hlf","Tef","Tefa","Tefb","Nfil3")
TF_fam[["BZIP_large"]] <- c(
  # Jun factors
  "Cjun","Jun","Junb","Jund",
  # NF-E2-like factors
  "Nfe2","Nfe2l1","Nfe2l2","Nfe2l3","Bach1","Bach2", 
  # ATF-2-like factors
  "Atf2","Atf7","Creb5",
  # Fos factors
  "Cfos","Fos","Fosb","Fosl1", "Fosl2","Fra1","Fra2", 
  # ATF-3-like factors
  "Atf3","Jdp2", 
  # Large Maf factors
  "Cmaf", "Maf","Mafa","Mafb","Nrl",
  # Small Maf factors
  "Maff","Mafg","Mafk","Batf","Batf2","Batf3","Xbp1","Xbp1u","Xbp1s","Atf4",
  "Atf5",
  # CREB-like factors
  "Creb","Atf1","Crem","Cremt","Crema","Cremb","Cremg","Cremtag","Cremt2a",
  "Cremt1a","Cremt1ag","Cremt1g","Cremag",
  # CREB-3-like factors
  "Creb3","Lzip","Creb3l1","Creb3l2","Creb3l3","Creb3l4",
  # ATF-6
  "Atf6a","Atf6b",
  # CREBZF-like factors
  "Crebzf",
  # CREBL2-like factors
  "Crebl2",
  # C/EBP
  "Cebpa","Cebpb","Cebpg","Cebpd","Cebpe","Ddit3",
  # PAR factors
  "Dbp","Hlf","Tef","Tefa","Tefb","Nfil3",
  # TSC22
  "Tsc22d1","Tsc22d2","Tsc22d3","Tsc22d4",
  # UTF
  "Utf1","Spib"
)

# Draw all PAR factors
DrawExpr(table, bin.number = 7, by.order = F, feature = TF_fam[["BZIP_large"]], scale = F, std = F, ylim = c(-1,1), col = c("#26C4EC","#357AB7", "#22427C"), lwd = 3,
         superposed = F, compare.with = F, write = T, dir = paste0(PATH, '/Figures/scRNA-seq/WT/BZIPs'))




#===============================================================================
## K-MEANS CLUSTERING ----------------------------------------------------------
#===============================================================================

old_var <- data@assays[["RNA"]]@var.features
data    <- FindVariableFeatures(data, selection.method = 'mean.var.plot')
new_var <- data@assays[["RNA"]]@var.features

length(new_var[new_var %in% old_var])


## SETUP TABLES ----------------------------------------------------------------

# Read TF file and produce OrderMatrix from scaled slot
TF_list      <- read.table('C:/Users/E15639P/Data/Lists/masterTFlist.txt')[,1]
scaled       <- OrderMatrix(data, "Slingshot_1", slot = 'scaledata')
pseudotime   <- scaled['Pseudotime',]
BIN          <- 12

# Filtering TF list to be in rownames and variables
TF           <- rownames(scaled)[rownames(scaled) %in% TF_list]
TF           <- TF[TF %in% data@assays[["RNA"]]@var.features]
#TF           <- new_var # get all varible genes


# Create tables
scaled.TF    <- scaled[rownames(scaled) %in% c('Pseudotime', TF),] # TF + pseudotime
scaled_kmean <- scaled.TF[rownames(scaled.TF) %!in% 'Pseudotime',] # only TF
matrix_expr  <- DrawExpr(scaled.TF, feature.list = TF, bin.number = BIN, scale = F,
                        scale.method = 2, std = F, superposed = T, by.order = F)

## NUMBER OF CLUSTERS ----------------------------------------------------------

set.seed(777)

# Determining optimal number of clusters
Clustering_list        <- list()
Tot.Wthnss             <- c()
for(i in 1:30){
  cluster_num_test     <- kmeans(matrix_expr, i)
  Clustering_list[[i]] <- cluster_num_test
  Tot.Wthnss           <- c(Tot.Wthnss, cluster_num_test[["tot.withinss"]])
}

ggplot(as.data.frame(Tot.Wthnss), aes(x= c(1:30), y = Tot.Wthnss)) +
  geom_line() + geom_point() + 
  #geom_vline(xintercept = 6, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(breaks = c(1:30)) +
  theme_classic() +
  xlab('Number of clusters k') + ylab('Total Within Sum of Square')

# Cluster number selection = 6
cluster_num          <- 8
scaled_kmean$Cluster <- Clustering_list[[cluster_num]][["cluster"]]
centroid_matrix      <- as.data.frame(Clustering_list[[cluster_num]][["centers"]])


## CLUSTERS VISUALIZATION ------------------------------------------------------

test = 'Hmga2'

# Plot all clusters
for(i in 1:cluster_num){
  # Draw all genes from cluster i
  group <- rownames(subset(scaled_kmean, Cluster %in% i))
  tab   <- scaled.TF[rownames(scaled.TF) %in% c('Pseudotime', group),]
  
  # Save gene list
  #write.table(rownames(tab)[rownames(tab) %!in% 'Pseudotime'], paste0(PATH_SAVE, '/Kmeans/K', i, '.txt'), row.names = F, col.names = F, quote = F)
  
  DrawExpr(tab, feature.list = TF, bin.number = BIN, scale = F, scale.method = 0,
           std = F, superposed = T, by.order = F, main = '', ylim = c(-1,1),
           col = c(rep('lightgrey', nrow(tab))))
  par(new=TRUE)
  # Add centroid curve
  plot(0:(BIN-1), centroid_matrix[i,], type = "l", lwd = 2, ylim = c(-1, 1),
       col = ColorBlind[14], ylab = '', xlab = '', main = paste('Cluster', i))
  par(new=TRUE)
  # Add famous genes
  
  if(test %in% rownames(tab)){
    DrawExpr(tab, feature.list = test, bin.number = BIN, scale = F, scale.method = 2,
             std = F, superposed = T, by.order = F, main = '', ylim = c(-1,1),
             col = ColorBlind[rep(8,100)])
  }
  par(new=FALSE)
  
}

scaled_kmean$Cluster[rownames(scaled_kmean) %in% 'Nfil3']
rownames(scaled_kmean[scaled_kmean$Cluster %in% 16,])



#===============================================================================
## LEAP RESULTS ANALYSIS -------------------------------------------------------
#===============================================================================

# Annotate Results - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#LEAP_WT  <- read.table(paste0(PATH_SAVE, '/LEAP/MAC_results_WT.txt'), header = T)
#index_WT <- read.table(paste0(PATH_SAVE, '/LEAP/gene_index_WT.txt'), header = T)
#thr      <- 0.1625
#Annotate <- AnnotateLeap(LEAP_WT, index_WT, write = T, dir = paste0(PATH_SAVE, '/LEAP'))
#Annotate <- Annotate[abs(Annotate$Correlation) > 0.1625,]
#write.table(Annotate, paste0(PATH_SAVE, '/LEAP/results_indexed_filtered.txt'), row.names = F, quote = F)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Annotate <- read.table(paste0(PATH_SAVE, '/LEAP/results_indexed_filtered.txt'), header = T)
TF.list  <- read.table('C:/Users/E15639P/Data/Lists/masterTFlist.txt')[,1]
table    <- OrderMatrix(data, "Slingshot_1", min.cells = 50, slot = 'scaledata')

NF       <- Annotate[Annotate$gene_row %in% 'Nfil3' & Annotate$Lag > 0,]
#write.table(NF, paste0(PATH_SAVE, '/LEAP/TABLE_Nfil3_Candidates.txt'), quote = F, row.names = F, col.names = T)
#write.table(NF$gene_col, paste0(PATH_SAVE, '/LEAP/Nfil3_Candidates.txt'), quote = F, row.names = F, col.names = F)


## CROSSING BULK RESULTS -------------------------------------------------------
Bulk_UP   <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Bulk-RNA/Candidates_UP.txt')[,1]
Bulk_DOWN <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Bulk-RNA/Candidates_DOWN.txt')[,1]

NF_UP   <- NF$gene_col[NF$Correlation > 0]
NF_DOWN <- NF$gene_col[NF$Correlation < 0]

Cross.UP   <- NF_UP[NF_UP %in% Bulk_UP]
Cross.DOWN <- NF_DOWN[NF_DOWN %in% Bulk_DOWN]


#write.table(NF_UP, paste0(PATH_SAVE, '/LEAP/Candidates_LEAP_UP.txt'), row.names = F, col.names = F, quote = F)
#write.table(NF_DOWN, paste0(PATH_SAVE, '/LEAP/Candidates_LEAP_DOWN.txt'), row.names = F, col.names = F, quote = F)

#write.table(Cross.UP, 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_UP.txt', row.names = F, col.names = F, quote = F)
#write.table(Cross.DOWN, 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_DOWN.txt', row.names = F, col.names = F, quote = F)
#write.table(c(Cross.UP, Cross.DOWN), 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_All.txt', row.names = F, col.names = F, quote = F)


for(g in c(NFIL3.am$gene_row, NFIL3.av$gene_col)){
  DrawExpr(table, feature.list = g, bin.number = 9, 
           by.order = F, compare.with = 'Nfil3',
           scale = T, scale.method = 0, 
           col = ColorBlind[c(2,8)],
           write = T, dir = paste0(output_dir, '/WT'))
}





################################################################################




