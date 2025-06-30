#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
# Pseudotemporal trajectories reconstruction from WT ILC model
#
# Load input files generated during previous steps
# Perform pseudotime reconstruction using Slingshot
# Save pseudotime cell order by trajectories in Seurat Object metadata
# Analyse samples distribution inside trajectories
# Draw gene expression variations across pseudotime
# Perform K-means clustering on all TF expression variations
# Save RDS file used as input for the next step



#===============================================================================
## SET UP ----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PROJECT INFO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

PATH <- 'C:/Users/E15639P/Desktop/NFIL3_dev_ILC'

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

setwd(PATH)

# Load Packages and custom functions
source('C:/Users/E15639P/Desktop/NFIL3_dev_ILC/Custom_Functions.R')
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(tidyverse))

# Create required directories
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq/WT'), 'Kmeans'))
dir.create(file.path(paste0(PATH, '/Saves/scRNA-seq/WT'), 'Kmeans'))

# Set up path for figures and saves
PATH_FIG  <- paste0(PATH, '/Figures/scRNA-seq/WT')
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/WT')

CellTypeOrder <- c('ALP', 'overLIP', 'TULIP', 'ILCpro')



#===============================================================================
## INPUT -----------------------------------------------------------------------
#===============================================================================

# Read files
data    <- readRDS(paste0(PATH_SAVE, '/4_Filtered.rds'))
TF_list <- read.table('C:/Users/E15639P/Data/Lists/masterTFlist.txt')[,1]
TF_fam  <- readRDS(paste0(PATH, '/Saves/DNAse-seq/TF_fam.RDS'))


# Set reduction to use and plot axis limits
red    <- 'rev_umap_cssz2'
xlimit <- c(min(data@reductions[[red]]@cell.embeddings[,1]),
            max(data@reductions[[red]]@cell.embeddings[,1]))
ylimit <- c(min(data@reductions[[red]]@cell.embeddings[,2]),
            max(data@reductions[[red]]@cell.embeddings[,2]))      



#===============================================================================
## SLINGSHOT -------------------------------------------------------------------
#===============================================================================
#Based on : https://bustools.github.io/BUS_notebooks_R/slingshot.html

DimPlot(data, reduction = red, group.by = 'seurat_clusters', 
        label = TRUE, repel = TRUE, 
        cols = c(ColorBlind[c(3,9,2,7,8,10)], "purple"), pt.size = 0.9) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')


## PSEUDOTIME RECONSTRUCTION ---------------------------------------------------
sce                                     <- as.SingleCellExperiment(data)     
colData(sce)$Seurat_clusters            <- as.character(data@active.ident)

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', 
                 reducedDim = toupper(red), start.clus = 2, end.clus = c(0,5))


## ADDING PSEUDOTIMES TO DATA --------------------------------------------------
pseudotime                              <- slingPseudotime(sce)
for(i in 1:ncol(pseudotime)){
  data@meta.data[[paste0('Slingshot_', i)]] <- pseudotime[,i]
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Saving SCE file after processing
saveRDS(sce, paste0(PATH_SAVE, '/Slingshot.rds'))
sce <- readRDS(paste0(PATH_SAVE, '/Slingshot.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## PLOTTING PSEUDOTIME TRAJECTORIES --------------------------------------------

# Specifying gradient colors for plots
gradient         <- hcl.colors(100)

# Plot with curves
curves <- slingCurves(sce, as.df = TRUE) %>%
  dplyr::rename("x" = paste0(str_replace_all(
    toupper(red), "[[:punct:]]", ""), '_1'), 
    "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

# First Trajectory
plot12 <- FeaturePlot(data, 'Slingshot_1', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            linewidth = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot12, PATH_FIG)

# Second Trajectory
plot13 <- FeaturePlot(data, 'Slingshot_2', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            linewidth = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot13, PATH_FIG)


## CELLTYPE DISTRIBUTION -------------------------------------------------------
# Calculate intermediates distribution inside ILC lineage

data.QC <- readRDS(paste0(PATH_SAVE, '/1_AfterQC.rds'))

compare <- data.frame(Slingshot = data$Slingshot_1, Sample = data$orig.ident)
compare <- compare[order(compare$Slingshot),]
compare <- compare[!is.na(compare$Slingshot),]
compare$Slingshot_order <- c(1:nrow(compare))


## TULIP proportion in ILC trajectory
TULIP_total   <- length(
  data.QC@meta.data$orig.ident[data.QC@meta.data$orig.ident %in% 'TULIP'])
TULIP_ILC     <- nrow(compare[compare$Sample %in% 'TULIP',])
fraction      <- TULIP_ILC/TULIP_total*100
fraction # 61/320 -> 19.0625%

overLIP_total <- length(
  data.QC@meta.data$orig.ident[data.QC@meta.data$orig.ident %in% 'overLIP'])
overLIP_ILC   <- nrow(compare[compare$Sample %in% 'overLIP',])
fraction      <- overLIP_ILC/overLIP_total*100 
fraction # 918/1073 -> 85.55452%


# Draw distributed density
plot14 <- ggplot(
  compare, aes(x=Slingshot, y=..scaled.., 
               fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(alpha=0.5, linewidth = 0.6) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,3,1,9)], name = '')
writePlot(plot14, PATH_FIG)

# Draw cumulated density
plot15 <- ggplot(
  data=compare, aes(x=Slingshot, 
                    fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(adjust=1.5, position="fill", linewidth = 0.5) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,3,1,9)], name = '')
writePlot(plot15, PATH_FIG)

# Draw linear density
plot16 <- ggplot(
  compare, aes(x=Slingshot, y = 1, 
               colour=factor(Sample, levels = CellTypeOrder))) +
  geom_point(shape=124) +
  theme_classic()+
  scale_color_manual(values = ColorBlind[c(10,3,1,9)], name = '')
writePlot(plot16, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Saving CDS file after processing and object data after adding metadata
saveRDS(data, paste0(PATH_SAVE, '/5_Pseudotime.rds'))
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## GENE EXPRESSION VISUALIZATION -----------------------------------------------
#===============================================================================

table <- OrderMatrix(data, 'Slingshot_1', slot = 'scaledata')

# Draw superposed curves - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DrawExpr(table, bin.number = 12, by.order = F, 
         feature = c('Nfil3','Tox','Zbtb16','Tcf7','Id2','Gata3'), scale = F, 
         ylim = c(-1,1.75), col = ColorBlind[c(1,2,3,8,9,10)], 
         lwd = 2, superposed = T)
plot17 <- recordPlot()
writePlot(plot17, PATH_FIG)

# Draw all BZIP factors
DrawExpr(table, bin.number = 12, by.order = F, feature = TF_fam[['BZIP_large']], 
         scale = F, std = F, ylim = c(-1,1), col = c("#26C4EC"), lwd = 2,
         superposed = F)



#===============================================================================
## K-MEANS CLUSTERING ----------------------------------------------------------
#===============================================================================

# Establish variable features after filtering
data    <- FindVariableFeatures(data, selection.method = 'mean.var.plot')
new_var <- data@assays[["RNA"]]@var.features


## SETUP TABLES ----------------------------------------------------------------
scaled       <- table
pseudotime   <- scaled[,'Pseudotime']
BIN          <- 12

# Filtering TF list to be in rownames and variables
TF           <- colnames(scaled)[colnames(scaled) %in% TF_list]
TF           <- TF[TF %in% data@assays[["RNA"]]@var.features]


# Create tables
scaled_TF    <- scaled[,colnames(scaled) %in% c('Pseudotime', TF)]
scaled_kmean <- scaled_TF[,colnames(scaled_TF) %!in% 'Pseudotime'] 
matrix_expr  <- DrawExpr(scaled_TF, feature.list = TF, bin.number = BIN, 
                         scale = F, std = F, superposed = T, by.order = F,
                         color = rep('lightgrey', 75))

## NUMBER OF CLUSTERS ----------------------------------------------------------

set.seed(777)

# Determining optimal number of clusters
Clustering_list        <- list()
Tot.Wthnss             <- c()
scaled_kmean_T         <- as.data.frame(t(scaled_kmean))
for(i in 1:30){
  cluster_num_test     <- kmeans(scaled_kmean_T, i, nstart = 25)
  Clustering_list[[i]] <- cluster_num_test
  Tot.Wthnss           <- c(Tot.Wthnss, cluster_num_test[["tot.withinss"]])
}

plot18 <- ggplot(as.data.frame(Tot.Wthnss), aes(x= c(1:30), y = Tot.Wthnss)) +
  geom_line() + geom_point() + 
  scale_x_continuous(breaks = c(1:30)) +
  theme_classic() +
  xlab('Number of clusters k') + ylab('Total Within Sum of Square')
writePlot(plot18, PATH_FIG)

# Cluster number selection = 8
cluster_num     <- 8
kmean_groups    <- data.frame(Clstr=Clustering_list[[cluster_num]][["cluster"]])
centroid_matrix <- data.frame(t(Clustering_list[[cluster_num]][["centers"]]))


## CLUSTERS VISUALIZATION ------------------------------------------------------

# Plot all clusters
for(i in 1:cluster_num){
  # Draw all genes from cluster i + centroid
  group        <- rownames(kmean_groups)[kmean_groups$Clstr %in% i]
  tab          <- scaled_TF[,colnames(scaled_TF) %in% c('Pseudotime', group)]
  tab$centroid <- centroid_matrix[,i]
  
  # Save gene list
  write.table(group, 
              paste0(PATH_SAVE, '/Kmeans/K', i, '.txt'), row.names = F, 
              col.names = F, quote = F)
  
  DrawExpr(tab, feature.list = c(group, 'centroid'), bin.number = BIN, 
           scale = F, superposed = T, by.order = F, main = paste0('Cluster',i), 
           ylim = c(-1,2), col = c(rep('lightgrey', length(group)), 'black'))
  
  plot <- recordPlot()
  writePlot(plot, paste0(PATH_FIG, '/Kmeans/'), paste0('K', i, '.png'))
}

# Look for the cluster number of a specific gene
kmean_groups['Nfil3',]

# Display all genes belonging to a specific cluster
rownames(kmean_groups)[kmean_groups$Clstr %in% 6]


