#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
# Based on : https://satijalab.org/seurat/archive/v3.1/integration.html
# https://github.com/quadbiolab/simspec/blob/master/vignette/vignette.md
#
# Integration CSS of WT ILC Il7LT+ intermediates
#
# Load datasets from distinct batches and convert them in a Seurat Object
# Perform QC and filter samples separately
# Lognormalize and scale merged datasets
# Apply cell cycle and mithochondrial genes regression during scaling
# Integrate experiments using Cluster Similarity Spectrum (CSS)
# Perform clustering and identify cell types using ImmGen.org
# Quantify and remove outliers based on this identification
# Visualize key genes expression
# Save RDS file used as input for the next step



#===============================================================================
## SETUP -----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PROJECT INFO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

PATH         <- 'C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC'
DATA_DIR     <- 'C:/Users/E15639P/Data/scRNA-seq/SingleCell_mm39'

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

setwd(PATH)

# Load Packages and custom functions
source('C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Scripts/Custom_Functions.R')
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(simspec))
suppressPackageStartupMessages(library(Nebulosa))

# Create required directories
dir.create(file.path(PATH, '/Figures'))
dir.create(file.path(PATH, '/Saves'))

dir.create(file.path(paste0(PATH, '/Figures'), 'scRNA-seq'))
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq'), 'WT'))
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq/WT'), 'Markers'))
dir.create(file.path(paste0(PATH, '/Saves'), 'scRNA-seq'))
dir.create(file.path(paste0(PATH, '/Saves/scRNA-seq'), 'WT'))

# Set up path for figures and saves
PATH_FIG  <- paste0(PATH, '/Figures/scRNA-seq/WT')
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/WT')

CellTypeOrder <- c('ALP', 'overLIP', 'TULIP', 'ILCpro')



#===============================================================================
## INPUT -----------------------------------------------------------------------
#===============================================================================

# Opening the datasets and turning them into Seurat object
# Old Experiments - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# All Lymphoid Progenitors
ALP       <- Read10X(
  data.dir = paste0(DATA_DIR, '/1_ABSC_CLP')) %>% 
  CreateSeuratObject(project = 'ALP', min.cells = 1, min.features = 0)

# sEILP, cEILP and ILC progenitors 
ILCpro    <- Read10X(
  data.dir = paste0(DATA_DIR, '/3_ABSC_EILP_WT')) %>%
  CreateSeuratObject(project = 'ILCpro', min.cells = 1, min.features = 0)

# a4b7 intermediates
TULIP     <- Read10X(
  data.dir = paste0(DATA_DIR, '/2_ABSC_TULIP')) %>% 
  CreateSeuratObject(project = 'TULIP', min.cells = 1, min.features = 0)

# Recent Experiment - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Il7rLT+ XT
overLIP   <- Read10X(
  data.dir = paste0(DATA_DIR, '/2019_scRNAseq_data/WT')) %>%
  CreateSeuratObject(project = 'overLIP', min.cells = 1, min.features = 0)


# Spread datasets in groups of experiment
ALP@meta.data[['Experiment']]     <- 'Experiment 1'
ILCpro@meta.data[['Experiment']]  <- 'Experiment 1'
TULIP@meta.data[['Experiment']]   <- 'Experiment 1'
overLIP@meta.data[['Experiment']] <- 'Experiment 2'

# Merge all the experiments to perform the quality check 
data  <- merge(TULIP, c(ILCpro, ALP, overLIP), 
               add.cell.ids = c('TULIP', 'ILCpro', 'ALP', 'overLIP'), 
               project = 'scRNA_Integration')
data@meta.data[['Project']] <- 'scRNA_WT_Integration'

# Organize cell type in the development order for the plots
data@active.ident <- factor(data@active.ident, levels = CellTypeOrder)



#===============================================================================
## QUALITY CHECK ---------------------------------------------------------------
#===============================================================================

# Calculate the percent of mitochondrial genes expression
data$Percent.mt <- PercentageFeatureSet(data, pattern = '^mt-')

plot0   <- VlnPlot(data, features = c('nFeature_RNA','Percent.mt'), 
                   group.by = 'Experiment', cols = ColorBlind[c(1, 3)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot0, PATH_FIG)

# Elimination of cells that over-expressed mitochondrial genes (dead cells)
data.QC <- subset(data, subset = ((Experiment == 'Experiment 1' & 
                                     nFeature_RNA > 200 & 
                                     nFeature_RNA < 5000 & Percent.mt < 2) |
                                  (Experiment != 'Experiment 1' & 
                                     nFeature_RNA > 1200 & 
                                     nFeature_RNA < 5000 & Percent.mt < 4))) 

plot1   <- VlnPlot(data.QC, features = c('nFeature_RNA', 'Percent.mt'), 
                   group.by = 'Experiment', cols = ColorBlind[c(1, 3)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot1, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after quality check
#saveRDS(data.QC, paste0(PATH_SAVE, '/1_AfterQC.rds'))
data.QC <- readRDS(paste0(PATH, '/1_AfterQC.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## NORMALIZATION AND SCALING ---------------------------------------------------
#===============================================================================

data.combined <- NormalizeData(data.QC, normalization.method = 'LogNormalize')

# Cell Cycle Genes identification and scoring
s.genes       <- capitalize(tolower(cc.genes$s.genes))
g2m.genes     <- capitalize(tolower(cc.genes$g2m.genes))
data.combined <- CellCycleScoring(data.combined, s.features = s.genes, 
                                  g2m.features = g2m.genes, set.ident = TRUE)


# Calculate difference between S and G2M is described as more relevant
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Alternate Workflow
data.combined$CC.Difference <- data.combined$S.Score - data.combined$G2M.Score
data.combined <- ScaleData(data.combined, 
                           vars.to.regress = c('CC.Difference','Percent.mt'), 
                           verbose = TRUE, features = rownames(data.combined), 
                           do.scale = TRUE)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after scaling step
#saveRDS(data.combined, paste0(PATH_SAVE, '/2_Scaled.rds'))
data.combined <- readRDS(paste0(PATH_SAVE, '/2_Scaled.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 


#===============================================================================
## CLUSTERING SIMILARITY SPECTRUM INTEGRATION ----------------------------------
#===============================================================================

# Before Integration - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

data.combined <- FindVariableFeatures(data.combined, nfeatures = 3000)
data.combined <- RunPCA(data.combined, npcs = 50)

plot2         <- ElbowPlot(data.combined, 30) +
  labs(title = '', y = 'Standard deviation', x = 'Component number') + 
  theme(plot.title = element_text(hjust = 0.5))
writePlot(plot2, PATH_FIG)

data.combined <- RunUMAP(data.combined, dims = 1:16, metric = 'euclidean', 
                         n.neighbors = 50)

plot3         <- DimPlot(data.combined, group.by = 'Experiment', 
                         reduction = 'umap', cols = ColorBlind[c(1,10)], 
                         pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot3, PATH_FIG)

plot4         <- DimPlot(data.combined, group.by = 'orig.ident', 
                         reduction = 'umap', cols = ColorBlind[c(10,3,1,9)], 
                         pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot4, PATH_FIG)


# After Integration Z-Transform  - - - - - - - - - - - - - - - - - - - - - - - -  

data          <- cluster_sim_spectrum(data.combined, label_tag='Experiment', 
                                      cluster_resolution = 0.6, 
                                      corr_method = 'spearman', lambda = 50, 
                                      reduction.name = 'cssz', 
                                      reduction.key = 'CSSZ_')
data          <- RunUMAP(data, reduction = 'cssz', 
                         dims = 1:ncol(Embeddings(data, 'cssz')), 
                         n.neighbors = 50, reduction.name='umap_cssz', 
                         reduction.key='UMAPCSSZ_')

plot5         <- DimPlot(data, group.by = 'Experiment', reduction = 'umap_cssz', 
                         cols = ColorBlind[c(1,10)], pt.size = 0.5, 
                         order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') + scale_x_reverse()
writePlot(plot5, PATH_FIG)

plot6         <- DimPlot(data, group.by = 'orig.ident', reduction = 'umap_cssz', 
                         cols = ColorBlind[c(10,3,1,9)], pt.size = 0.5, 
                         order = rev(CellTypeOrder)) +
  labs(title = "", x = 'UMAP1', y = 'UMAP2') + scale_x_reverse()
writePlot(plot6, PATH_FIG)

# Check cell cycle distribution
FeaturePlot(data, 'S.Score', reduction = 'umap_cssz') + scale_x_reverse()
FeaturePlot(data, 'G2M.Score', reduction = 'umap_cssz') + scale_x_reverse()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after integration
#saveRDS(data, paste0(PATH_SAVE, '/3_Integrated.rds'))
data <- readRDS(paste0(PATH_SAVE, '/3_Integrated.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## CLUSTERING AND OUTLIER ELIMINATION ------------------------------------------
#===============================================================================

data          <- FindNeighbors(data, reduction = 'umap_cssz', 
                               dims = 1:ncol(Embeddings(data, 'umap_cssz')))
data          <- FindClusters(data)

plot7         <- DimPlot(data, reduction = 'umap_cssz', 
                         group.by = 'seurat_clusters', label = TRUE, 
                         repel = TRUE, cols = c(rep(ColorBlind,10)), 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2') + scale_x_reverse()
writePlot(plot7, PATH_FIG)


# Identification of cell types
clusterMarker <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
top30         <- clusterMarker %>% group_by(cluster) %>% top_n(30, avg_log2FC) 
#write.csv(top30, paste0(PATH_SAVE, '/WT_Top30markers.csv'), row.names = FALSE)

# - - - - - - - - - - - - - - - - #
# Based on ImmGen.org annotation  #
# - - - - - - - - - - - - - - - - #
# 18 = ILC2 (Icos, Bcl11b)        #
# 19 = B cells                    #
# 20 = Mast Basophils             #
# 21 = ILC2 (Icos, Bcl11b, Rora)  #
# 22 = Macro/Mono/Granulo         #
# 23 = Macro/Mono/Granulo         #
# - - - - - - - - - - - - - - - - #


## KEY GENES REPRESENTATION ----------------------------------------------------

markers_genes <- c(
  'Nfil3',
  # B Cells
  'Pax5', 'Ebf1', 
  # Mast Cells
  'Cpa3', 'Gzmb', 'Kit', 
  # Mono/Macrophages
  'Lgals3', 'Lyz2', 
  # ILC2
  'Il12ra', 'Icos', 'Stab2', 'Bcl11b', 'Tox2'
)

for(g in markers_genes){
  if(g %in% rownames(data)){
    density <- plot_density(data, g, pal = 'magma', 
                            size = 2) + scale_x_reverse()
    writePlot(density, paste0(PATH_FIG, '/Markers'), 
              paste0('integrated_', g, '_DensityPlot.png'))
    
    feature <- FeaturePlot(data, g, order = T, pt.size = 1, 
                           reduction = 'umap_cssz') + scale_x_reverse()
    writePlot(feature, paste0(PATH_FIG, '/Markers'), 
              paste0('integrated_', g, '_FeaturePlot.png'))
  }
}


## CONTAMINANT PERCENTAGE ------------------------------------------------------

# Overall outliers
Outliers  <- subset(data, seurat_clusters %in% c(18,19,20,22,23))
Inliers   <- subset(data, seurat_clusters %!in% c(18,19,20,21,22,23))

prop_out  <- table(factor(Outliers$orig.ident, levels = CellTypeOrder))
prop_in   <- table(factor(Inliers$orig.ident, levels = CellTypeOrder))
prop_all  <- (prop_out/(prop_in+prop_out))*100

plot8     <- barplot(prop_all[c(2,3)], col = ColorBlind[c(3,1)], 
                     border = 'white', ylab = 'Outlier percentage', 
                     main = 'Percentage of outliers' , ylim = c(0, 60))
y         <- prop_all[c(2,3)]
text(plot8, y+2, labels = paste(format(y, scientific=FALSE, digits = 2 ), '%'))
plot8     <- recordPlot()
writePlot(plot8, PATH_FIG)


# Outliers by cellype
prop_table           <- as.data.frame(matrix(0,4,4))
colnames(prop_table) <- c('ILC2', 'B_cells', 'Mast_Baso', 'Macro_Mono_Granulo')
rownames(prop_table) <- names(table(Outliers$orig.ident))

for(i in 1:length(Outliers$seurat_clusters)){
  if(Outliers$seurat_clusters[i] %in% c(18,21)){
    prop_table[Outliers$orig.ident[i],'ILC2'] <- prop_table[
      Outliers$orig.ident[i],'ILC2']+1
  }else if(Outliers$seurat_clusters[i] %in% 19){
    prop_table[Outliers$orig.ident[i],'B_cells'] <- prop_table[
      Outliers$orig.ident[i],'B_cells']+1
  }else if(Outliers$seurat_clusters[i] %in% 20){
    prop_table[Outliers$orig.ident[i],'Mast_Baso'] <- prop_table[
      Outliers$orig.ident[i],'Mast_Baso']+1
  }else if(Outliers$seurat_clusters[i] %in% c(22,23)){
    prop_table[Outliers$orig.ident[i],'Macro_Mono_Granulo'] <- prop_table[
      Outliers$orig.ident[i],'Macro_Mono_Granulo']+1
  }
}

prop_table <- rbind(prop_table, Total = colSums(prop_table))
perc_table <- prop_table[1:4,]/prop_table[rep(5,4),]*100
bar_table  <- data.frame(Percentage = c(perc_table[,1],perc_table[,2], 
                                        perc_table[,3], perc_table[,4]),
                         CellType = rep(colnames(perc_table), each = 4),
                         Sample = rep(rownames(perc_table)))

plot9 <- ggplot(bar_table[bar_table$Sample %in% c('TULIP', 'overLIP'),], 
                aes(fill=Sample, y=Percentage, x=CellType)) + 
  geom_bar(position='stack', stat='identity') + 
  ggtitle('Proportions of celltypes in outleirs') +
  scale_fill_manual(values= ColorBlind[c(3,1)]) + 
  theme_classic()
writePlot(plot9, PATH_FIG)


## DATA FILTERING --------------------------------------------------------------

# Remove mature cells (18), added B cells (19) and contaminants (20,22,23)
data  <- Inliers

data  <- RunUMAP(data, reduction = 'cssz', 
                 dims = 1:ncol(Embeddings(data, 'cssz')), n.neighbors = 50, 
                 reduction.name='umap_cssz2', reduction.key='UMAPCSSZ2_')
data  <- FindNeighbors(data, reduction = 'umap_cssz2', 
                      dims = 1:ncol(Embeddings(data, 'umap_cssz2')))
data  <- FindClusters(data, resolution = 0.1)


# Create inverted embedding to simplify representation
embed <- as.matrix(data.frame(
  -data@reductions[['umap_cssz2']]@cell.embeddings[,1],
  data@reductions[['umap_cssz2']]@cell.embeddings[,2]))

colnames(embed)           <- c('REVUMAPCSSZ2_1', 'REVUMAPCSSZ2_2')
data[['rev_umap_cssz2']]  <- CreateDimReducObject(embeddings = embed, 
                                                  key = 'REVUMAPCSSZ2_', 
                                                  assay = DefaultAssay(data))
data@reductions[['rev_umap_cssz2']]@global <- T


plot10 <- DimPlot(data, group.by = 'orig.ident', reduction = 'rev_umap_cssz2', 
                 cols = ColorBlind[c(10,3,1,9)], pt.size = 1, 
                 order = rev(CellTypeOrder)) + 
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot10, PATH_FIG)

plot11 <- DimPlot(data, reduction = 'rev_umap_cssz2', 
                  group.by = 'seurat_clusters', 
                  label = TRUE, repel = TRUE, 
                  cols = c(ColorBlind[c(3,9,2,7,8,10)], "purple"), 
                  pt.size = 1) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot11, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint of subsets
#saveRDS(data, paste0(PATH_SAVE, '/4_Filtered.rds'))
data <- readRDS(paste0(PATH_SAVE, '/4_Filtered.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## GENE VISUALIZATION ----------------------------------------------------------
#===============================================================================
# Based on markers used in : doi:10.1038/s41590-019-0445-7

## DC SCORING ------------------------------------------------------------------

# Loading a subset of DC and pre-DC markers
DC_table    <- read.csv('C:/Users/E15639P/Data/Lists/DC_enrichement_list.csv', 
                        header = TRUE, sep = ';', row.names = 1)
DC_genelist <- list(row.names(DC_table[(DC_table$p_val < 0.01 & 
                                          DC_table$avg_logFC > 0.4),]))

# Adding DC_score as metadata
data        <- AddModuleScore(data, DC_genelist, name = 'DC_score')

# Score spatial repartition in cells
plot_density(data, 'DC_score1', pal = 'magma', size = 2)


## GENE VISUALIZATION ----------------------------------------------------------

markers_genes <- c(
  # Genes of interest
  'Bcl11b', 'Cd74', 'Flt3', 'Tcf7', 'Zbtb16', 'Tox', 'Id2', 'Gata3', 'Batf', 
  'Zeb2',
  # Expressed BZIPs
  'Nfil3', 'Cebpa', 'Batf3', 'Creb3l2',
  # Other BZIPs
  'Atf1', 'Bach2', 'Cebpb', 'Jun', 'Maf', 'Mafg', 'Mafk', 'Tsc22d1',
  # K-means cluster
  'Aff3', 'Nfkb1', 'Runx3', 'Sox5', 'Dach1', 'Mycn'
)

for(g in markers_genes){
  density <- plot_density(data, g, pal = 'magma', size = 2)
  writePlot(density, paste0(PATH_FIG, '/Markers'), 
            paste0('filtered_', g, '_DensityPlot.png'))
  
  feature <- FeaturePlot(data, g, order = T, reduction = 'rev_umap_cssz2', 
                         pt.size = 2)
  writePlot(feature, paste0(PATH_FIG, '/Markers'), 
            paste0('filtered_', g, '_FeaturePlot.png'))
}



