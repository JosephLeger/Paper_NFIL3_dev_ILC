#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
# Based on : https://satijalab.org/seurat/archive/v3.1/integration.html
# https://github.com/quadbiolab/simspec/blob/master/vignette/vignette.md
#
# Integration CSS of WT and NFIL3-KO ILC Il7LT+ intermediates
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
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq'), 'KO'))
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq/KO'), 'Markers'))
dir.create(file.path(paste0(PATH, '/Saves/scRNA-seq'), 'KO'))

# Set up path for figures and saves
PATH_FIG  <- paste0(PATH, '/Figures/scRNA-seq/KO')
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/KO')

CellTypeOrder     <- c('ALP', 'TULIP', 'ILCpro', 'overLIP', 'NF-KO', 'TOX-KO')



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

# Recent Experiments - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Il7rLT+ WT 
overLIP   <- Read10X(
  data.dir = paste0(DATA_DIR, '/2019_scRNAseq_data/WT')) %>%
  CreateSeuratObject(project = 'overLIP', min.cells = 1, min.features = 0)

# Il7rLT+ NFIL3-KO
NF.KO     <- Read10X(
  data.dir = paste0(DATA_DIR, '/2019_scRNAseq_data/NFIL3_KO')) %>%
  CreateSeuratObject(project = 'NF-KO', min.cells = 1, min.features = 0)

# Il7rLT+ TOX-KO
TOX.KO    <- Read10X(
  data.dir = paste0(DATA_DIR, '/2019_scRNAseq_data/TOX_KO')) %>%
  CreateSeuratObject(project = 'TOX-KO', min.cells = 1, min.features = 0)


# Spread datasets in groups of experiment
ALP@meta.data[['Experiment']]         <- 'Experiment 1'
ILCpro@meta.data[['Experiment']]      <- 'Experiment 1'
TULIP@meta.data[['Experiment']]       <- 'Experiment 1'
overLIP@meta.data[['Experiment']]     <- 'Experiment 2'
NF.KO@meta.data[['Experiment']]       <- 'Experiment 2'
TOX.KO@meta.data[['Experiment']]      <- 'Experiment 2'

# Merge all the experiments to perform the quality check 
data  <- merge(TULIP, c(ILCpro, ALP, overLIP, NF.KO, TOX.KO), add.cell.ids = c(
  'TULIP', 'ILCpro', 'ALP', 'overLIP', 'NF-KO', 'TOX-KO'), 
               project = 'scRNA_Integration')
data@meta.data[['Project']] <- 'scRNA_ALL_Integration'

# Organize cell type in the development order for the plots
data@active.ident    <- factor(data@active.ident, levels = CellTypeOrder)
data[['orig.ident']] <- factor(data@meta.data$orig.ident, 
                               levels = CellTypeOrder)



#===============================================================================
## QUALITY CHECK ---------------------------------------------------------------
#===============================================================================

# Calculate the percent of mitochondrial genes expression
data$Percent.mt <- PercentageFeatureSet(data, pattern = '^mt-')

plot0   <- VlnPlot(data, features = c('nFeature_RNA','Percent.mt'), 
                   group.by = 'Experiment', cols = ColorBlind[c(1,3)]) + 
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
                   group.by = 'Experiment', cols = ColorBlind[c(1,3)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot1, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after quality check
saveRDS(data.QC, paste0(PATH_SAVE, '/1_AfterQC.rds'))
data.QC <- readRDS(paste0(PATH, '/1_AfterQC.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## NORMALIZATION AND SCALING ---------------------------------------------------
#===============================================================================

data.combined <- NormalizeData(data.QC, normalization.method = 'LogNormalize')

# Cell Cycle Genes identification and scoring
suppressPackageStartupMessages(library(Hmisc))
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

## BEFORE INTEGRATION ---------------------------------------------------------- 

data.combined <- FindVariableFeatures(data.combined, nfeatures = 3000)
data.combined <- RunPCA(data.combined, npcs = 50)

plot2         <- ElbowPlot(data.combined, 30) +
  labs(title = "", y ='Standard deviation', x = 'Component number') + 
  theme(plot.title = element_text(hjust = 0.5))
writePlot(plot2, PATH_FIG)

data.combined <- RunUMAP(data.combined, dims = 1:17, metric = 'euclidean', 
                         n.neighbors = 50)

plot3         <- DimPlot(data.combined, group.by = 'Experiment', 
                         reduction = 'umap', cols = ColorBlind[c(1,3)], 
                         pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot3, PATH_FIG)

plot4         <- DimPlot(data.combined, group.by = 'orig.ident', 
                         reduction = 'umap', cols = ColorBlind[c(10,1,9,3,6,8)], 
                         pt.size = 0.9, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot4, PATH_FIG)


## AFTER Z-TRANSFORM INTEGRATION -----------------------------------------------  

data  <- cluster_sim_spectrum(data.combined, label_tag='Experiment', 
                              cluster_resolution = 0.6, 
                              corr_method = 'spearman', lambda = 50, 
                              reduction.name = 'cssz', reduction.key = 'CSSZ_')
data  <- RunUMAP(data, reduction = "cssz", 
                 dims = 1:ncol(Embeddings(data, "cssz")), n.neighbors = 50, 
                 reduction.name='umap_cssz', reduction.key='UMAPCSSZ_')

plot5 <- DimPlot(data, group.by = 'Experiment', dim =c(2,1), 
                 reduction = 'umap_cssz', cols = ColorBlind[c(1,3)], 
                 pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP2', y = 'UMAP1')
writePlot(plot5, PATH_FIG)

plot6 <- DimPlot(data, group.by = 'orig.ident', dim =c(2,1), 
                 reduction = 'umap_cssz', cols = ColorBlind[c(10,1,9,3,8,6)], 
                 pt.size = 0.9, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP2', y = 'UMAP1')
writePlot(plot6, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after integration
saveRDS(data, paste0(PATH_SAVE, '/3_Integrated.rds'))
data <- readRDS(paste0(PATH_SAVE, '/3_Integrated.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## CLUSTERING AND OUTLIER ELIMINATION ------------------------------------------
#===============================================================================

data <- FindNeighbors(data, reduction = 'umap_cssz', 
                      dims = 1:ncol(Embeddings(data, 'umap_cssz')))
data <- FindClusters(data)

plot7         <- DimPlot(data, dim =c(2,1), reduction = 'umap_cssz', 
                         group.by = 'seurat_clusters', label = TRUE, 
                         repel = TRUE, cols = c(rep(ColorBlind,5)), 
                         pt.size = 0.5) +
  labs(title = '', x = 'UMAP2', y = 'UMAP1')
writePlot(plot7, PATH_FIG)

# Identification of cell types
cluster.marker <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, 
                                 thresh.use = 0.25)
top30          <- cluster.marker %>% group_by(cluster) %>% top_n(30, avg_log2FC) 
write.csv(top30, paste0(PATH_SAVE, '/ALL_Top30markers.csv'), row.names = FALSE)

# - - - - - - - - - - - - - - - - #
# Based on ImmGen.org annotation  #
# - - - - - - - - - - - - - - - - #
# 17 = Mature ILC2                #
# 19 = Macro/Mono/Granulo         #
# 20 = DC                         #
# 22 = Macro/Mono/Granulo         #
# 23 = Mast Basophils             #
# 24 = B cells                    #
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
    density <- plot_density(data, g, pal = 'magma', size = 2, dim = c(2,1))
    writePlot(density, paste0(PATH_FIG, '/Markers'), 
              paste0('integrated_', g, '_DensityPlot.png'))
    
    feature <- FeaturePlot(data, g, order = T, pt.size = 1, 
                           reduction = 'umap_cssz', dim =c(2,1))
    writePlot(feature, paste0(PATH_FIG, '/Markers'), 
              paste0('integrated_', g, '_FeaturePlot.png'))
  }
}


## DATA FILTERING : NF ONLY ----------------------------------------------------

NF <- subset(data, orig.ident %!in% 'TOX-KO')
NF <- subset(NF, seurat_clusters %!in% c(17,19,22,23,24))

NF <- RunUMAP(NF, reduction = "cssz", dims = 1:ncol(Embeddings(NF, 'cssz')), 
              n.neighbors = 50, reduction.name='umap_cssz2', 
              reduction.key='UMAPCSSZ2_')
NF <- FindNeighbors(NF, reduction = 'umap_cssz2', 
                    dims = 1:ncol(Embeddings(NF, 'umap_cssz2')))
NF <- FindClusters(NF, resolution = 0.06)

# Create inverted embedding to simplify representation
embed           <- as.matrix(data.frame(
  NF@reductions[['umap_cssz2']]@cell.embeddings[,1],
  -NF@reductions[['umap_cssz2']]@cell.embeddings[,2]))

colnames(embed) <- c('REVUMAPCSSZ2_1', 'REVUMAPCSSZ2_2')
NF[['rev_umap_cssz2']]  <- CreateDimReducObject(embeddings = embed, 
                                                key = 'REVUMAPCSSZ2_', 
                                                assay = DefaultAssay(NF))
NF@reductions[['rev_umap_cssz2']]@global <- T


plot8 <- DimPlot(NF, group.by = 'orig.ident', reduction = 'rev_umap_cssz2', 
                 cols = ColorBlind[c(10,1,9,3,6,8)], pt.size = 0.9, 
                 order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot8, PATH_FIG)

plot8b <- DimPlot(NF, group.by = 'orig.ident', reduction = 'rev_umap_cssz2', 
                  cols = ColorBlind[c(11,11,11,3,6)], pt.size = 0.9, 
                  order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot8b, PATH_FIG)

plot9 <- DimPlot(NF, reduction = 'rev_umap_cssz2', group.by = 'seurat_clusters', 
                 label = TRUE, repel = TRUE, 
                 cols = ColorBlind[c(1,7,2,10,9,12)], pt.size = 0.9) + 
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot9, PATH_FIG)

# Draw population markers
markers_genes <- c('Nfil3', 'Flt3', 'Tcf7', 'Zbtb16', 'Tfc7', 'Bcl11b', 'Cd74')

for(g in markers_genes){
  if(g %in% rownames(NF)){
    density <- plot_density(NF, g, pal = 'magma', size = 2)
    writePlot(density, paste0(PATH_FIG, '/Markers'), 
              paste0('NF_filtered_', g, '_DensityPlot.png'))
    
    feature <- FeaturePlot(NF, g, order = T, pt.size = 1, 
                           reduction = 'rev_umap_cssz2')
    writePlot(feature, paste0(PATH_FIG, '/Markers'), 
              paste0('NF_filtered_', g, '_FeaturePlot.png'))
  }
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint of subsets
saveRDS(NF, paste0(PATH_SAVE, '/4_Filtered_NF.rds'))
NF <- readRDS(paste0(PATH_SAVE, '/4_Filtered_NF.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## DATA FILTERING : ALL DATASETS -----------------------------------------------

data <- subset(data, seurat_clusters %!in% c(17,19,22,23,24))

data <- RunUMAP(data, reduction = 'cssz', 
                dims = 1:ncol(Embeddings(data, 'cssz')), n.neighbors = 50, 
                reduction.name='umap_cssz2', reduction.key='UMAPCSSZ2_')
data <- FindNeighbors(data, reduction = 'umap_cssz2', 
                      dims = 1:ncol(Embeddings(data, 'umap_cssz2')))
data <- FindClusters(data, resolution = 0.06)

# Create inverted embedding to simplify representation
embed           <- -data@reductions[['umap_cssz2']]@cell.embeddings
colnames(embed) <- c('REVUMAPCSSZ2_1', 'REVUMAPCSSZ2_2')
data[['rev_umap_cssz2']]  <- CreateDimReducObject(embeddings = embed, 
                                                  key = 'REVUMAPCSSZ2_',
                                                  assay = DefaultAssay(data))
data@reductions[['rev_umap_cssz2']]@global <- T

plot10 <- DimPlot(data, group.by = 'orig.ident', reduction = 'rev_umap_cssz2', 
                  cols = ColorBlind[c(10,1,9,3,6,8)], pt.size = 0.9, 
                  order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot10, PATH_FIG)

plot10b <- DimPlot(data, group.by = 'orig.ident', reduction = 'rev_umap_cssz2', 
                   cols = ColorBlind[c(11,11,11,3,6,8)], pt.size = 0.9, 
                   order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot10b, PATH_FIG)

plot11 <- DimPlot(data, reduction = 'rev_umap_cssz2', 
                  group.by = 'seurat_clusters', label = TRUE, repel = TRUE, 
                  cols = ColorBlind[c(1,7,2,10,9)], pt.size = 0.9) + 
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot11, PATH_FIG)

# Draw population markers
markers_genes <- c('Nfil3', 'Flt3', 'Tcf7', 'Zbtb16', 'Tfc7', 'Bcl11b', 'Cd74')

for(g in markers_genes){
  if(g %in% rownames(NF)){
    density <- plot_density(NF, g, pal = 'magma', size = 2)
    writePlot(density, paste0(PATH_FIG, '/Markers'), 
              paste0('ALL_filtered_', g, '_DensityPlot.png'))
    
    feature <- FeaturePlot(NF, g, order = T, pt.size = 1, 
                           reduction = 'rev_umap_cssz2')
    writePlot(feature, paste0(PATH_FIG, '/Markers'), 
              paste0('ALL_filtered_', g, '_FeaturePlot.png'))
  }
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint of subsets
saveRDS(data, paste0(PATH_SAVE, '/4_Filtered_ALL.rds'))
data <- readRDS(paste0(PATH_SAVE, '/4_Filtered_ALL.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## DC SCORING ------------------------------------------------------------------
#===============================================================================
# Based on markers used in : doi:10.1038/s41590-019-0445-7

## DC SCORING ------------------------------------------------------------------

# Loading a subset of DC and pre-DC markers
DC_table <- read.csv('C:/Users/E15639P/Data/Lists/DC_enrichement_list.csv', 
                     header = TRUE, sep = ';', row.names = 1)
DC_genelist <- list(row.names(DC_table[(DC_table$p_val < 0.01 & 
                                          DC_table$avg_logFC > 0.4),]))

# Adding DC_score as metadata
data <- AddModuleScore(data, DC_genelist, name = 'DC_score')

# Score spatial repartition in cells
plot_density(data, 'DC_score1', pal = 'magma', size = 2)



