#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
# Based on : https://satijalab.org/seurat/archive/v3.1/integration.html
#
# Exploration of datasets containing already described ILC early precursors
#
# Load datasets from a same batch and convert them in a Seurat Object
# Perform QC and filter samples separately
# Lognormalize and scale merged datasets
# Apply cell cycle genes regression during scaling
# Perform clustering and identify cell types using ImmGen.org
# Remove outliers based on this identification



#===============================================================================
## SET UP ----------------------------------------------------------------------
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
source('C:/Users/E15639P/NFIL3_dev_ILC/Custom_Functions.R')
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
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq'), 'WT_No_Integration'))
dir.create(file.path(paste0(PATH, '/Saves'), 'scRNA-seq'))
dir.create(file.path(paste0(PATH, '/Saves/scRNA-seq'), 'WT_No_Integration'))

# Set up path for figures and saves
PATH_FIG  <- paste0(PATH, '/Figures/scRNA-seq/WT_No_Integration')
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/WT_No_Integration')

CellTypeOrder <- c('ALP', 'ILCpro')



#===============================================================================
## INPUT -----------------------------------------------------------------------
#===============================================================================

# Opening the datasets and turning them into Seurat object
# Old Experiments :
ALP       <- Read10X(data.dir = paste0(DATA_DIR, '/1_ABSC_CLP')) %>% 
  CreateSeuratObject(project = 'ALP', min.cells = 1, min.features = 0)
ILCpro    <- Read10X(data.dir = paste0(DATA_DIR, '/3_ABSC_EILP_WT')) %>%
  CreateSeuratObject(project = 'ILCpro', min.cells = 1, min.features = 0)

# Merge samples to perform the quality check 
data  <- merge(ALP, ILCpro, add.cell.ids = c('ALP', 'ILCpro'), 
               project = 'scRNA_No_Integration')
data@meta.data[['Project']] <- 'WT_No_Integration'

# Organize cell type in the development order for the plots
data@active.ident <- factor(data@active.ident, levels = CellTypeOrder)



#===============================================================================
## QUALITY CHECK ---------------------------------------------------------------
#===============================================================================

# Calculate the percent of mitochondrial genes expression
data$Percent.mt <- PercentageFeatureSet(data, pattern = '^mt-')

plot0   <- VlnPlot(data, features = c('nFeature_RNA','Percent.mt'), 
                   group.by = 'orig.ident', cols = ColorBlind[c(10, 9)]) + 
  theme(axis.title.x = element_blank())
writePlot(plot0, PATH_FIG)

# Elimination of cells that over-expressed mitochondrial genes (dead cells)
data.QC <- subset(data, subset = (nFeature_RNA > 200 & 
                                    nFeature_RNA < 5000 & Percent.mt < 2)) 

plot1   <- VlnPlot(data.QC, features = c('nFeature_RNA', 'Percent.mt'), 
                   group.by = 'orig.ident', cols = ColorBlind[c(10, 9)]) + 
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
s.genes       <- capitalize(tolower(cc.genes$s.genes))
g2m.genes     <- capitalize(tolower(cc.genes$g2m.genes))
data.combined <- CellCycleScoring(data.combined, s.features = s.genes, 
                                  g2m.features = g2m.genes, set.ident = TRUE)


# Calculate difference between S and G2M is described as more relevant
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Alternate Workflow
data.combined$CC.Difference <- data.combined$S.Score - data.combined$G2M.Score
data.combined <- ScaleData(data.combined, vars.to.regress = c('CC.Difference'), 
                           verbose = TRUE, features = rownames(data.combined), 
                           do.scale = TRUE)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after scaling step
saveRDS(data.combined, paste0(PATH_SAVE, '/2_Scaled.rds'))
data.combined <- readRDS(paste0(PATH_SAVE, '/2_Scaled.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 


#===============================================================================
## VISUALIZATION ---------------------------------------------------------------
#===============================================================================

data.combined <- FindVariableFeatures(data.combined, nfeatures = 3000)
data.combined <- RunPCA(data.combined, npcs = 50)

plot2         <- ElbowPlot(data.combined, 30) +
  labs(title = "", y = 'Standard deviation', x = 'Component number') + 
  theme(plot.title = element_text(hjust = 0.5))
writePlot(plot2, PATH_FIG)

data.combined <- RunUMAP(data.combined, dims = 1:16, metric = 'euclidean', 
                         n.neighbors = 50)

plot3         <- DimPlot(data.combined, group.by = 'orig.ident', 
                         reduction = 'umap', cols = ColorBlind[c(10,9)], 
                         pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot3, PATH_FIG)



#===============================================================================
## CLUSTERING AND OUTLIER ELIMINATION ------------------------------------------
#===============================================================================

data  <- FindNeighbors(data.combined)
data  <- FindClusters(data)

plot4 <- DimPlot(data, reduction = 'umap', group.by = 'seurat_clusters', 
                 label = TRUE, repel = TRUE, cols = c(rep(ColorBlind,3)), 
                 pt.size = 0.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot4, PATH_FIG)

data  <- subset(data, seurat_clusters %!in% 10)
data  <- RunUMAP(data, dims = 1:16, metric = 'euclidean', n.neighbors = 50)

plot5         <- DimPlot(data, group.by = 'orig.ident', reduction = 'umap', 
                         cols = ColorBlind[c(10,9)], pt.size = 0.5, 
                         order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')
writePlot(plot5, PATH_FIG)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint of subsets
saveRDS(data, paste0(PATH_SAVE, '/3_Filtered.rds'))
data <- readRDS(paste0(PATH_SAVE, '/3_Filtered.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


