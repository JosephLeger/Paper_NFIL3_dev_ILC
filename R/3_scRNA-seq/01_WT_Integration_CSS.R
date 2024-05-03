#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
#Based on : https://satijalab.org/seurat/archive/v3.1/integration.html
#           https://github.com/quadbiolab/simspec/blob/master/vignette/vignette.md

# First step for a new and better integration of all WT datasets
# Integration used here is based on Clustering Similarity Spectrum (CSS) method
# After QC, Integration and outlier elimination, data are clustered 
# Result file is ready for pseudotime analysis on Step_02



#===============================================================================
## SETUP -----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

PATH      <- 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC'
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/WT')
setwd(PATH)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(simspec))
#suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(Nebulosa))

source("C:/Users/E15639P/Doctorat/scRNA-seq/Usefull_functions.R")
CellTypeOrder     <- c("ALP", "overLIP", "TULIP", "ILCpro")



#===============================================================================
## INPUT -----------------------------------------------------------------------
#===============================================================================

dir.data <- "C:/Users/E15639P/Data/scRNA-seq/SingleCell_mm39"

# Opening the datasets and changing them into Seurat object
# Old Experiments :
ALP       <- Read10X(data.dir = paste(dir.data, "/1_ABSC_CLP", sep = '')) %>% 
  CreateSeuratObject(project = "ALP", min.cells = 1, min.features = 0)
ILCpro    <- Read10X(data.dir = paste(dir.data, "/3_ABSC_EILP_WT", sep = '')) %>%
  CreateSeuratObject(project = "ILCpro", min.cells = 1, min.features = 0)
TULIP     <- Read10X(data.dir = paste(dir.data, "/2_ABSC_TULIP", sep = '')) %>% 
  CreateSeuratObject(project = "TULIP", min.cells = 1, min.features = 0)

# Recent Experiment :    
overLIP   <- Read10X(data.dir = paste(dir.data, "/2019_scRNAseq_data/WT", sep = '')) %>%
  CreateSeuratObject(project = "overLIP", min.cells = 1, min.features = 0)


# Spread datasets in groups of experiment
ALP@meta.data[["Experiment"]]         <- "Experiment 1"
ILCpro@meta.data[["Experiment"]]      <- "Experiment 1"
TULIP@meta.data[["Experiment"]]       <- "Experiment 1"
overLIP@meta.data[["Experiment"]]     <- "Experiment 2"

# Merge all the experiments to perform the quality check 
data  <- merge(TULIP, c(ILCpro, ALP, overLIP), add.cell.ids = c("TULIP", "ILCpro", "ALP", "overLIP"), project = "scRNA_Integration")
data@meta.data[["Project"]] <- "scRNA_WT_Integration"

# Organize cell type in the development order for the plots
data@active.ident <- factor(data@active.ident, levels = CellTypeOrder)



#===============================================================================
## QUALITY CHECK ---------------------------------------------------------------
#===============================================================================

# Calculate the percent of mitochondrial genes expression
data$Percent.mt <- PercentageFeatureSet(data, pattern = "^mt-")

plot0   <- VlnPlot(data, features = c("nFeature_RNA","Percent.mt"), group.by = "Experiment", cols = ColorBlind[c(1, 3)]) & 
  theme(axis.title.x = element_blank())
plot0

# Elimination of cells that over-expressed mitochondrial genes (dead cells)
data.QC <- subset(data, subset = ((Experiment == 'Experiment 1' & nFeature_RNA > 200 & nFeature_RNA < 5000 & Percent.mt < 2) |
                                  (Experiment != 'Experiment 1' & nFeature_RNA > 1200 & nFeature_RNA < 5000 & Percent.mt < 4))) 

plot1   <- VlnPlot(data.QC, features = c("nFeature_RNA", "Percent.mt"), group.by = "Experiment", cols = ColorBlind[c(1, 3)]) & 
  theme(axis.title.x = element_blank())
plot1


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
data.combined <- CellCycleScoring(data.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


# Calculate difference between S and G2M is described as more relevant
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Alternate Workflow
data.combined$CC.Difference <- data.combined$S.Score - data.combined$G2M.Score
data.combined <- ScaleData(data.combined, vars.to.regress = c("CC.Difference","Percent.mt"), verbose = TRUE, features = rownames(data.combined), do.scale = TRUE)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after scaling step
#saveRDS(data.combined, paste0(PATH_SAVE, '/2_Scaled.rds'))
data.combined <- readRDS(paste0(PATH_SAVE, '/2_Scaled.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 


#===============================================================================
## CLUSTERING SIMILARITY SPECTRUM INTEGRATION ----------------------------------
#===============================================================================
# To reorder color plotting
CellTypeOrder     <- c("ILCpro", "ALP", "TULIP", "overLIP")
data.combined@active.ident <- factor(data.combined@active.ident, levels = CellTypeOrder)

# Before Integration - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

data.combined <- FindVariableFeatures(data.combined, nfeatures = 3000)
data.combined <- RunPCA(data.combined, npcs = 50)

plot2         <- ElbowPlot(data.combined, 30) +
  labs(title = "", y = "Écart type", x = "Numéro de la composante") + 
  theme(plot.title = element_text(hjust = 0.5))
plot2

data.combined <- RunUMAP(data.combined, dims = 1:16, metric = "euclidean", n.neighbors = 50)

plot3         <- DimPlot(data.combined, group.by = "Experiment", reduction = "umap", cols = ColorBlind[c(1,10)], pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = "", x = "UMAP1", y = "UMAP2")
plot3

plot4         <- DimPlot(data.combined, group.by = "orig.ident", reduction = "umap", cols = ColorBlind[c(9,10,1,3)], pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = "", x = "UMAP1", y = "UMAP2")
plot4


# After Integration Z-Transform  - - - - - - - - - - - - - - - - - - - - - - - -  

data <- cluster_sim_spectrum(data.combined, label_tag="Experiment", cluster_resolution = 0.6, corr_method = "spearman", lambda = 50, reduction.name = "cssz", reduction.key = "CSSZ_")
data <- RunUMAP(data, reduction = "cssz", dims = 1:ncol(Embeddings(data, "cssz")), n.neighbors = 50, reduction.name="umap_cssz", reduction.key="UMAPCSSZ_")

plot5         <- DimPlot(data, group.by = "Experiment", reduction = "umap_cssz", cols = ColorBlind[c(1,10)], pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = "", x = "UMAP1", y = "UMAP2") + scale_x_reverse()
plot5

plot6         <- DimPlot(data, group.by = "orig.ident", reduction = "umap_cssz", cols = ColorBlind[c(9,10,1,3)], pt.size = 0.5, order = rev(CellTypeOrder)) +
  labs(title = "", x = "UMAP1", y = "UMAP2") + scale_x_reverse()
plot6

FeaturePlot(data, "S.Score", reduction = 'umap_cssz') + scale_x_reverse()
FeaturePlot(data, "G2M.Score", reduction = 'umap_cssz') + scale_x_reverse()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint after integration
#saveRDS(data, paste0(PATH_SAVE, '/3_Integrated.rds'))
data <- readRDS(paste0(PATH_SAVE, '/3_Integrated.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## CLUSTERING AND OUTLIER ELIMINATION ------------------------------------------
#===============================================================================

data <- FindNeighbors(data, reduction = "umap_cssz", dims = 1:ncol(Embeddings(data, "umap_cssz")))
data <- FindClusters(data)

plot7         <- DimPlot(data, reduction = "umap_cssz", group.by = "seurat_clusters", label = TRUE, repel = TRUE, cols = c(rep(ColorBlind,10)), pt.size = 0.5) +
  labs(title = "", x = "UMAP1", y = "UMAP2") + scale_x_reverse()
plot7


# Identification of cell types
cluster.marker <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top30          <- cluster.marker %>% group_by(cluster) %>% top_n(30, avg_log2FC) 
#write.csv(top30, paste0(PATH_SAVE, '/WT_Top30markers.csv'), row.names = FALSE)

# - - - - - - - - - - - - - - - - #
# 18 = ILC2 (Icos, Bcl11b)        #
# 19 = B cells                    #
# 20 = Mast Basophils             #
# 21 = ILC2 (Icos, Bcl11b, Rora)  #
# 22 = Macro/Mono/Granulo         #
# 23 = Macro/Mono/Granulo         #
# - - - - - - - - - - - - - - - - #


## KEY GENES REPRESENTATION ----------------------------------------------------

plot_density(data, 'Nfil3', pal = "magma") + scale_x_reverse()

# B cells
plot_density(data, 'Pax5', pal = "magma") + scale_x_reverse()
plot_density(data, 'Ebf1', pal = "magma") + scale_x_reverse()
# Mast cells
plot_density(data, 'Cpa3', pal = "magma") + scale_x_reverse()
plot_density(data, 'Gzmb', pal = "magma") + scale_x_reverse()
plot_density(data, 'Kit', pal = "magma") + scale_x_reverse()
# Mono/Macro
plot_density(data, 'Lgals3', pal = "magma") + scale_x_reverse()
plot_density(data, 'Lyz2', pal = "magma") + scale_x_reverse()
# ILC2
plot_density(data, 'Il2ra', pal = "magma") + scale_x_reverse()
plot_density(data, 'Icos', pal = "magma") + scale_x_reverse()
plot_density(data, 'Stab2', pal = "magma") + scale_x_reverse()
plot_density(data, 'Bcl11b', pal = "magma") + scale_x_reverse()
plot_density(data, 'Tox2', pal = "magma") + scale_x_reverse()

FeaturePlot(data, 'Nfil3', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Pax5', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Ebf1', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Cpa3', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Gzmb', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Kit', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Lgals3', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Lyz2', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Il2ra', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Icos', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Stab2', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Bcl11b', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Tox2', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()
FeaturePlot(data, 'Il7r', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()

FeaturePlot(data, 'Foxo1', order = T, pt.size = 1, reduction = "umap_cssz") + scale_x_reverse()


## CONTAMINANT PERCENTAGE ------------------------------------------------------

# Overall outliers
Outliers <- subset(data, seurat_clusters %in% c(18,19,20,22,23))
Inliers  <- subset(data, seurat_clusters %!in% c(18,19,20,21,22,23))

prop_out <- table(factor(Outliers$orig.ident, levels = CellTypeOrder))
prop_in  <- table(factor(Inliers$orig.ident, levels = CellTypeOrder))
prop_all <- (prop_out/(prop_in+prop_out))*100

plot8    <- barplot(prop_all[c(2,3)], col = ColorBlind[c(3,1)], border = 'white', 
                    ylab = 'Outlier percentage', main = 'Percentage of outliers' , ylim = c(0, 60))
y        <- prop_all[c(2,3)]
text(plot8, y+2, labels = paste(format(y, scientific=FALSE, digits = 2 ), "%"))
plot8    <- recordPlot()


# Outliers by cellype
prop_table <- as.data.frame(matrix(0,4,4))
colnames(prop_table) <- c('ILC2', 'B_cells', 'Mast_Baso', 'Macro_Mono_Granulo')
rownames(prop_table) <- names(table(Outliers$orig.ident))

for(i in 1:length(Outliers$seurat_clusters)){
  if(Outliers$seurat_clusters[i] %in% c(18,21)){
    prop_table[Outliers$orig.ident[i],'ILC2'] <- prop_table[Outliers$orig.ident[i],'ILC2']+1
  }else if(Outliers$seurat_clusters[i] %in% 19){
    prop_table[Outliers$orig.ident[i],'B_cells'] <- prop_table[Outliers$orig.ident[i],'B_cells']+1
  }else if(Outliers$seurat_clusters[i] %in% 20){
    prop_table[Outliers$orig.ident[i],'Mast_Baso'] <- prop_table[Outliers$orig.ident[i],'Mast_Baso']+1
  }else if(Outliers$seurat_clusters[i] %in% c(22,23)){
    prop_table[Outliers$orig.ident[i],'Macro_Mono_Granulo'] <- prop_table[Outliers$orig.ident[i],'Macro_Mono_Granulo']+1
  }
}

prop_table <- rbind(prop_table, Total = colSums(prop_table))
perc_table <- prop_table[1:4,]/prop_table[rep(5,4),]*100
bar_table  <- data.frame(Percentage = c(perc_table[,1],perc_table[,2], perc_table[,3], perc_table[,4]),
                         CellType = rep(colnames(perc_table), each = 4),
                         Sample = rep(rownames(perc_table)))

plot9 <- ggplot(bar_table[bar_table$Sample %in% c('TULIP', 'overLIP'),], aes(fill=Sample, y=Percentage, x=CellType)) + 
  geom_bar(position='stack', stat='identity') + ggtitle('Proportions of celltypes in outleirs') +
  scale_fill_manual(values= ColorBlind[c(3,1)]) +
  theme_classic()
plot9


## DATA FILTERING --------------------------------------------------------------

# Remove mature cells (18), added B cells (19) and contaminants (20,22,23)
data <- Inliers

data <- RunUMAP(data, reduction = "cssz", dims = 1:ncol(Embeddings(data, "cssz")), n.neighbors = 50, reduction.name="umap_cssz2", reduction.key="UMAPCSSZ2_")
data <- FindNeighbors(data, reduction = "umap_cssz2", dims = 1:ncol(Embeddings(data, "umap_cssz2")))
data <- FindClusters(data, resolution = 0.1)


# Create inverted embedding to simplify representation
embed           <- as.matrix(data.frame(-data@reductions[["umap_cssz2"]]@cell.embeddings[,1],
                                        data@reductions[["umap_cssz2"]]@cell.embeddings[,2]))
colnames(embed) <- c('REVUMAPCSSZ2_1', 'REVUMAPCSSZ2_2')
data[['rev_umap_cssz2']]  <- CreateDimReducObject(embeddings = embed, key = "REVUMAPCSSZ2_", assay = DefaultAssay(data))
data@reductions[["rev_umap_cssz2"]]@global <- T


plot10 <- DimPlot(data, group.by = "orig.ident", reduction = 'rev_umap_cssz2', 
                 cols = ColorBlind[c(9,10,1,3)], pt.size = 1, 
                 order = rev(CellTypeOrder)) + 
  labs(title = "", x = "UMAP1", y = "UMAP2")
plot10

plot11 <- DimPlot(data, reduction = 'rev_umap_cssz2', group.by = "seurat_clusters", 
                 label = TRUE, repel = TRUE, 
                 cols = c(ColorBlind[c(3,9,2,7,8,10)], "purple"), pt.size = 0.9) +
  labs(title = "", x = "UMAP1", y = "UMAP2")
plot11

plot_density(data, 'Nfil3', pal = "magma")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# CheckPoint of subsets
#saveRDS(data, paste0(PATH_SAVE, '/4_Filtered.rds'))
data <- readRDS(paste0(PATH_SAVE, '/4_Filtered.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## DC SCORING AND GENE VISUALIZATION -------------------------------------------
#===============================================================================
# Based on markers used in : doi:10.1038/s41590-019-0445-7

## DC SCORING ------------------------------------------------------------------

# Loading a subset of DC and pre-DC markers
DC_table <- read.csv("C:/Users/E15639P/Data/Lists/DC_enrichement_list.csv", header = TRUE, sep = ";", row.names = 1)
DC_genelist <- list(row.names(DC_table[(DC_table$p_val < 0.01 & DC_table$avg_logFC > 0.4),]))

# Adding DC_score as metadata
data <- AddModuleScore(data, DC_genelist, name = 'DC_score')

# Score spatial repartition in cells
plot_density(data, 'DC_score1', pal = 'magma', size = 2)


## GENE VISUALIZATION ----------------------------------------------------------

plot_density(data, 'Nfil3', pal = 'magma', size = 2)
plot_density(data, 'Bcl11b', pal = 'magma', size = 2)
plot_density(data, 'Cd74', pal = 'magma', size = 2)
plot_density(data, 'Flt3', pal = 'magma', size = 2)
plot_density(data, 'Tcf7', pal = 'magma', size = 2)
plot_density(data, 'Zbtb16', pal = 'magma', size = 2)

# Expressed BZIPs
plot_density(data, 'Nfil3', pal = "magma", size = 2)
plot_density(data, 'Cebpa', pal = "magma", size = 2)
plot_density(data, 'Batf3', pal = "magma", size = 2)
plot_density(data, 'Creb3l2', pal = "magma", size = 2)
FeaturePlot(data, 'Nfil3', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Cebpa', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Batf3', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Creb3l2', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)

# Other Bzips
plot_density(data, 'Atf1', pal = 'magma', size = 2)
plot_density(data, 'Bach2', pal = 'magma', size = 2)
plot_density(data, 'Cebpb', pal = 'magma', size = 2)
plot_density(data, 'Jun', pal = 'magma', size = 2)
plot_density(data, 'Maf', pal = 'magma', size = 2)
plot_density(data, 'Mafg', pal = 'magma', size = 2)
plot_density(data, 'Mafk', pal = 'magma', size = 2)
plot_density(data, 'Tsc22d1', pal = 'magma', size = 2)

# Same K-mean cluster
plot_density(data, 'Aff3', pal = "magma", size = 2)
plot_density(data, 'Nfkb1', pal = "magma", size = 2)
plot_density(data, 'Runx3', pal = "magma", size = 2)
plot_density(data, 'Sox5', pal = "magma", size = 2)
plot_density(data, 'Dach1', pal = "magma", size = 2)
plot_density(data, 'Mycn', pal = "magma", size = 2)
FeaturePlot(data, 'Aff3', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Nfkb1', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Runx3', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Sox5', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Dach1', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Mycn', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)




plot_density(data, 'Id2', pal = "magma", size = 2)
plot_density(data, 'Tox', pal = "magma", size = 2)
plot_density(data, 'Gata3', pal = "magma", size = 2)
plot_density(data, 'Batf', pal = "magma", size = 2)
plot_density(data, 'Zeb2', pal = "magma", size = 2)
FeaturePlot(data, 'Id2', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Tox', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Gata3', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Batf', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)
FeaturePlot(data, 'Zeb2', reduction = 'rev_umap_cssz2', pt.size = 2, order = T)

