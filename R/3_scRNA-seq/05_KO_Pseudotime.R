#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
# Pseudotime analysis integrating NFIL3-KO and TOX-KO ILC intermediates
#
# Load input files generated during previous steps
# Perform pseudotime reconstruction using Slingshot
# Save pseudotime cell order by trajectories in Seurat Object metadata
# Analyse samples distribution inside trajectories
# Draw gene expression variations across pseudotime



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
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq/KO'), 'Curves'))

# Set up path for figures and saves
PATH_FIG  <- paste0(PATH, '/Figures/scRNA-seq/KO')
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/KO')

CellTypeOrder     <- c('ALP', 'TULIP', 'ILCpro', 'overLIP', 'NF-KO', 'TOX-KO')



#===============================================================================
## SLINGSHOT : WT & NF-KO ------------------------------------------------------
#===============================================================================
#Based on : https://bustools.github.io/BUS_notebooks_R/slingshot.html

NF  <- readRDS(paste0(PATH_SAVE, '/4_Filtered_NF.rds'))
red <- 'rev_umap_cssz2'
xlimit <- c(min(NF@reductions[[red]]@cell.embeddings[,1]),
            max(NF@reductions[[red]]@cell.embeddings[,1]))
ylimit <- c(min(NF@reductions[[red]]@cell.embeddings[,2]),
            max(NF@reductions[[red]]@cell.embeddings[,2]))  

DimPlot(NF, reduction = red, group.by = 'seurat_clusters', 
        label = TRUE, repel = TRUE, 
        cols = ColorBlind[c(3,9,2,7,8,10)], pt.size = 0.9) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')

DimPlot(NF, group.by = "orig.ident", reduction = 'rev_umap_cssz2', 
        cols = ColorBlind[c(11,11,11,3,6)], pt.size = 0.9, 
        order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')


## PSEUDOTIME RECONSTRUCTION ---------------------------------------------------
sce                                     <- as.SingleCellExperiment(NF)     
colData(sce)$Seurat_clusters            <- as.character(NF@active.ident)

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', 
                 reducedDim = toupper(red), start.clus = 0, end.clus = c(5,3))


## ADDING PSEUDOTIMES TO DATA --------------------------------------------------
pseudotime                              <- slingPseudotime(sce)
for(i in 1:ncol(pseudotime)){
  NF@meta.data[[paste0('Slingshot_', i)]] <- pseudotime[,i]
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Saving SCE file after processing
saveRDS(sce, paste0(PATH_SAVE, '/Slingshot_NF.rds'))
sce <- readRDS(paste0(PATH_SAVE, '/Slingshot_NF.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## PLOTTING PSEUDOTIME ---------------------------------------------------------

# Specifying gradient colors for plots
gradient         <- hcl.colors(100)

# Plot with curves
curves <- slingCurves(sce, as.df = TRUE) %>%
  dplyr::rename("x" = paste0(str_replace_all(
    toupper(red), "[[:punct:]]", ""), '_1'), 
    "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

# First Trajectory
plot12 <- FeaturePlot(NF, 'Slingshot_1', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            linewidth = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot12, PATH_FIG)

# Second Trajectory
plot13 <- FeaturePlot(NF, 'Slingshot_2', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            linewidth = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot13, PATH_FIG)


## CELLTYPE DISTRIBUTION -------------------------------------------------------

compare <- data.frame(Slingshot = NF$Slingshot_1, Sample = NF$orig.ident)
compare <- compare[order(compare$Slingshot),]
compare <- compare[!is.na(compare$Slingshot),]
compare$Slingshot_order <- c(1:nrow(compare))

# Draw distributed density - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot14 <- ggplot(
  compare, aes(x=Slingshot, y=..scaled.., 
               fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(alpha=0.5, linewidth = 0.6) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,6,8)], name = '')
writePlot(plot14, PATH_FIG)

# Draw cumulated density - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot15 <- ggplot(
  data=compare, aes(x=Slingshot, 
                    fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(adjust=1.5, position="fill", linewidth = 0.5) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,6,8)], name = '')
writePlot(plot15, PATH_FIG)

# Draw linear density - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot16 <- ggplot(
  compare, aes(x=Slingshot, y = 1, 
               colour=factor(Sample, levels = CellTypeOrder))) +
  geom_point(shape=124) +
  theme_classic()+
  scale_color_manual(values = ColorBlind[c(10,1,9,3,6,8)], name = '')
writePlot(plot16, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Saving CDS file after processing and object data after adding metadata
saveRDS(NF, paste0(PATH_SAVE, '/5_Pseudotime_NF.rds'))
NF <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime_NF.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## SLINGSHOT : WT & NF-KO & TOX-KO ---------------------------------------------
#===============================================================================
#Based on : https://bustools.github.io/BUS_notebooks_R/slingshot.html

ALL <- readRDS(paste0(PATH_SAVE, '/4_Filtered_ALL.rds'))
red <- 'rev_umap_cssz2'
xlimit <- c(min(ALL@reductions[[red]]@cell.embeddings[,1]),
            max(ALL@reductions[[red]]@cell.embeddings[,1]))
ylimit <- c(min(ALL@reductions[[red]]@cell.embeddings[,2]),
            max(ALL@reductions[[red]]@cell.embeddings[,2]))      

DimPlot(ALL, reduction = red, group.by = 'seurat_clusters', 
        label = TRUE, repel = TRUE, 
        cols = ColorBlind[c(3,9,2,7,8,10)], pt.size = 0.9) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')

DimPlot(ALL, group.by = "orig.ident", reduction = 'rev_umap_cssz2', 
        cols = ColorBlind[c(11,11,11,3,6,8)], pt.size = 0.9, 
        order = rev(CellTypeOrder)) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2')


## PSEUDOTIME RECONSTRUCTION ---------------------------------------------------
sce                                     <- as.SingleCellExperiment(ALL)     
colData(sce)$Seurat_clusters            <- as.character(ALL@active.ident)

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', 
                 reducedDim = toupper(red), start.clus = 0, end.clus = c(4,3))


## ADDING PSEUDOTIMES TO DATA --------------------------------------------------
pseudotime                              <- slingPseudotime(sce)
for(i in 1:ncol(pseudotime)){
  ALL@meta.data[[paste0('Slingshot_', i)]] <- pseudotime[,i]
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Saving SCE file after processing
saveRDS(sce, paste0(PATH_SAVE, '/Slingshot_ALL.rds'))
sce <- readRDS(paste0(PATH_SAVE, '/Slingshot_ALL.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## PLOTTING PSEUDOTIME ---------------------------------------------------------

# Specifying gradient colors for plots
gradient         <- hcl.colors(100)

# Plot with curves
curves <- slingCurves(sce, as.df = TRUE) %>%
  dplyr::rename("x" = paste0(str_replace_all(
    toupper(red), "[[:punct:]]", ""), '_1'), 
    "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

# First Trajectory
plot17 <- FeaturePlot(ALL, 'Slingshot_1', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            linewidth = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot17, PATH_FIG)

# Second Trajectory
plot18 <- FeaturePlot(ALL, 'Slingshot_2', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            linewidth = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot18, PATH_FIG)


## CELLTYPE DISTRIBUTION -------------------------------------------------------

compare <- data.frame(Slingshot = ALL$Slingshot_2, Sample = ALL$orig.ident)
compare <- compare[order(compare$Slingshot),]
compare <- compare[!is.na(compare$Slingshot),]
compare$Slingshot_order <- c(1:nrow(compare))


# Draw distributed density - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot19 <- ggplot(
  compare, aes(x=Slingshot, y=..scaled.., 
               fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(alpha=0.5, linewidth = 0.6) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,6,8)], name = '')
writePlot(plot19, PATH_FIG)

# Draw cumulated density - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot20 <- ggplot(
  data=compare, aes(x=Slingshot, 
                    fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(adjust=1.5, position="fill", linewidth = 0.5) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,6,8)], name = '')
writePlot(plot20, PATH_FIG)

# Draw linear density - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plot21 <- ggplot(
  compare, aes(x=Slingshot, y = 1, 
               colour=factor(Sample, levels = CellTypeOrder))) +
  geom_point(shape=124) +
  theme_classic()+
  scale_color_manual(values = ColorBlind[c(10,1,9,3,6,8)], name = '')
writePlot(plot21, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Saving CDS file after processing and object data after adding metadata
saveRDS(ALL, paste0(PATH_SAVE, '/5_Pseudotime_ALL.rds'))
ALL <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime_ALL.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## GENE EXPRESSION VISUALIZATION : WT & NF-KO ----------------------------------
#===============================================================================

red <- 'rev_umap_cssz2'

## WT vs NF-KO -----------------------------------------------------------------
NF <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime_NF.rds'))

for(i in 1:length(NF@meta.data[["orig.ident"]])){
  if(NF@meta.data[["orig.ident"]][i] %!in% 'NF-KO'){
    NF@meta.data[["Sample"]][i] <- 'WT'
  }else{
    NF@meta.data[["Sample"]][i] <- 'NF-KO'
  }
}

plot22 <- DimPlot(NF, reduction = red, group.by='Sample', 
                  cols = ColorBlind[c(6,11)], pt.size = 1)
writePlot(plot22, PATH_FIG)

# Generate table with cells ordered by pseudotime
table_NF   <- OrderMatrix(NF, 'Slingshot_1', min.cells = 50, slot = 'scaledata')


# Draw gene expression across pseudotime - - - - - - - - - - - - - - - - - - - -
DrawExpr(table_NF, bin.number = 11, by.order = F, feature = 'Nfil3', 
         scale = T, std = F, ylim = c(-0.5,1.75), col = 'black', lwd = 2,
         superposed = F)

# Draw gene expression subsetted by sample - - - - - - - - - - - - - - - - - - -
# Add Sample row 
Sample <- c()
Ann <- data.frame(Cell.ID = NF@assays[["RNA"]]@counts@Dimnames[[2]], 
                  Sample = NF@meta.data[["Sample"]])
for(i in 1:nrow(table_NF)){
  x <- ifelse(Ann$Sample[Ann$Cell.ID %in% rownames(table_NF)[i]] %in% 'NF-KO', 'NF-KO', 'WT')
  Sample <- c(Sample, x)
}
table_NF <- cbind(Sample = as.vector(Sample),table_NF)

gene = 'Nfil3'
plot <- CompareExpr(x = table_NF, bin.number = 11, feature = gene, by.order = F, 
                    scale = F, col = ColorBlind[c(8,3)], main = gene, 
                    ylim = c(-0.5,1.75))
writePlot(plot, paste0(PATH_FIG, '/Curves'), filename = paste0('NF_', gene))



#===============================================================================
## GENE EXPRESSION VISUALIZATION : ALL -----------------------------------------
#===============================================================================

red <- 'rev_umap_cssz2'

## ALL LT+ POPULATIONS ---------------------------------------------------------
ALL     <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime_ALL.rds'))
overALL <- subset(ALL, orig.ident %in% c('overLIP', 'NF-KO', 'TOX-KO'))

plot23 <- DimPlot(overALL, reduction = red, group.by='orig.ident', 
                  cols = ColorBlind[c(3,6,8)], pt.size = 1)
writePlot(plot23, PATH_FIG)

# Generate table with cells ordered by pseudotime and add bins 
table_ALL   <- OrderMatrix(overALL, 'Slingshot_2' ,min.cells = 50,  slot = 'scaledata')


# Draw gene expression by sample - - - - - - - - - - - - - - - - - - - - - - - -
# Add Sample row 
Sample <- c()
Ann <- data.frame(Cell.ID = overALL@assays[["RNA"]]@counts@Dimnames[[2]], Sample = overALL@meta.data[["orig.ident"]])
for(i in 1:nrow(table_ALL)){
  if(Ann$Sample[Ann$Cell.ID %in% rownames(table_ALL)[i]] %in% 'NF-KO'){
    Sample <- c(Sample, 'NF-KO')
  }else if(Ann$Sample[Ann$Cell.ID %in% rownames(table_ALL)[i]] %in% 'TOX-KO'){
    Sample <- c(Sample, 'TOX-KO')
  }else{
    Sample <- c(Sample, 'WT')
  }
}
table_ALL  <- cbind(Sample = as.vector(Sample), table_ALL )


genelist <- c('Id2', 'Itgb7', 'Gata3', 'Nfkb1', 'Pdcd1', 'Runx3', 'Sox5', 'Tcf7', 'Tox',
          'Zbtb16', 'Hoxa9', 'Irf8', 'Lmo2', 'Zeb2', 'Nfil3')

for(gene in genelist){
  CompareExpr(x = table_ALL , bin.number = 11, feature = gene, by.order = F, 
              scale = T, col = ColorBlind[c(6,8,3)], main = gene, ylim = c(-0.5,2.5))
  plot <- recordPlot()
  writePlot(plot, paste0(PATH_FIG, '/Curves'), filename = paste0('overALL_', gene),
            height = 8, width = 9 )
}


