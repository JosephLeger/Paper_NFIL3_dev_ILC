#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================




#===============================================================================
## SET UP ----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PROJECT INFO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

PATH <- 'C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC'

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

setwd(PATH)

# Load Packages and custom functions
source('C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Scripts/Custom_Functions.R')
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
## SLINGSHOT : NF --------------------------------------------------------------
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
            size = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot12, PATH_FIG)

# Second Trajectory
plot13 <- FeaturePlot(NF, 'Slingshot_2', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            size = 1.5) +
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
## SLINGSHOT : ALL -------------------------------------------------------------
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
            size = 1.5) +
  labs(title = '', x = 'UMAP1', y = 'UMAP2', pt.size = 0.9) + 
  scale_color_gradientn(name = 'Pseudotime', colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
writePlot(plot17, PATH_FIG)

# Second Trajectory
plot18 <- FeaturePlot(ALL, 'Slingshot_2', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), 
            size = 1.5) +
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
## GENE EXPRESSION VISUALIZATION : NF ------------------------------------------
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
         scale = F, std = F, ylim = c(-0.5,1.75), col = c("#26C4EC"), lwd = 3,
         superposed = F, compare.with = F, write = T, 
         dir = paste0(PATH_FIG, '/Curves'))


# Draw gene expression subsetted by sample - - - - - - - - - - - - - - - - - - -
# Add Sample row 
Sample <- c()
Ann <- data.frame(Cell.ID = NF@assays[["RNA"]]@counts@Dimnames[[2]], Sample = NF@meta.data[["Sample"]])
for(i in 1:ncol(table_NF)){
  x <- ifelse(Ann$Sample[Ann$Cell.ID %in% colnames(table_NF)[i]] %in% 'NF-KO', 1, 0)
  Sample <- c(Sample, x)
}
table_NF <- rbind(Sample = as.vector(Sample),table_NF)

gene = 'Nfil3'
plot <- CompareExpr(x = table_NF, bin.number = 11, feature = gene, by.order = F, 
                    scale = F, col = ColorBlind[c(11, 6)], main = gene, 
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
  plot <- CompareExpr(x = table_ALL , bin.number = 11, feature = gene, by.order = F, 
                      scale = F, col = ColorBlind[c(6,8,3)], main = gene, ylim = c(-0.5,2.5))
  writePlot(plot, paste0(PATH_FIG, '/Curves'), filename = paste0('overALL_', gene),
            height = 8, width = 9 )
}


###
## ALTERNATIVE FIGURES
###

gene <- 'Tox'

# Draw curves with corresponding points - - - - - - - - - - - - - - - - - - - -

table_rep   <- PseudotimeRepartition(table_ALL, 11) 

# With exact points
CompareExpr(x = table_ALL , bin.number = 11, feature = gene, by.order = F, lwd =3, 
            scale = F, col = ColorBlind[c(6,8,3)], main = gene, ylim = c(-1,3))
par(new=TRUE)
plot(x=as.numeric(table_rep$Pseudotime), y=as.numeric(table_rep[,gene]), 
     pch=19, col = ColorBlind[c(6,8,3)][factor(table_rep$Sample)], 
     cex = 0.33, xlab='', ylab='', ylim = c(-1,3))
par(new=TRUE)
CompareExpr(x = table_ALL , bin.number = 11, feature = gene, by.order = F, lwd =3, 
            scale = F, col = ColorBlind[c(6,8,3)], main = gene, ylim = c(-1,3))

# With binned points
CompareExpr(x = table_ALL , bin.number = 11, feature = gene, by.order = F, lwd =3, 
            scale = F, col = ColorBlind[c(6,8,3)], main = gene, ylim = c(-1,3))
par(new=TRUE)
plot(x=as.numeric(table_rep$Bin), y=as.numeric(table_rep[,gene]), 
     pch=19, col = ColorBlind[c(6,8,3)][factor(table_rep$Sample)], 
     cex = 0.33, xlab='', ylab='', ylim = c(-1,3),
     xlim = c(min(table_rep$Pseudotime), max(table_rep$Pseudotime)))
par(new=TRUE)
CompareExpr(x = table_ALL , bin.number = 11, feature = gene, by.order = F, lwd =3, 
            scale = F, col = ColorBlind[c(6,8,3)], main = gene, ylim = c(-1,3))


# Get statistical confidence for comparison - - - - - - - - - - - - - - - - - -

table_rep     <- PseudotimeRepartition(table_ALL, 11) 
gene          <- 'Tox'

# Generate empty final table
columns       <- c()
samp_list     <- factor(c('WT', 'NF-KO', 'TOX-KO'), levels = c('WT', 'NF-KO', 'TOX-KO'))
for(i in 1:(length(samp_list)-1)){
    for(j in (i+1):length(samp_list)){
      columns <- c(columns, paste0(samp_list[i], '_vs_', samp_list[j], '_FC'))
      columns <- c(columns, paste0(samp_list[i], '_vs_', samp_list[j], '_pval'))
    }
}

table_stat    <- as.data.frame(matrix(NA, length(unique(table_rep$Bin)), length(columns)))
colnames(table_stat) <- columns


# Fill created table
# For each bin present in the table
for(b in 1:length(unique(table_rep$Bin))){
  bin         <- sort(unique(table_rep$Bin))[b]
  sub         <- table_rep[table_rep$Bin %in% bin,]
  newrow      <- c()
  # For each non-mirrored pair of sample
  for(i in 1:(length(samp_list)-1)){
    for(j in (i+1):length(samp_list)){
      # If both sample are present in current bin, get their values
      if(samp_list[i] %in% sub$Sample & samp_list[j] %in% sub$Sample){
        i_val <- sub[sub$Sample %in% samp_list[i], gene]
        j_val <- sub[sub$Sample %in% samp_list[j], gene]
        # If enough values are present, calculate FC and pval using Mann-Whitney
        if(length(i_val) > 5 & length(j_val) > 5 ){
          wtest <- wilcox.test(as.numeric(i_val), as.numeric(j_val), conf.level = 0.99)
          newrow <- c(newrow, mean(i_val)/mean(j_val))
          newrow <- c(newrow, wtest[["p.value"]])
        }else{
          newrow <- c(newrow, NA, NA)
        }
      }else{
        newrow <- c(newrow, NA, NA)
      }
      
    }
  }
  # Fill empty table row by row
  table_stat[b,] <- newrow
}









#################################################################################

# Draw all candidates
Candidates.UP   <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_UP.txt')[,1]
Candidates.DOWN <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_DOWN.txt')[,1]

for(g in Candidates.UP){
  #pdf(paste0('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Figures/scRNA-seq/KO/Candidates/UP/', g, '.pdf'),
  #width = 9.00, height = 9.00)
  CompareExpr(x = table_ALL , bin.number = 11, feature = g, by.order = F, 
              scale = F, col = ColorBlind[c(3,6,8)], main = g, ylim = c(-0.5,1.75))
  #dev.off()
}

for(g in Candidates.DOWN){
  #pdf(paste0('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Figures/scRNA-seq/KO/Candidates/DOWN/', g, '.pdf'),
  #    width = 9.00, height = 9.00)
  CompareExpr(x = table_ALL , bin.number = 11, feature = g, by.order = F, 
              scale = F, col = ColorBlind[c(3,6,8)], main = g, ylim = c(-0.5,1.75))
  #dev.off()
}


## Candidates
TF_list <- read.table('C:/Users/E15639P/Data/Lists/masterTFlist.txt', header = T)[,1]


for(g in c(Candidates.UP, Candidates.DOWN, 'Nfil3')){
  if(g %in% c(TF_list, List)){
    pdf(paste0('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Figures/scRNA-seq/KO/Candidates/TF/', g, '.pdf'),
        width = 9.00, height = 9.00)
    CompareExpr(x = table_ALL , bin.number = 11, feature = g, by.order = F, 
                scale = F, col = ColorBlind[c(3,6,8)], main = g, ylim = c(-0.5,1.75))
    dev.off()
  }
  
}



                                                                                             


