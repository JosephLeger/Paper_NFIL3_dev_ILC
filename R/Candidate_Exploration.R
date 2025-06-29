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

PATH <- 'C:/Users/E15639P/Desktop/Review'

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

setwd(PATH)

# Load Packages and custom functions
source('C:/Users/E15639P/Doctorat/Custom_Functions.R')
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(VennDiagram))


# Create required directories
dir.create(file.path(paste0(PATH, '/Figures'), 'Candidates'))
dir.create(file.path(paste0(PATH, '/Figures/Candidates'), 'Curves'))
dir.create(file.path(paste0(PATH, '/Figures/Candidates'), 'Dotplots'))
dir.create(file.path(paste0(PATH, '/Figures/Candidates/Curves'), 'WT'))
dir.create(file.path(paste0(PATH, '/Figures/Candidates/Curves'), 'KO'))
dir.create(file.path(paste0(PATH, '/Saves'), 'Candidates'))

# Set up path for figures and saves
PATH_FIG  <- paste0(PATH, '/Figures/Candidates')
PATH_SAVE <- paste0(PATH, '/Saves/Candidates')

CellTypeOrder     <- c('ALP', 'TULIP', 'ILCpro', 'overLIP', 'NF-KO', 'TOX-KO')



#===============================================================================
## DEFINE USEFULL FUNCIONS -----------------------------------------------------
#===============================================================================

display_venn <- function(x, ...){
  
  # Draw a Venn Diagram based on a list of categorized elements
  #
  # x = list of categorized elements to cross
  
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

#-------------------------------------------------------------------------------

Regions_PieChart <- function(x){
  
  # Draw a piechart of regions distribution from a peaks annotated table
  #
  # x = HOMER peaks annotated table
  
  regions     <- capitalize(unlist(lapply(
    strsplit(as.character(x$Annotation), " \\("),"[[",1)))
  y           <- as.data.frame(table(regions))
  colnames(y) <- c('Region', 'value')
  plot <- ggplot(y, aes(x='', y=value, fill=Region)) +
    geom_col(color='black') +
    geom_text(aes(x = 1.6, label=value),
              position=position_stack(vjust=0.5)) +
    scale_fill_manual(values=ColorBlind[c(3,2,7,8,9,10)]) +
    coord_polar(theta = 'y') +
    theme_void()
  print(plot)
  return(plot)
}



#===============================================================================
## LOAD NECESSARY FILES --------------------------------------------------------
#===============================================================================

# Bulk RNA-seq NFIL3 table (GSE291075_Table_Norm_Stats_Batch1-2.txt)
NF_Table   <- read.table('C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Saves/Bulk_RNA-seq/Nfil3/Table_Results.txt')
NF_META    <- read.table('C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Data/Bulk_RNA-seq/Nfil3/SampleSheet_Bulk_RNA.csv', sep = ',', header = T)
colnames(NF_META)[3] <- 'Group'

# Bulk RNA-seq TOX table (GSE291316_Table_Norm_Stats_Batch4.txt)
TOX_Table  <- read.table('C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Saves/Bulk_RNA-seq/Tox/Table_Results.txt')
TOX_META   <- read.table('C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Data/Bulk_RNA-seq/Tox/SampleSheet_Bulk_RNA.csv', sep = ',', header = T)
colnames(TOX_META)[3] <- 'Group'

# scRNA-seq saved object for all datasets
scRNAseq   <- readRDS('C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Saves/scRNA-seq/KO/5_Pseudotime_ALL.rds')


# scRNA-seq LEAP Candidates (Available in GitHub)
LEAP_DOWN <- read.table(
  paste0(PATH,'/Saves/scRNA-seq/WT/LEAP/LEAP_Candidates_DOWN.txt'))[,1]
LEAP_UP   <- read.table(
  paste0(PATH,'/Saves/scRNA-seq/WT/LEAP/LEAP_Candidates_UP.txt'))[,1]

# Bulk RNA-seq Candidates (Available in GitHub)
Bulk_DOWN <- read.table(
  paste0(PATH,'/Saves/Bulk_RNA-seq/Nfil3/Bulk_Candidates_DOWN.txt'))[,1]
Bulk_UP   <- read.table(
  paste0(PATH,'/Saves/Bulk_RNA-seq/Nfil3/Bulk_Candidates_UP.txt'))[,1]


# DNase-seq transiently open regions (Available in GitHub)
Transient_regions <- read.table(paste0(
  PATH,'/Data/DNase-seq/Transiently_Open_Clustered_K3_Annotated.txt'), 
  sep = '\t', header = T)

# CUT&RUN peak annotation results (GSE291321_Table_HOMER_NFIL3_peaks.txt)
HA_Annotated <- read.table(paste0(
  PATH,'/Data/CUT&RUN/HOMER/ALP_HA_filtered_merged_peaks_annotated_input.txt'), 
  sep='\t', header = T)
colnames(HA_Annotated)[1] <- 'Peak_ID'

# Manually re-annotated peaks (Available in GitHub)
peaks_manual <- read.table(paste0(
  PATH,'/Data/CUT&RUN/HOMER/Peak_Manual_Annotation.csv'), sep = ',', header = T)




#===============================================================================
## CANDIDATES in vitro / ex vivo -----------------------------------------------
#===============================================================================

# Get shared candidates
Shared_DOWN <- LEAP_DOWN[LEAP_DOWN %in% Bulk_DOWN]
Shared_UP   <- LEAP_UP[LEAP_UP %in% Bulk_UP]
Shared_ALL  <- c(Shared_DOWN, Shared_UP)

# Draw Venn Diagrams - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
display_venn(x = list(LEAP_DOWN, Bulk_DOWN), 
             category.names = c('ex vivo','in vitro'))
display_venn(x = list(LEAP_UP, Bulk_UP), 
             category.names = c('ex vivo','in vitro'))

# Draw heatmap of shared candidates - - - - - - - - - - - - - - - - - - - - - - 
DrawHeatmap(NF_Table[,c(1:28)], NF_META, c(Shared_DOWN, Shared_UP), 
            by.group = T, cluster_rows = T,
            groups =c('GFP', 'NFIL3','ALP', 'sEILP', 'cEILP', 'ILCP'),
            gaps_col = 2, treeheight_row = 0, fontsize_row = 5)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save shared candidates lists
write.table(Shared_DOWN, paste0(PATH_SAVE, '/Shared_Candidates_DOWN.txt'), 
            quote = F, col.names = F, row.names = F)
write.table(Shared_UP, paste0(PATH_SAVE, '/Shared_Candidates_UP.txt'), 
            quote = F, col.names = F, row.names = F)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



################################################################################
#            READY TO USE  CANDIDATES LISTS FOR METASCAPE ANNOTATION           #
#            https://metascape.org/gp/index.html#/main/step1                   #
################################################################################
#
# After annotation, full results in one .zip file are downloaded and stored.
# Following function are working directly from this unique file to generate 
# plots from resulting annotation.
# 
################################################################################
#                                                                              #
################################################################################
source(paste0('https://raw.githubusercontent.com/JosephLeger/Metascape_Toolkit',
              '/refs/heads/main/Custom_Functions.R'))

Meta_down <- meta_Read(paste0(
  PATH, '/Data/Metascape/Metascape_30_Candidates_DOWN.zip'))
Meta_up   <- meta_Read(paste0(
  PATH, '/Data/Metascape/Metascape_72_Candidates_UP.zip'))

meta_DrawEnrichment(Meta_down, stat = 'P', title = '30 Candidates Down')
meta_DrawEnrichment(Meta_up, stat = 'P', title = '72 Candidates Up', 
                    stat.cutoff = 3)



#===============================================================================
## MISREGULATED GENES IN TOX-KO ------------------------------------------------
#===============================================================================

## SINGLE-CELL MISREGULATED GENES ----------------------------------------------

# Subset scRNA-seq data to keep only Il7rLT+ cells (WT, NFIL3-KO or TOX-KO)
ALL_Il7r_LT <- subset(scRNAseq, orig.ident %in% c('overLIP', 'NF-KO', 'TOX-KO'))
table_ALL   <- OrderMatrix(ALL_Il7r_LT, 'Slingshot_2' , min.cells = 50,  
                            slot = 'scaledata')

# Annotate sample genotype
Sample <- c()
Ann <- data.frame(Cell.ID = ALL_Il7r_LT@assays[["RNA"]]@counts@Dimnames[[2]], 
                  Sample = ALL_Il7r_LT@meta.data[["orig.ident"]])
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

# Compare expression patterns across samples - - - - - - - - - - - - - - - - - -
stat_list <- list()
for(gene in Shared_ALL){
  stat_list[[gene]] <- CompareExpr(x = table_ALL, bin.number = 11, 
                                   feature = gene, by.order = F, scale = T, 
                                   col = ColorBlind[c(6,3,8)], ylim = NULL, 
                                   min.cells = 5, stat = 't')
  ggsave(paste0(PATH_FIG, '/Curves/', gene, '.png'))
}

# Identify significant genes between TOX-KO and WT
TOX_KO_sig <- c()
for(i in 1:length(stat_list)){
  x <- data.frame(stat_list[[i]])[c(6:11),]
  if(min(x$WT_vs_TOX.KO, na.rm = T) < 0.05){
    TOX_KO_sig <- c(TOX_KO_sig, names(stat_list)[i])
  }
}
TOX_KO_sig



## VALIDATION IN TOX INDUCED RNA-SEQ -------------------------------------------

# Looking for genes differentially expressed after TOX ectopic expression in
# WT or NFIL3-KO ALP

TOX_in_WT <- c()
TOX_in_KO <- c()

Stats_WT <- 'TOX_WT_vs_GFP_WT_pval'
Stats_KO <- 'TOX_KO_vs_GFP_KO_pval'

for(gene in Shared_ALL){
  is_not_sig <- T
  # Test significance in WT cells
  if(gene %in% TOX_Table$Symbol && 
     !is.na(TOX_Table[,Stats_WT][TOX_Table$Symbol %in% gene]) &&
     TOX_Table[,Stats_WT][TOX_Table$Symbol %in% gene] < 0.05){
    is_not_sig <- F
    TOX_in_WT  <- c(TOX_in_WT, gene)
    p          <- DrawDotplot(TOX_Table, TOX_META, gene, 
                              c('GFP_WT', 'TOX_WT', 'GFP_KO', 'TOX_KO'), 
                              repeated = NA, colstat = Stats_WT, show.fc = T, 
                              main = paste(gene, 'in WT'))
    ggsave(paste0(PATH_FIG, '/Dotplots/WT/', gene, '.png'), p)
  }
  # Test significance in NFIL3-KO cells
  if(gene %in% TOX_Table$Symbol && 
     !is.na(TOX_Table[,Stats_KO][TOX_Table$Symbol %in% gene]) &&
     TOX_Table[,Stats_KO][TOX_Table$Symbol %in% gene] < 0.05){
    is_not_sig <- F
    TOX_in_KO  <- c(TOX_in_KO, gene)
    p          <- DrawDotplot(TOX_Table, TOX_META, gene, 
                              c('GFP_WT', 'TOX_WT', 'GFP_KO', 'TOX_KO'), 
                              repeated = NA, colstat = Stats_KO, show.fc = T, 
                              main = paste(gene, 'in NFIL3-KO'))
    ggsave(paste0(PATH_FIG, '/Dotplots/KO/', gene, '.png'), p)
  }
  # Otherwise draw plot in global Dotplots directory
  if(is_not_sig && gene %in% TOX_Table$Symbol){
    p          <- DrawDotplot(TOX_Table, TOX_META, gene, 
                              c('GFP_WT', 'TOX_WT', 'GFP_KO', 'TOX_KO'), 
                              repeated = NA, colstat = NA, show.fc = T, 
                              main = paste(gene))
    ggsave(paste0(PATH_FIG, '/Dotplots/', gene, '.png'), p)
  }
}



#===============================================================================
## NFIL3 CUT&RUN ANALYSIS ------------------------------------------------------
#===============================================================================

## REFINED CANDIDATES ----------------------------------------------------------

# Add manual annotation for peaks in known regulatory regions
peak_added <- HA_Annotated[HA_Annotated$Peak_ID %in% peaks_manual$Peak_ID,]
peak_added <- merge(peak_added, peaks_manual, 'Peak_ID')
peak_added$Gene.Name <- peak_added$Gene
peak_added$Gene.Type <- 'manually_annotated'
peak_added$Gene      <- NULL

HA_Rescued <- HA_Annotated[HA_Annotated$Peak_ID %!in% peak_added$Peak_ID,]
HA_Rescued <- rbind(HA_Rescued, peak_added)

# Get intersection between CUT&RUN peaks and ex vivo/in vitro shared candidates
HA_cand_enriched_peaks <- HA_Rescued[HA_Rescued$Gene.Name %in% Shared_ALL,] 
HA_cand_enriched_bed   <- HA_cand_enriched_peaks[,c(2,3,4,1)]
HA_cand_enriched       <- unique(HA_cand_enriched_peaks$Gene.Name)

# Draw heatmap of refined candidates
DrawHeatmap(NF_Table[,c(1:28)], NF_META, HA_cand_enriched, by.group = T, 
            groups =c('GFP', 'NFIL3','ALP', 'sEILP', 'cEILP', 'ILCP'),
            gaps_col = 2)


## UP/DOWN REGULATED CANDIDATES ANALYSIS ---------------------------------------

HA_cand_enriched_down_peaks <- HA_cand_enriched_peaks[HA_cand_enriched_peaks$Gene.Name %in% Shared_DOWN,]
HA_cand_enriched_up_peaks   <- HA_cand_enriched_peaks[HA_cand_enriched_peaks$Gene.Name %in% Shared_UP,]
HA_cand_enriched_down_bed   <- HA_cand_enriched_down_peaks[,c(2,3,4,1)]
HA_cand_enriched_up_bed     <- HA_cand_enriched_up_peaks[,c(2,3,4,1)]

HA_cand_enriched_down <- unique(HA_cand_enriched_down_peaks$Gene.Name)
HA_cand_enriched_up   <- unique(HA_cand_enriched_up_peaks$Gene.Name)

write.table(HA_cand_enriched_down_bed, 'C:/Users/E15639P/Desktop/Review/Saves/CUT&RUN/HA_cand_enriched_down_peaks_input.bed',
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(HA_cand_enriched_up_bed, 'C:/Users/E15639P/Desktop/Review/Saves/CUT&RUN/HA_cand_enriched_up_peaks_input.bed',
            sep = '\t', quote = F, row.names = F, col.names = F)

# Draw PieCharts of regions distribution
Regions_PieChart(HA_cand_enriched_down_peaks)
Regions_PieChart(HA_cand_enriched_up_peaks)


## SUPPLEMENTARY FIGURES -------------------------------------------------------

# Draw PieChart of overall 17572 peaks regions distribution 
Regions_PieChart(HA_Annotated)


# Draw BoxPlots of gene expressions 
# All genes annotated from called peaks
box_gfp   <- NF_Table$mean_GFP[
  NF_Table$Symbol %in% unique(HA_Annotated$Gene.Name)]
box_nf    <- NF_Table$mean_NFIL3[
  NF_Table$Symbol %in% unique(HA_Annotated$Gene.Name)]

stats     <- t.test(box_gfp, box_nf)
box_table <- data.frame(Expression = log(c(box_gfp, box_nf)+1),
                        Group = c(rep(c('GFP', 'NFIL3'), each = length(box_gfp))))

ggplot(box_table, aes(x=Group, y=Expression, fill=Group)) +
  geom_boxplot(outlier.size = 0) +
  scale_fill_manual(values=ColorBlind[c(16,3)]) +
  geom_jitter(color="black", size=0.5, alpha=0.9) +
  theme_classic() +
  theme(legend.position="none") +
  ggtitle(paste('p-value =', stats$p.value )) +
  xlab("")

# Focus on DNase-seq transiently open regions by cluster
for(i in 1:3){
  box_gfp <- NF_Table$mean_GFP[
    NF_Table$Symbol %in% unique(
      Transient_regions$Symbol[Transient_regions$Cluster %in% paste0('cluster_', i)])]
  box_nf  <- NF_Table$mean_NFIL3[
    NF_Table$Symbol %in% unique(
      Transient_regions$Symbol[Transient_regions$Cluster %in% paste0('cluster_', i)])]
  
  stats     <- t.test(box_gfp, box_nf)
  box_table <- data.frame(Expression = log(c(box_gfp, box_nf)+1),
                          Group = c(rep(c('GFP', 'NFIL3'), each = length(box_gfp))))
  
  p <- ggplot(box_table, aes(x=Group, y=Expression, fill=Group)) + 
    geom_boxplot(outlier.size = 0) +
    scale_fill_manual(values=ColorBlind[c(16,3)]) +
    geom_jitter(color="black", size=0.5, alpha=0.9) +
    theme_classic() +
    theme(legend.position="none") +
    ggtitle(paste('Cluster',i,'| p-value =', stats$p.value )) +
    xlab("")
  
  print(p)
}





