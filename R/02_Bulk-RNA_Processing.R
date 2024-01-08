#!/usr/bin/env Rscript


#===============================================================================
# DESCRIPTION ------------------------------------------------------------------
#===============================================================================

# Based on CS023454 and CS021176 experiments and corresponding formated Tables
# previously generated in '/Doctorat/Bulk_RNA-seq/Data_Format.'
# Here we filter Tables to keep every sample (ex-vivo or tranduced in vitro) 
# that could be meaningful for ILC development analysis :
# ALP, CLP, sEILP, PLZFneg, cEILP, PLZFpos, Tulip, EILPWT, ALP+GFP, CLP+GFP and
# ALP+NFIL3
# To ask biological questions, further analysis require more precise filtering



#===============================================================================
# SETUP ------------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

PATH         <- 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC'
DATA_DIR     <- 'C:/Users/E15639P/Data/Bulk/NFIL3_dev_ILC/genes.results'
SAMPLE_SHEET <- 'C:/Users/E15639P/Data/Bulk/NFIL3_DEV_ILC/SampleSheet_Bulk_RNA.csv'
setwd(PATH)

source('C:/Users/E15639P/Doctorat/Bulk_RNA-seq/Custom_functions.R')
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(org.Mm.eg.db))

# Create required directories
dir.create(file.path(paste0(PATH, '/Figures'), 'Bulk-RNA'))
dir.create(file.path(paste0(PATH, '/Saves'), 'Bulk-RNA'))

PATH_FIG  <- paste0(PATH, '/Figures/Bulk-RNA')
PATH_SAVE <- paste0(PATH, '/Saves/Bulk-RNA')



#===============================================================================
# READING FILES ----------------------------------------------------------------
#===============================================================================

# Read metadata sample sheet
METADATA        <- read.table(SAMPLE_SHEET, header = T, sep = ',')

# Import files
TXI             <- tximport(paste(DATA_DIR, METADATA$FileName, sep = '/'), 
                            type = 'rsem')
Table           <- as.data.frame(TXI$counts)
colnames(Table) <- METADATA$SampleName

# Save aggregated Table_Raw available in GEO
write.table(Table, 
            paste0(PATH_SAVE, '/Table_Raw.txt'), quote = F)



#===============================================================================
# DATA DISTRIBUTION BEFORE NORMALIZATION  --------------------------------------
#===============================================================================

# Make heatmaps by celltype
for(n in names(table(METADATA$CellType))){
  heatdata <- Table[colnames(Table) %in% 
                      METADATA$SampleName[METADATA$CellType %in% n]]

  # Eliminate low expressed genes for heatmap
  heatdata$count <- apply(heatdata, 1, sum)
  heatdata       <- subset(heatdata, count > 10)
  heatdata$count <- NULL

  # Plot and save heatmap
  pheatmap(log2(heatdata+1), 
           show_rownames = FALSE, 
           treeheight_row = 50, 
           treeheight_col = 50,
           cluster_cols = FALSE,
           main = n) 
  dev.print(png, file=paste0(PATH_FIG, '/Heatmap_', n,'.png'), 
            width=9, height=9, units='in', res=100)
}


# Make complete heatmap
pheatdata       <- Table
pheatdata$count <- apply(pheatdata, 1, sum)
pheatdata       <- subset(pheatdata, count > 10)
pheatdata$count <- NULL

# Plot and save heatmap
pheatmap(log2(pheatdata+1), 
         show_rownames = FALSE, 
         treeheight_row = 50,
         treeheight_col = 50, 
         luster_cols = FALSE)
dev.print(png, file=paste0(PATH_FIG, '/Complete_Heatmap.png'), 
          width=9, height=9, units='in', res=100)

# Make dendrogram of raw data
dendodata       <- t(Table)
dist            <- dist(dendodata[ ,c(1:ncol(dendodata))], 
                        diag=TRUE, 
                        method = 'euclidian')
hc              <- hclust(dist,
                          method = 'complete')
dendro_raw      <- plot(hc, 
                        main = 'Raw Samples', 
                        xlab = 'Samples', 
                        sub = '')
dev.print(png, file=paste0(PATH_FIG, '/Dendro_Raw.png'), 
          width=9, height=9, units='in', res=100)



#===============================================================================
# QUANTILE NORMALIZATION & BATCH EFFECT REMOVING -------------------------------
#===============================================================================

# Eliminate zero expressed genes
cutoff         <- 1
eliminate      <- which(apply(cpm(Table), 1, max) < cutoff)
Table_filtered <- Table[-eliminate,]

# Save aggregated Table after non-expressed genes removal
write.table(Table_filtered, 
            paste0(PATH_SAVE, '/Table_Filtered.txt'), quote = F)

# Plot distribution before quantile-normalization
Table_filtered %>% 
  gather(Sample, Count) %>% 
  ggplot(aes(Sample, Count)) + 
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 60)
dev.print(png, file=paste0(PATH_FIG, '/Distrib_Raw.png'), 
          width=9, height=9, units='in', res=100)

# Use custom function for quantile naormalization
Table_norm  <- Quantile_Normalization(Table_filtered)

# Plot distribution after quantile normalization
Table_norm %>% 
  gather(Sample, Count) %>% 
  ggplot(aes(Sample, Count)) + 
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 60)
dev.print(png, file=paste0(PATH_FIG, '/Distrib_QuantileNorm.png'), 
          width=9, height=9, units='in', res=100)


# Batch effect removal
Table_batched <- ComBat_seq(Table_norm, batch = METADATA$Batch)

# Plot dendrogram after normalization and batch effect removal
dendodata     <- t(Table_batched)
dist          <- dist(dendodata[,c(1:ncol(dendodata))],
                      diag=TRUE, 
                      method = 'euclidian')
hc            <- hclust(dist, 
                        method = 'complete')
plot(hc, main = '', xlab = 'Samples', sub = '')
dev.print(png, file=paste0(PATH_FIG, '/Dendro_Norm.png'), 
          width=9, height=9, units='in', res=100)

# Save aggregated Table_Norm available in GEO
write.table(Table_batched, 
            paste0(PATH_SAVE, '/Table_Norm.txt'), quote = F)



#===============================================================================
# GENE SYMBOL ANNOTATION -------------------------------------------------------
#===============================================================================

ID     <- rownames(Table_batched)
mapped <- mapIds(org.Mm.eg.db, 
                 keys = ID,
                 keytype = 'ENSEMBL', 
                 column = 'SYMBOL')
Symbol <- c()
for(j in 1:length(ID)){
  if(!is.na(mapped[j])){
    Symbol <- c(Symbol, mapped[j])
  }else{
    Symbol <- c(Symbol, ID[j])
  }
}
Table_annotated <- cbind(Table_batched, Symbol = Symbol)
Table_annotated <- Table_annotated[,c(ncol(Table_annotated), 
                                      1:(ncol(Table_annotated)-1))]

# Save table with annotated gene symbol
write.table(Table_annotated, 
            paste0(PATH_SAVE, '/Table_Annotated.txt'), quote = F)












