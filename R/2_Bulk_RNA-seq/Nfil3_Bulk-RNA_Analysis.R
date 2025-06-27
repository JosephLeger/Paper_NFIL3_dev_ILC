#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
# Script used to perform differential gene expression analysis of Bulk RNA-seq
#
# Load files using TXImport
# Check sample distribution to identify potential outliers
# Apply quantile normalization
# Remove batch effect using CombatSeq
# Annotate gene symbols using org.Mm.eg.db package
# Perform statistical DEG analysis using Limma



#===============================================================================
## SETUP -----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PROJECT INFO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

PATH         <- 'C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC'
DATA_DIR     <- 'C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Data/Bulk_RNA-seq/Nfil3/genes.results'
SAMPLE_SHEET <- 'C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Data/Bulk_RNA-seq/Nfil3/SampleSheet_Bulk_RNA_Nfil3.csv'
COMP_TO_MAKE <- 'C:/Users/E15639P/Desktop/GitHub_NF_dev_ILC/Data/Bulk_RNA-seq/Nfil3/Comparisons_to_make_Nfil3.csv' 

################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

setwd(PATH)

# Load Packages and custom functions
source('C:/Users/E15639P/Doctorat/Custom_Functions.R')
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(EnhancedVolcano))

# Create required directories
dir.create(file.path(PATH, '/Figures'))
dir.create(file.path(PATH, '/Saves'))

dir.create(file.path(paste0(PATH, '/Figures'), 'Bulk_RNA-seq'))
dir.create(file.path(paste0(PATH, '/Figures/Bulk_RNA-seq'), 'Nfil3'))
dir.create(file.path(paste0(PATH, '/Saves'), 'Bulk_RNA-seq'))
dir.create(file.path(paste0(PATH, '/Saves/Bulk_RNA-seq'), 'Nfil3'))

# Set up path for figures and saves
PATH_FIG  <- paste0(PATH, '/Figures/Bulk_RNA-seq/Nfil3')
PATH_SAVE <- paste0(PATH, '/Saves/Bulk_RNA-seq/Nfil3')



#===============================================================================
## READING FILES ---------------------------------------------------------------
#===============================================================================

# Read metadata sample sheet and comparisons to make
METADATA        <- read.table(SAMPLE_SHEET, header = T, sep = ',')
COMPARISONS     <- read.csv(COMP_TO_MAKE)

# Import files
TXI             <- tximport(paste(DATA_DIR, METADATA$FileName, sep = '/'), 
                            type = 'rsem')
Table           <- as.data.frame(TXI$counts)
colnames(Table) <- METADATA$SampleName

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save Table_Raw available in GEO
write.table(Table, paste0(PATH_SAVE, '/Table_Raw.txt'), quote = F)
Table <- read.table(paste0(PATH_SAVE, '/Table_Raw.txt'), header = T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## SAMPLE DISTRIBUTION BEFORE NORMALIZATION  -----------------------------------
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
## QUANTILE NORMALIZATION & BATCH EFFECT REMOVAL -------------------------------
#===============================================================================

# Eliminate zero expressed genes
cutoff         <- 1
eliminate      <- which(apply(cpm(Table), 1, max) < cutoff)
Table_filtered <- Table[-eliminate,]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save Table after non-expressed genes removal
write.table(Table_filtered, paste0(PATH_SAVE, '/Table_Filtered.txt'), quote = F)
Table_filtered <- read.table(paste0(PATH_SAVE, '/Table_Filtered.txt'), header=T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


# Plot distribution before quantile-normalization
Table_filtered %>% 
  gather(Sample, Count) %>% 
  ggplot(aes(Sample, Count)) + 
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 60)
dev.print(png, file=paste0(PATH_FIG, '/Distrib_Raw.png'), 
          width=9, height=9, units='in', res=100)

# Use custom function for quantile normalization
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

# Plot PCA of sample distribution
pca_res <- prcomp(t(Table_batched))
pca_sum <- summary(pca_res)

label_x = paste0('PC1 ', '(', round(pca_sum[["importance"]][2]*100), '%)')
label_y = paste0('PC2 ', '(', round(pca_sum[["importance"]][5]*100), '%)')

pca_t   <- data.frame(
  pca_res$x, 
  Sample=factor(METADATA$CellType, 
                levels = c('ALP', 'sEILP', 'cEILP', 'ILCP', 'GFP', 'NFIL3')))

ggplot(pca_t, aes(x=PC1, y=PC2, color=Sample)) +
  theme_classic() +
  geom_point(size = 5, shape = 20) +
  scale_color_manual(
    values=c(ColorBlind[c(10,9,8)], 'darkred', ColorBlind[c(3,5)])) +
  xlab(label_x) + ylab(label_y)
dev.print(png, file=paste0(PATH_FIG, '/PCA_plot.png'), 
          width=9, height=9, units='in', res=100)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save Table_Norm
write.table(Table_batched, paste0(PATH_SAVE, '/Table_Norm.txt'), quote = F)
Table_batched <- read.table(paste0(PATH_SAVE, '/Table_Norm.txt'), header = T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## GENE SYMBOL ANNOTATION ------------------------------------------------------
#===============================================================================

ID     <- rownames(Table_batched)
mapped <- mapIds(org.Mm.eg.db, 
                 keys = ID,
                 keytype = 'ENSEMBL', 
                 column = 'SYMBOL')
mapped <- data.frame(ID = names(mapped), Symbol = mapped)

# Define function to fill in faster than iteration
fill_if_NA <- function(line){
  if(is.na(line[2])){
    return(line[1])
  }else{
    return(line[2])
  }
}
Symbol           <-  apply(mapped, 1, fill_if_NA)
Table_annotated  <- cbind(Symbol = Symbol, Table_batched)
Table_annotated  <- Table_annotated[order(rownames(Table_annotated)),]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save table with annotated gene symbol
write.table(Table_annotated, paste0(PATH_SAVE, '/Table_Annotated.txt'), quote=F)
Table_annotated <- read.table(paste0(PATH_SAVE,'/Table_Annotated.txt'),header=T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 



#===============================================================================
## FITTING MODEL ---------------------------------------------------------------
#===============================================================================

# Create DGEList object
DGE_0     <- DGEList(Table_filtered)

## PREPROCESSGING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DGE_0     <- calcNormFactors(DGE_0, method = 'TMM')

# Remove low expressed genes for the following analysis
cutoff    <- 10
eliminate <- which(apply(cpm(DGE_0), 1, max) < cutoff)
DGE       <- DGE_0[-eliminate,] 

# Show sample distribution before VOOM
samples   <- METADATA$SampleName
group     <- METADATA$CellType
batch     <- as.character(METADATA$Batch)

colors <- c()
for(i in group){
  for(g in 1:length(names(table(group)))){
    if(i == names(table(group))[g]){
      colors <- c(colors, g)
    }
  }
}
plotMDS(DGE, col = colors)
dev.print(png, file=paste0(PATH_FIG, '/plotMDS_before.png'), 
          width=9, height=9, units='in', res=100)


## VOOM TRANSFORMATION AND CALCULATION OF VARIANCE WEIGHT - - - - - - - - - - -

model <- model.matrix(~0 + group + batch)

y_0   <- voom(DGE_0, model, plot = TRUE)
y     <- voom(DGE, model, plot = TRUE)
plotMDS(y, col = colors)
dev.print(png, file=paste0(PATH_FIG, '/plotMDS_after.png'), 
          width=9, height=9, units='in', res=100)
fit   <- lmFit(y, model)



#===============================================================================
## MAKE COMPARISONS ------------------------------------------------------------
#===============================================================================

Results_list <- list()

for(i in 1:nrow(COMPARISONS)){
  
  ## SELECT COMPARISON ---------------------------------------------------------
  # Set current comparison members
  control <- COMPARISONS$Control[i] 
  tested  <- COMPARISONS$Tested[i]
  title   <- paste(tested, "vs", control, sep = "_")
  pattern <- paste0("group", tested, " - ", "group", control)
  print(pattern)
  
  # Run comparison
  contr   <- makeContrasts(pattern,  levels = colnames(coef(fit)))
  tmp     <- contrasts.fit(fit, contr)
  tmp     <- eBayes(tmp)
  
  result  <- topTable(tmp, sort.by = "P", n = Inf)
  
  # Adding gene symbol
  symbol <- mapIds(org.Mm.eg.db, keys = rownames(result), keytype = "ENSEMBL", 
                   column = "SYMBOL")
  symbol.completed <- c()
  for(j in 1:nrow(result)){
    if(!is.na(symbol[j])){
      symbol.completed <- c(symbol.completed, symbol[j])
    }else{
      symbol.completed <- c(symbol.completed, rownames(result)[j])
    }
  }
  result$Symbol <- symbol.completed
  result <- result[,c("Symbol", 
                      colnames(result)[colnames(result) %!in% "Symbol"])]
  result <- result[order(rownames(result)),]
  Results_list[[title]] <- result

  ## SAVE RESULT TABLE ---------------------------------------------------------
  write.table(result, paste0(PATH_SAVE,'/LimmaStats_',title,'.txt'), quote = F)
  
  
  ## PLOT MA -------------------------------------------------------------------
  ggplot(as.data.frame(result), 
         aes(x = as.numeric(AveExpr), y = as.numeric(logFC))) +
    geom_point() +
    ggtitle(title) +
    labs(x = 'Log10(baseMean)', y = 'Log2 FoldChange') +
    scale_color_manual('Significant', values = ColorBlind[c(1,6)])
  ggsave(paste0(PATH_FIG, '/', title, '_MAplot.png'),  
         width = 3000, height = 2500, units = "px")
  
  
  ## VOLCANO PLOT WITH PVALUE --------------------------------------------------
  keyvals <- ifelse(
    result$logFC < -1 & result$P.Value < 0.05, 'royalblue',
    ifelse(result$P.Value < 0.05 & result$logFC > 1, 'red',
           'grey'))
  result$keyvals <- keyvals
  
  
  keyvals[is.na(keyvals)]                <- 'grey'
  names(keyvals)[keyvals == 'red']       <- 'Up'
  names(keyvals)[keyvals == 'grey']      <- 'NS'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down'
  # Draw Volcano plot
  EnhancedVolcano(result , lab = result$Symbol, x = 'logFC', y = 'P.Value', 
                  title = title, subtitle = '',
                  legendPosition = 'right',
                  selectLab = rownames(result$Symbol)[
                    which(names(keyvals) %in% c('Up', 'Down'))],
                  pCutoff = 0.05, FCcutoff = 1, 
                  pointSize = 1.2, labSize = 3,
                  colCustom = keyvals, colAlpha = 1,
                  ylab = bquote(~-Log[10] ~ italic(P.Value)))
  ggsave(paste0(PATH_FIG, '/', title, '_Volcano_pval.png'), 
         width = 3000, height = 2500, units = "px")
  
  ## SAVE GENE LIST FOR PVALUE -------------------------------------------------
  UP   <- unique(result$Symbol[result$keyvals == "red"])
  DOWN <- unique(result$Symbol[result$keyvals == "royalblue"])

  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
  write.table(UP, paste0(PATH_SAVE, '/', title, '_UP.txt'), 
              row.names = F, col.names = F, quote = F)
  write.table(DOWN, paste0(PATH_SAVE, '/', title, '_DOWN.txt'),
              row.names = F, col.names = F, quote = F)
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
}



#===============================================================================
## FINAL TABLE RESULT ----------------------------------------------------------
#===============================================================================
# Generate full table annotated with expression and stats
FULL <- read.table(paste0(PATH_SAVE, '/Table_Annotated.txt'), header = T)

# Read results files 
Results_list <- list()
files        <- list.files(PATH_SAVE)
files        <- files[grepl('LimmaStats', files)]
for(f in files){
  name                 <- str_remove(str_remove(f, 'LimmaStats_'), '.txt')
  Results_list[[name]] <- read.table(paste0(PATH_SAVE, '/', f))
}

# Calculate mean expression for each celltype
for(n in names(table(METADATA$CellType))){
  subx <- FULL[,METADATA$SampleName[METADATA$CellType %in% n]]
  mean <- rowMeans(subx)
  FULL <- cbind(FULL, newcol = mean)
  colnames(FULL)[ncol(FULL)] <- paste0('mean_', n)
}

# Separate genes with and without statistical results
plus <- FULL[rownames(FULL) %!in% rownames(Results_list[[1]]),]
FULL <- FULL[rownames(FULL) %in% rownames(Results_list[[1]]),]

# Add statistics 
for(i in 1:length(Results_list)){
  Stats <- Results_list[[i]]
  Stats <- Stats[,c('logFC', 'P.Value', 'adj.P.Val')]
  colnames(Stats) <- paste0(names(Results_list[i]), 
                            c('_logFC', '_pval', '_padj'))
  FULL <- cbind(FULL, Stats)
}

# Add genes with no statistical results
for(i in 1:(ncol(FULL)-ncol(plus))){
  plus <- cbind(plus, NA)
  colnames(plus)[ncol(plus)] <- colnames(FULL)[ncol(plus)]
}
FULL <- rbind(FULL, plus)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save final normalized table with stats available in GEO
write.table(FULL, paste0(PATH_SAVE, '/Table_Results.txt'), quote=F, sep='\t')
FULL <- read.table(paste0(PATH_SAVE, '/Table_Results.txt'), header = T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## ANALYZE AND PLOT RESULTS ----------------------------------------------------
#===============================================================================

Stats  <- read.table(paste0(PATH_SAVE, '/LimmaStats_sEILP_vs_ALP.txt'))
DEG_UP <- Stats$Symbol[Stats$P.Value < 0.05 & Stats$logFC > 1]
TF_fam <- readRDS(paste0(PATH, '/Saves/DNAse-seq/TF_fam.RDS'))


## BZIP FACTORS ALP VS EILP - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Significant genes considering only Subfamilies
TF_fam$BZIP[TF_fam$BZIP %in% DEG_UP]
TF_fam$PAR[TF_fam$PAR %in% DEG_UP]

# Significant genes while extending analysis to all bZIP factors
TF_fam$BZIP_large[TF_fam$BZIP_large %in% DEG_UP]


## TRANSCRIPTION FACTORS INVOLVED IN ILC DEV - - - - - - - - - - - - - - - - - -

celltype_order <- c("ALP", "GFP", "sEILP", "cEILP","NFIL3", "ILCP")
gene_list      <- c("Tcf4","Lmo2","Mef2c","Irf8","Spi1","Nfil3","Tcf7","Zbtb16",
                    "Tox","Gata3","Batf3","Maf","Id2")

x           <- FULL[FULL$Symbol %in% gene_list, 
                    c('Symbol', paste0('mean_', celltype_order))]
colnames(x) <- c('Symbol', celltype_order)
rownames(x) <- x$Symbol
x           <- x[gene_list,]

# Drawing heatmap of factors expression along ILC dev / NFIL3 induction
pheatmap(x[,2:7], show_rownames = T, treeheight_row = 50, treeheight_col = 50,
         cluster_rows = F, cluster_cols = F, scale = 'row',
         color = colorRampPalette(c("blue","white","red"))(50), main = "", 
         silent=F)
dev.print(png, file=paste0(PATH_FIG, '/TF_dev_ILC.png'), 
          width=9, height=9, units='in', res=100)


## NFIL3 CANDIDATE TARGET GENES  - - - - - - - - - - - - - - - - - - - - - - - -
# Cross list of genes
Genelist <- list()

for(n in names(Results_list)){
  Genelist[[paste0(n, '_UP')]]   <- read.table(paste0(PATH_SAVE, '/', n, 
                                                      '_UP.txt'))[,1]
  Genelist[[paste0(n, '_DOWN')]] <- read.table(paste0(PATH_SAVE, '/', n, 
                                                      '_DOWN.txt'))[,1]
}

Candidates.UP   <- Genelist[['NFIL3_vs_GFP_UP']][
  Genelist[['NFIL3_vs_GFP_UP']] %in% Genelist[['NFIL3_vs_ALP_UP']]]
Candidates.DOWN <- Genelist[['NFIL3_vs_GFP_DOWN']][
  Genelist[['NFIL3_vs_GFP_DOWN']] %in% Genelist[['NFIL3_vs_ALP_DOWN']]]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
write.table(Candidates.UP, paste0(PATH_SAVE, '/Candidates_UP.txt'), 
            row.names = F, col.names = F, quote = F)
write.table(Candidates.DOWN, paste0(PATH_SAVE, '/Candidates_DOWN.txt'), 
            row.names = F, col.names = F, quote = F)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =






