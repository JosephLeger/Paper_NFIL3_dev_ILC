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

source('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Custom_functions.R')
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(pheatmap))

PATH_FIG  <- paste0(PATH, '/Figures/Bulk-RNA')
PATH_SAVE <- paste0(PATH, '/Saves/Bulk-RNA')


#===============================================================================
# READING FILES ----------------------------------------------------------------
#===============================================================================

METADATA <- read.table(SAMPLE_SHEET, header = T, sep = ',')
Table    <- read.table(paste0(PATH_SAVE, '/Table_Filtered.txt'))



#===============================================================================
# FITTING MODEL ----------------------------------------------------------------
#===============================================================================

# Create DGEList object
DGE_0     <- DGEList(Table)

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
# MAKE COMPARISONS -------------------------------------------------------------
#===============================================================================

# Run comparison
contrast <- makeContrasts('groupsEILP - groupALP', 
                          levels = colnames(coef(fit)))
tmp      <- contrasts.fit(fit, contrast)
tmp      <- eBayes(tmp)
result   <- topTable(tmp, sort.by = 'P', n = Inf)

# Adding gene symbol
symbol   <- mapIds(org.Mm.eg.db, 
                   keys = rownames(result), 
                   keytype = 'ENSEMBL', 
                   column = 'SYMBOL')
symbol.completed <- c()
for(j in 1:nrow(result)){
  if(!is.na(symbol[j])){
    symbol.completed <- c(symbol.completed, symbol[j])
  }else{
    symbol.completed <- c(symbol.completed, rownames(result)[j])
  }
}
result$Symbol <- symbol.completed
result <- result[,c('Symbol', colnames(result)[
  colnames(result) %!in% 'Symbol'])]
result <- result[order(rownames(result)),]

# Save stat results
write.table(result, 
            paste0(PATH_SAVE, '/LimmaStats_sEILP_vs_ALP.txt'))

# Extract significant gene lists
result$Sig <- ifelse(
  result$logFC < -1 & result$P.Value < 0.05, 'DOWN',
  ifelse(result$P.Value < 0.05 & result$logFC > 1, 'UP','NS'))

genes_down <- result$Symbol[result$Sig %in% 'DOWN']
genes_up   <- result$Symbol[result$Sig %in% 'UP']

# Save significant gene lists
write.table(genes_down, paste0(PATH_SAVE, '/sEILP_vs_ALP_DOWN.txt'), 
            col.names = F, 
            row.names = F,
            quote = F)
write.table(genes_up, paste0(PATH_SAVE, '/sEILP_vs_ALP_UP.txt'),
            col.names = F, 
            row.names = F,
            quote = F)



#===============================================================================
# FINAL TABLE RESULT -----------------------------------------------------------
#===============================================================================

# Generate full table annotated with expression and stats
FULL <- Table_annotated[rownames(Table_annotated) %in% rownames(Stats),]

# Calculate mean expression for each celltype
for(n in names(table(METADATA$CellType))){
  subx <- FULL[,METADATA$SampleName[METADATA$CellType %in% n]]
  mean <- rowMeans(subx)
  FULL <- cbind(FULL, newcol = mean)
  colnames(FULL)[ncol(FULL)] <- paste0('mean_', n)
}

# Add statistics 
FULL <- cbind(FULL, 
              ALP_vs_EILP_logFC = Stats$logFC,
              ALP_vs_EILP_pval = Stats$P.Value, 
              ALP_vs_EILP_padj = Stats$adj.P.Val)

# Save final table
write.table(FULL, paste0(PATH_SAVE, '/FULL_Table_Results.txt'), 
            quote = F, sep = '\t')


#===============================================================================
# ANALYZE AND PLOT RESULTS -----------------------------------------------------
#===============================================================================

FULL   <- read.table(paste0(PATH_SAVE, '/FULL_Table_Results.txt'))
Stats  <- read.table(paste0(PATH_SAVE, '/LimmaStats_sEILP_vs_ALP.txt'))
DEG_UP <- Stats$Symbol[Stats$P.Value < 0.05 & Stats$logFC > 2]
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
         color = colorRampPalette(c("blue","white","red"))(50), main = "", silent=F)
dev.print(png, file=paste0(PATH_FIG, '/TF_dev_ILC.png'), 
          width=9, height=9, units='in', res=100)







