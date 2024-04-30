#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================




#===============================================================================
## SET UP ----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

PATH      <- 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC'
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/KO')
setwd(PATH)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

source("C:/Users/E15639P/Doctorat/scRNA-seq/Usefull_functions.R")

data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime_ALL.rds'))
red = 'rev_umap_cssz2'



#===============================================================================
## LEAP MATRIX PREPARATION -----------------------------------------------------
#===============================================================================

DimPlot(data, reduction = red, group.by='orig.ident', cols = ColorBlind[c(10,1,9,3,8,6)],
        pt.size = 1)
gradient         <- hcl.colors(100)

## OVERLIP WT ------------------------------------------------------------------

overWT <- subset(data, orig.ident %in% 'overLIP' & Slingshot_2 %!in% NA)
FeaturePlot(overWT, 'Slingshot_2', pt.size = 2, reduction = red) + 
  labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient)

result_S1 <- LeapMatrix(overWT, trajectory = "Slingshot_2", suffix = "overWT",  min.cells = 50, write = F, dir = paste0(PATH_SAVE, '/LEAP/overWT'))


## OVERLIP NF ------------------------------------------------------------------

overNF <- subset(data, orig.ident %in% 'NF-KO' & Slingshot_2 %!in% NA)
FeaturePlot(overNF, 'Slingshot_2', pt.size = 2, reduction = red) + 
  labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient)

result_S1 <- LeapMatrix(overNF, trajectory = "Slingshot_2", suffix = "overNF",  min.cells = 50, write = F, dir = paste0(PATH_SAVE, '/LEAP/overNF'))


## OVERLIP TOX ------------------------------------------------------------------

overTOX <- subset(data, orig.ident %in% 'TOX-KO' & Slingshot_2 %!in% NA)
FeaturePlot(overTOX, 'Slingshot_2', pt.size = 2, reduction = red) + 
  labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient)

result_S1 <- LeapMatrix(overNF, trajectory = "Slingshot_2", suffix = "overTOX",  min.cells = 50, write = F, dir = paste0(PATH_SAVE, '/LEAP/overTOX'))



#===============================================================================
## GENE EXPRESSION VISUALIZATION -----------------------------------------------
#===============================================================================

## WT vs NF-KO -----------------------------------------------------------------
NF <- subset(data, orig.ident %!in% 'TOX-KO')

for(i in 1:length(NF@meta.data[["orig.ident"]])){
  if(NF@meta.data[["orig.ident"]][i] %!in% 'NF-KO'){
    NF@meta.data[["Sample"]][i] <- 'WT'
  }else{
    NF@meta.data[["Sample"]][i] <- 'NF-KO'
  }
}

DimPlot(NF, reduction = red, group.by='Sample', cols = ColorBlind[c(6,11)],
        pt.size = 1)

table_NF   <- OrderMatrix(NF, "Slingshot_2", min.cells = 50, slot = 'scaledata')


# Add Sample row 
Sample <- c()
Ann <- data.frame(Cell.ID = NF@assays[["RNA"]]@counts@Dimnames[[2]], Sample = NF@meta.data[["Sample"]])
for(i in 1:ncol(table_NF)){
  x <- ifelse(Ann$Sample[Ann$Cell.ID %in% colnames(table_NF)[i]] %in% 'NF-KO', 1, 0)
  Sample <- c(Sample, x)
}
table_NF <- rbind(Sample = as.vector(Sample),table_NF)

gene = 'Nfil3'
test <- CompareExpr(x = table_NF, bin.number = 11, feature = gene, by.order = F, 
                    scale = F, col = ColorBlind[c(16, 6)], main = gene, ylim = c(-0.5,1.5))


test <- subset(NF, orig.ident %!in% 'NF-KO')


## ALL overLIPs ----------------------------------------------------------------


# Initialize
overALL <- subset(data, orig.ident %in% c('overLIP', 'NF-KO', 'TOX-KO'))

DimPlot(overALL, reduction = red, group.by='orig.ident', cols = ColorBlind[c(3,6,8)],
        pt.size = 1)
DimPlot(overALL, reduction = red, group.by='seurat_clusters', cols = ColorBlind[c(10,1,9,3,8,6)],
        pt.size = 1)

table   <- OrderMatrix(overALL, "Slingshot_2", min.cells = 50, slot = 'scaledata')

# Add Sample row 
Sample <- c()
Ann <- data.frame(Cell.ID = overALL@assays[["RNA"]]@counts@Dimnames[[2]], Sample = overALL@meta.data[["orig.ident"]])
for(i in 1:ncol(table)){
  if(Ann$Sample[Ann$Cell.ID %in% colnames(table)[i]] %in% 'NF-KO'){
    Sample <- c(Sample, 1)
  }else if(Ann$Sample[Ann$Cell.ID %in% colnames(table)[i]] %in% 'TOX-KO'){
    Sample <- c(Sample, 2)
  }else{
    Sample <- c(Sample, 0)
  }
}
table <- rbind(Sample = as.vector(Sample), table)

gene = 'Nfil3'
test <- CompareExpr(x = table, bin.number = 11, feature = gene, by.order = F, 
                    scale = F, col = ColorBlind[c(3,6,8)], main = gene, ylim = c(-0.5,1.75))

# Draw all candidates
Candidates.UP   <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_UP.txt')[,1]
Candidates.DOWN <- read.table('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Saves/Candidates_Cross_DOWN.txt')[,1]

for(g in Candidates.UP){
  #pdf(paste0('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Figures/scRNA-seq/KO/Candidates/UP/', g, '.pdf'),
      #width = 9.00, height = 9.00)
  CompareExpr(x = table, bin.number = 11, feature = g, by.order = F, 
              scale = F, col = ColorBlind[c(3,6,8)], main = g, ylim = c(-0.5,1.75))
  #dev.off()
}

for(g in Candidates.DOWN){
  #pdf(paste0('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Figures/scRNA-seq/KO/Candidates/DOWN/', g, '.pdf'),
  #    width = 9.00, height = 9.00)
  CompareExpr(x = table, bin.number = 11, feature = g, by.order = F, 
              scale = F, col = ColorBlind[c(3,6,8)], main = g, ylim = c(-0.5,1.75))
  #dev.off()
}


## Candidates
TF_list <- read.table('C:/Users/E15639P/Data/Lists/masterTFlist.txt', header = T)[,1]
List <- c('Id2', 'Itgb7', 'Gata3', 'Nfkb1', 'Pdcd1', 'Runx3', 'Sox5', 'Tcf7', 'Tox',
          'Zbtb16', 'Hoxa9', 'Irf8', 'Lmo2', 'Zeb2')

for(g in c(Candidates.UP, Candidates.DOWN, 'Nfil3')){
  if(g %in% c(TF_list, List)){
      pdf(paste0('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Figures/scRNA-seq/KO/Candidates/TF/', g, '.pdf'),
          width = 9.00, height = 9.00)
      CompareExpr(x = table, bin.number = 11, feature = g, by.order = F, 
                  scale = F, col = ColorBlind[c(3,6,8)], main = g, ylim = c(-0.5,1.75))
      dev.off()
  }
  
}



#===============================================================================
## LEAP RESULTS ANALYSIS -------------------------------------------------------
#===============================================================================

## ANNOTATE ALL RESULTS --------------------------------------------------------

alpha        <- 0.05
LEAP_Results <- list()

for(r in c('overNF', 'overTOX', 'overWT')){
  # Read files
  res   <- read.table(paste0(PATH_SAVE, '/LEAP/', r, '/MAC_results_', r, '.txt'), header = T)
  index <- read.table(paste0(PATH_SAVE, '/LEAP/', r, '/gene_index_', r, '.txt'), header = T)
  perm  <- read.table(paste0(PATH_SAVE, '/LEAP/', r, '/perm_', r, '.csv'), header = T, sep = ',')
  
  # Get significant threshold using Perm
  perm  <- perm[perm$fdr < 0.05,]
  thr   <- perm$cors[nrow(perm)]
  
  # Filter and save res
  res   <- res[res$Correlation >= thr,]
  #LEAP_Results[[r]] <- AnnotateLeap(res, index, write = T, dir = paste0(PATH_SAVE, '/LEAP/', r))
  res   <- NULL
}


## GET NF-SPECIFIC NETWORK -----------------------------------------------------

# NF regulator
NF_in_TOX <- LEAP_Results$overTOX[LEAP_Results$overTOX$gene_row %in% 'Nfil3',]
NF_in_NF <- LEAP_Results$overNF[LEAP_Results$overNF$gene_row %in% 'Nfil3',]
NF_in_WT <- LEAP_Results$overWT[LEAP_Results$overWT$gene_row %in% 'Nfil3',]

# NF regulated
NF_in_TOX <- LEAP_Results$overTOX[LEAP_Results$overTOX$gene_col %in% 'Nfil3',]
NF_in_NF <- LEAP_Results$overNF[LEAP_Results$overNF$gene_col %in% 'Nfil3',]
NF_in_WT <- LEAP_Results$overWT[LEAP_Results$overWT$gene_col %in% 'Nfil3',]



#===============================================================================
## GENE EXPRESSION VISUALIZATION -----------------------------------------------
#===============================================================================

overALL <- subset(data, orig.ident %in% c('overLIP', 'NF-KO', 'TOX-KO'))

table   <- OrderMatrix(overALL, "Slingshot_2", min.cells = 50)
table   <- rbind(Sample = unlist(strsplit(colnames(table), '_'))[c(T,F)], table)

gene

test <- CompareExpr(x = table, bin.number = 20, feature = 'Nfil3', 
                    by.order=F, scale = T, col = ColorBlind[c(14,16,15)])





################################################################################
################################### STOP HERE ##################################
################################################################################


################################################################################
### GENE EXPRESSION VISUALIZATION ----------------------------------------------
################################################################################

table <- OrderMatrix(data, "Slingshot_1", min.cells = 50)

#gene_list = 'Cd74'
#DrawExpr(table, bin.number = 9, by.order = T, scale = T, feature = gene_list, ylim = c(0, 1), std = T, col = c("darkblue","yellow", "red", "orange"), compare.with="Nfil3", main = "test",
#          write = F, dir = "")


#gene_list <- read.table('C:/Users/E15639P/Downloads/Liste_gènes.txt')[,1]

## FUNCTIONNAL !
DrawExpr(table, bin.number = 9, by.order = T, feature = c("Il7r", "Nfil3", "Tox", "Tcf7", "Id2"), scale = T, ylim = c(0,1), col = c("#26C4EC","#357AB7", "#22427C"), lwd = 3,
         superposed = T, compare.with = F, write = F, dir = "C:/Users/E15639P/Desktop/plot_florent")
##

gene_list <- "Nfil3"

# Working from provided gene_list and saving resulting plots

output_dir <- "C:/Users/E15639P/Desktop/plot_florent"

absent_genes <- c()
for(g in gene_list){
  if(g %in% rownames(table)){
    png(paste(output_dir, "/", g, ".png", sep = ""))
    DrawExpr(table, bin.number = 30, by.order = F, feature = g, ylim = c(0, 1), std = F, col = "darkblue")
    dev.off()
  }else{
    absent_genes <- c(absent_genes, g)
  }
}




 # # # # # # # # # # # # # # # # #
gene_list <- c("Tcf7", "Tox", "Tox2", "Bcl11b", "Thy1", "Dntt", "Il7r", "Cd74", "Id2", "Irf8", "Sp1", "Gata3")

gene_all <- c("Nfil3")
n_cells  <- 3456 


table <- OrderMatrix(data, "Slingshot_2", min.cells = 10)
DrawExpr(table, bin.number = 21, feature = "Uba52", ylim = c(0, 3))


## SLOT = data -----------------------------------------------------------------

## Trajectory Slingshot_2

# Get corresponding ordered matrix with pseudotime kept
table <- OrderMatrix(data, "Slingshot_2", min.cells = 10)

for(g in gene_list){
  #png(filename = paste(PATH, "/Figures/Step_03/S2_data/", g, ".png", sep = ""), width = 600, height = 600)
  DrawExpr(table, bin.number = 21, feature = g)
  #dev.off()
}

# For more precise curves
for(g in gene_all){
  #png(filename = paste(PATH, "/Figures/Step_03/all_cells/", g, ".png", sep = ""), width = 600, height = 600)
  DrawExpr(table, bin.number = 800, feature = g) #bin.number = 21
  #dev.off()
}

## Trajectory Monocle_2

# Get corresponding ordered matrix with pseudotime kept
table <- OrderMatrix(data, "Slingshot_1", min.cells = 100)

for(g in gene_list){
  #png(filename = paste(PATH, "/Figures/Step_03/M2_min100/", g, ".png", sep = ""), width = 600, height = 600)
  DrawExpr(table, bin.number = 21, feature = g)
  #dev.off()
}



################################################################################
### LEAP RESULTS ANALYSIS ------------------------------------------------------
################################################################################

## ILC LINEAGE -----------------------------------------------------------------

# Results annotation
S1.results   <- read.table(paste(PATH, '/Saves/Step_03/20230511/S1/MAC_results_S1.txt', sep = ''), check.names = F, header = T)
S1.index     <- read.table(paste(PATH, '/Saves/Step_03/20230511/S1/gene_index_S1.txt', sep = ''), check.names = F, header = T)
S1.annotated <- AnnotateLeap(S1.results, S1.index)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#write.table(S1.annotated, paste(PATH, '/Saves/Step_03/20230511/S1/MAC_results_S1_annotated.txt', sep = ''), row.names = F)
S1.annotated <- read.table(paste(PATH, '/Saves/Step_03/20230511/S1/MAC_results_S1_annotated.txt', sep = ''), check.names = F, header = T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


# Subset on NFIL3
NF <- S1.annotated[S1.annotated$gene_row == "Nfil3",]
NF[NF$gene_col == "Tox",]
NF <- NF[abs(NF$Correlation) >= 0.1125, ] # fdr < 0.5
NF <- NF[abs(NF$Correlation) >= 0.1625, ] # fdr < 0.05
NF <- NF[abs(NF$Correlation) >= 0.18, ] # fdr < 0.01

#write.table(NF$gene_col, paste(PATH, '/Saves/Step_03/20230511/S1/NFIL3_Candidates_fdr_0.01.txt', sep = ''), quote = F, col.names = F, row.names = F)

NF[NF$gene_col == "Zeb2",]


## DC LINEAGE ------------------------------------------------------------------

# Results annotation
S2.results   <- read.table(paste(PATH, '/Saves/Step_03/20230102/S2/MAC_results_S2.txt', sep = ''), check.names = F, header = T)
S2.index     <- read.table(paste(PATH, '/Saves/Step_03/20230102/S2/gene_index_S2.txt', sep = ''), check.names = F, header = T)
S2.annotated <- AnnotateLeap(S2.results, S2.index)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#write.table(S2.annotated, paste(PATH, '/Saves/Step_03/20230102/S2/MAC_results_S2_annotated.txt', sep = ''), row.names = F)
S2.annotated <- read.table(paste(PATH, '/Saves/Step_03/20230102/S2/MAC_results_S2_annotated.txt', sep = ''), check.names = F, header = T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


# Subset on NFIL3
NF <- S1.annotated[S1.annotated$gene_row == "Nfil3",]
NF[NF$gene_col == "Tox",]
#NF <- NF[abs(NF$Correlation) >= 0.1125, ]



NF[NF$gene_col == "Zeb2",]


#===============================================================================
## SCRNA2BONI TABLE PREPARATION  -----------------------------------------------
#===============================================================================

# Create matrix with evalues
BoNI  <- as.data.frame(t(as.data.frame(data@assays[["RNA"]]@data)))

# Get metadata needed for BoNI
BoNI2 <- data.frame(Name=rownames(BoNI),
                    Embryo=data@meta.data$orig.ident,
                    clusterUmap=data@meta.data$seurat_clusters, 
                    Pseudotime=data@meta.data$Slingshot_1, 
                    Dataset=data@meta.data$Experiment, 
                    EmbryoDay=rep(1, nrow(BoNI)), 
                    Stage=data@meta.data$orig.ident,
                    clusterEmbryoCells=data@meta.data$seurat_clusters)

# Merge both table and filter them on S1 traj
BoNI_full <- cbind(BoNI2, BoNI)
BoNI_full <- BoNI_full[!is.na(BoNI_full$Pseudotime),]

# Write csv table for BoNI
#write.csv(BoNI_full, paste0(PATH, '/Saves/Step_03/20230511/S1/BoNI_matrix.csv'), row.names = F)



################################################################################
# SCALE DATA / DATA COMPARISON ------------------------------------------------
################################################################################

results <- read.table(paste(PATH, "/Saves/Step_03/S2_data_20221027/results.txt", sep = ""))
results_scaled <- read.table(paste(PATH, "/Saves/Step_03/S2_scaled_compare_20221027/results.txt", sep = ""), row.names = NULL)
results_scaled <- results_scaled[,c(2,3,4,5)]
index   <- read.table(paste(PATH, "/Saves/Step_03/S2_data_20221027/gene_index_20221027_S2.txt", sep = ""))

annotated <- AnnotateLeap(results, index)
annotated_scaled <- AnnotateLeap(results_scaled, index)

#write.table(annotated, paste(PATH, "/Saves/Step_03/S2_data_20221027/results_indexed.txt", sep = ""), row.names = FALSE)
#write.table(annotated_scaled, paste(PATH, "/Saves/Step_03/S2_scaled_compare_20221027/results_indexed.txt", sep = ""), row.names = FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

annotated <- read.table(paste(PATH, "/Saves/Step_03/S2_data_20221027/results_indexed.txt", sep = ""), header = TRUE)
annotated_scaled <- read.table(paste(PATH, "/Saves/Step_03/S2_scaled_compare_20221027/results_indexed.txt", sep = ""), header = TRUE)

# Subset on NFIL3
NF <- annotated[annotated$gene_row == "Nfil3",]
NF_scaled <- annotated_scaled[annotated_scaled$gene_row == "Nfil3",]

# Apply the correlation cutoff defined by respective Perm_Matrix
NF <- NF[abs(NF$Correlation) >= 0.175,] #0.165 = 0.094 FDR / 0.175 = 0.046 FDR 
NF_scaled <- NF_scaled[abs(NF_scaled$Correlation) >= 0.1975,] #0.1875 = 0.044 FDR

# Data slot is better



#===============================================================================
## WGCNA RESULTS ---------------------------------------------------------------
#===============================================================================

TF_list <- read.table("C:/Users/E15639P/Data/Lists/masterTFlist.txt", header = T)[,1]
gene.id <- read.table(paste(PATH, "/Saves/Step_03/20230511/S1/gene_index_S1.txt", sep = ""), check.names = F, header = T)[,2]
wgcna <- read.csv(paste(PATH, "/Saves/Step_03/20230511/S1/MAC_symmetric_S1.csv", sep = ""))

row.names(wgcna) <- gene.id
colnames(wgcna) <- gene.id

sub_matrix <- wgcna[rownames(wgcna) %in% TF_list, colnames(wgcna) %in% TF_list]

write.table(sub_matrix, paste(PATH, "/Saves/Step_03/20230511/S1/WGCNA/GRN_TF_matrix.txt", sep = ""), quote = F)

test <- read.table(paste(PATH, "/Saves/Step_03/20230511/S1/WGCNA/GRN_TF_matrix.txt", sep = ""))





#   /)/)   zZ
#  ( u.u) z
# *(__()()






