#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
# Gene Regulatory Network analysis based on WT ILC developmental model
#
# Prepare matrix of cells ordered by pseudotime
# Launch LEAP algorithm on a calculation cluster
# Explore LEAP results
# Prepare the resulting gene regulatory network in a .gml file for Gephi  



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
suppressPackageStartupMessages(library(igraph))

# Create required directories
dir.create(file.path(paste0(PATH, '/Saves/scRNA-seq/WT'), 'LEAP'))
dir.create(file.path(paste0(PATH, '/Saves/scRNA-seq/WT'), 'Kmeans'))
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq/WT'), 'Curves'))
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq/WT'), 'Kmeans'))
dir.create(file.path(paste0(PATH, '/Figures/scRNA-seq/WT/Curves'), 'Targets'))

# Set up path for figures and saves
PATH_FIG  <- paste0(PATH, '/Figures/scRNA-seq/WT')
PATH_SAVE <- paste0(PATH, '/Saves/scRNA-seq/WT')



#===============================================================================
## INPUT -----------------------------------------------------------------------
#===============================================================================

data     <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime.rds'))



#===============================================================================
## LEAP MATRIX PREPARATION -----------------------------------------------------
#===============================================================================

## Generate matrix and index for ILC trajectory
result_S1 <- LeapMatrix(data, trajectory = 'Slingshot_1', suffix = 'WT',  
                        min.cells = 50, write = T, 
                        dir = paste0(PATH_SAVE, '/LEAP'), slot = 'data')
index  <- as.data.frame(result_S1[2])
matrix <- as.data.frame(result_S1[1])


################################################################################
#    READY TO USE GENERATED TABLE FOR LEAP ANALYSIS ON CALCULATION CLUSTER     #
################################################################################
#
# MAC_results <- MAC_counter(data = as.matrix(matrix), max_lag_prop = 1/3, 
#                            MAC_cutoff = 0.1, lag_matrix = TRUE, 
#                            file_name = 'WT')
# write.table(x = MAC_results, 
#             paste0(PATH, '/', 'MAC_results_', SUFFIX, '.txt'), 
#             row.names = FALSE)
#
# MAC_perm(data = as.matrix(matrix), MACs_observ = MAC_matrix, num_perms = 100, 
#          max_lag_prop = 1/3, FDR_cutoffs = 401, perm_file_name = 'WT')
#
################################################################################
#                                                                              #
################################################################################


# Read LEAP results from computation cluster 
MAC_results <- read.table(paste0(PATH_SAVE, '/LEAP/MAC_results_WT.txt'), 
                          header = T)
gene_index  <- read.table(paste0(PATH_SAVE, '/LEAP/gene_index_WT.txt'), 
                          header = T)
perm        <- read.table(paste0(PATH_SAVE, '/LEAP/perm_WT.csv'), sep = ',', 
                          header =T)

Annotate    <- AnnotateLeap(MAC_results, gene_index, write = T, 
                            dir = paste0(PATH_SAVE, '/LEAP'), 
                            filename = 'results_indexed')

# Apply correlation cutoff defined by permutation
perm_cutoff <- 0.1625
Annotate    <- Annotate[abs(Annotate$Correlation) > perm_cutoff,]


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save filtered interactions
write.table(Annotate, paste0(PATH_SAVE, '/LEAP/results_indexed_filtered.txt'),
            row.names = F, col.names = T, quote = F)
Annotate <- read.table(paste0(PATH_SAVE, '/LEAP/results_indexed_filtered.txt'),
                       header = T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## LEAP RESULTS ANALYSIS -------------------------------------------------------
#===============================================================================

NF      <- Annotate[Annotate$gene_row %in% 'Nfil3' & Annotate$Lag > 0,]
NF_UP   <- NF$gene_col[NF$Correlation > 0]
NF_DOWN <- NF$gene_col[NF$Correlation < 0]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save table and candidates
write.table(NF, paste0(PATH_SAVE, '/LEAP/Supplementary_Table_2.txt'), 
            quote = F, row.names = F, col.names = T)
write.table(NF_UP, paste0(PATH_SAVE, '/LEAP/LEAP_Candidates_UP.txt'), 
            quote = F, row.names = F, col.names = F)
write.table(NF_DOWN, paste0(PATH_SAVE, '/LEAP/LEAP_Candidates_DOWN.txt'), 
            quote = F, row.names = F, col.names = F)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#===============================================================================
## NETWORK ---------------------------------------------------------------------
#===============================================================================

Annotate <- read.table(paste0(PATH_SAVE, '/LEAP/results_indexed_filtered.txt'), 
                       header = T)

# Generate network for TF known to be involved during ILC development
Factors <- c('Nfil3', 'Id2', 'Tox', 'Gata3', 'Tcf7', 'Zbtb16')

Network <- Annotate[Annotate$gene_row %in% Factors & 
                      Annotate$gene_col %in% Factors,]
Network <- Network[Network$Lag %!in% 0,]
Network <- Network[Network$Correlation > perm_cutoff,]


## iGraph generation for Candidates only and saving to GML file
Genes      <- unique(c(Network$gene_row, Network$gene_col))

nodes <- data.frame(Name   = Genes)
edges <- data.frame(From   = Network$gene_row,
                    To     = Network$gene_col,
                    Cor    = Network$Correlation,
                    AbsCor = abs(Network$Correlation),
                    DirCor = ifelse(Network$Correlation > 0, 'Up', 'Down'),
                    Lag    = Network$Lag,
                    Neglag = 1.5-as.numeric(Network$Lag/max(Network$Lag))) 

graph  <- graph_from_data_frame(edges, directed = T, vertices = nodes)
plot(graph, layout = layout.fruchterman.reingold, vertex.size = 45,
     edge.arrow.size = 0.5, vertex.label=Genes, 
     vertex.color = ColorBlind[3], vertex.label.color='black',)
plot19 <- recordPlot()
writePlot(plot19, PATH_FIG)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Save Network files for GEPHI visualization 
write_graph(graph, paste0(PATH_SAVE, '/Network.gml'), format='gml')
write.table(Network, paste0(PATH_SAVE, '/LEAP/Supplementary_Table_1.txt'),
            quote = F, row.names = F, col.names = T)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


################################################################################
#      READY FOR DEEPER ANALYSIS OF GENE REGULATORY NETWORK USING GEPHI        #
################################################################################

                                                                        
