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
suppressPackageStartupMessages(library(slingshot))

source("C:/Users/E15639P/Doctorat/scRNA-seq/Usefull_functions.R")
CellTypeOrder     <- c("ALP", "TULIP", "ILCpro", "overLIP", "NF-KO", "TOX-KO")

#data <- readRDS(paste0(PATH_SAVE, '/4_Filtered.rds'))
#red = 'rev_umap_cssz2'
#xlimit <- c(min(data@reductions[[red]]@cell.embeddings[,1]),
#            max(data@reductions[[red]]@cell.embeddings[,1]))
#ylimit <- c(min(data@reductions[[red]]@cell.embeddings[,2]),
#            max(data@reductions[[red]]@cell.embeddings[,2]))      


#===============================================================================
## SLINGSHOT : NF --------------------------------------------------------------
#===============================================================================
#Based on : https://bustools.github.io/BUS_notebooks_R/slingshot.html

data <- readRDS(paste0(PATH_SAVE, '/4_Filtered_NF.rds'))
red = 'rev_umap_cssz2'
xlimit <- c(min(data@reductions[[red]]@cell.embeddings[,1]),
            max(data@reductions[[red]]@cell.embeddings[,1]))
ylimit <- c(min(data@reductions[[red]]@cell.embeddings[,2]),
            max(data@reductions[[red]]@cell.embeddings[,2]))  

DimPlot(data, reduction = red, group.by = "seurat_clusters", 
        label = TRUE, repel = TRUE, 
        cols = c(ColorBlind[c(3,9,2,7,8,10)], "purple"), pt.size = 0.9) +
  labs(title = "", x = "UMAP1", y = "UMAP2")

DimPlot(data, group.by = "orig.ident", reduction = "rev_umap_cssz2", cols = ColorBlind[c(11,11,11,3,6)], pt.size = 0.9, order = rev(CellTypeOrder)) +
  labs(title = "", x = "UMAP1", y = "UMAP2")


## PSEUDOTIME RECONSTRUCTION ---------------------------------------------------
sce                                     <- as.SingleCellExperiment(data)     
colData(sce)$Seurat_clusters            <- as.character(data@active.ident)

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = toupper(red), start.clus = 0, end.clus = c(5,3))


## ADDING PSEUDOTIMES TO DATA --------------------------------------------------
pseudotime                              <- slingPseudotime(sce)
for(i in 1:ncol(pseudotime)){
  data@meta.data[[paste("Slingshot_", i, sep = "")]] <- pseudotime[,i]
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Saving SCE file after processing
#saveRDS(sce, paste0(PATH_SAVE, '/Slingshot_NF.rds'))
sce <- readRDS(paste0(PATH_SAVE, '/Slingshot_NF.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## PLOTTING PSEUDOTIME ---------------------------------------------------------

# Specifying gradient colors for plots
gradient         <- hcl.colors(100)

# First Trajectory
plot1 <- FeaturePlot(data, 'Slingshot_1', pt.size = 2, reduction = red) + labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient) +
  xlim(xlimit) + ylim(ylimit)
plot1


# Second
plot2 <- FeaturePlot(data, 'Slingshot_2', pt.size = 2, reduction = red) + labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient) +
  xlim(xlimit) + ylim(ylimit)
plot2


# Plot with curves
curves <- slingCurves(sce, as.df = TRUE) %>%
  dplyr::rename("x" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_1'), "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

plot2 <- FeaturePlot(data, 'Slingshot_1', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), size = 1.5) +
  labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
plot2



## TEST LINEAR REPRESENTATION

compare <- data.frame(Slingshot = data$Slingshot_1, Sample = data$orig.ident)
compare <- compare[order(compare$Slingshot),]
compare <- compare[!is.na(compare$Slingshot),]
compare$Slingshot_order <- c(1:nrow(compare))


ggplot(compare, aes(x=Slingshot, fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(alpha=0.5, linewidth = 1) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,8,6)], name = "")

ggplot(data=compare, aes(x=Slingshot, fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(adjust=1.5, position="fill", linewidth = 0.5) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,8,6)], name = "")
#dev.print(png, file=paste(PATH, "/Figures/Step_02/Density_Plot.png", sep = ""), width=12, height=4, units="in", res=200)


ggplot(compare, aes(x=Slingshot, y = 1, colour=factor(Sample, levels = CellTypeOrder))) +
  geom_point(shape=124) +
  theme_classic()+
  scale_color_manual(values = ColorBlind[c(11,11,11,3,6)], name = "")




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Saving CDS file after processing and object data after adding metadata
#DO NOT SAVE CDS FROM MONOCLE3 < 20GB
#saveRDS(data, paste0(PATH_SAVE, '/5_Pseudotime_NF.rds'))
#data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime_NF.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =







#===============================================================================
## SLINGSHOT : ALL -------------------------------------------------------------
#===============================================================================
#Based on : https://bustools.github.io/BUS_notebooks_R/slingshot.html

data <- readRDS(paste0(PATH_SAVE, '/4_Filtered_ALL.rds'))
red = 'rev_umap_cssz2'
xlimit <- c(min(data@reductions[[red]]@cell.embeddings[,1]),
            max(data@reductions[[red]]@cell.embeddings[,1]))
ylimit <- c(min(data@reductions[[red]]@cell.embeddings[,2]),
            max(data@reductions[[red]]@cell.embeddings[,2]))      

DimPlot(data, reduction = red, group.by = "seurat_clusters", 
        label = TRUE, repel = TRUE, 
        cols = c(ColorBlind[c(3,9,2,7,8,10)], "purple"), pt.size = 0.9) +
  labs(title = "", x = "UMAP1", y = "UMAP2")


## PSEUDOTIME RECONSTRUCTION ---------------------------------------------------
sce                                     <- as.SingleCellExperiment(data)     
colData(sce)$Seurat_clusters            <- as.character(data@active.ident)

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = toupper(red), start.clus = 0, end.clus = c(4,3))


## ADDING PSEUDOTIMES TO DATA --------------------------------------------------
pseudotime                              <- slingPseudotime(sce)
for(i in 1:ncol(pseudotime)){
  data@meta.data[[paste("Slingshot_", i, sep = "")]] <- pseudotime[,i]
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Saving SCE file after processing
#saveRDS(sce, paste0(PATH_SAVE, '/Slingshot_ALL.rds'))
sce <- readRDS(paste0(PATH_SAVE, '/Slingshot_ALL.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## PLOTTING PSEUDOTIME ---------------------------------------------------------

# Specifying gradient colors for plots
gradient         <- hcl.colors(100)

# First Trajectory
plot1 <- FeaturePlot(data, 'Slingshot_1', pt.size = 2, reduction = red) + labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient) +
  xlim(xlimit) + ylim(ylimit)
plot1


# Second
plot2 <- FeaturePlot(data, 'Slingshot_2', pt.size = 2, reduction = red) + labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient) +
  xlim(xlimit) + ylim(ylimit)
plot2


# Plot with curves
curves <- slingCurves(sce, as.df = TRUE) %>%
  dplyr::rename("x" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_1'), "y" = paste0(str_replace_all(toupper(red), "[[:punct:]]", ""), '_2'))

plot2 <- FeaturePlot(data, 'Slingshot_2', pt.size = 2, reduction = red) + 
  geom_path(data = curves %>% arrange(Order), aes(x, y, group = Lineage), size = 1.5) +
  labs(title = "", x = "UMAP1", y = "UMAP2", pt.size = 0.9) + 
  scale_color_gradientn(name = "Pseudotime", colours = gradient) +
  xlim(xlimit) + ylim(ylimit) 
plot2



## TEST LINEAR REPRESENTATION

compare <- data.frame(Slingshot = data$Slingshot_2, Sample = data$orig.ident)
compare <- compare[order(compare$Slingshot),]
compare <- compare[!is.na(compare$Slingshot),]
compare$Slingshot_order <- c(1:nrow(compare))


ggplot(compare, aes(x=Slingshot, fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(alpha=0.5, linewidth = 1) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,8,6)], name = "")

ggplot(data=compare, aes(x=Slingshot, fill=factor(Sample, levels = CellTypeOrder))) +
geom_density(adjust=1.5, position="fill", linewidth = 0.5) +
theme_classic() +
scale_fill_manual(values = ColorBlind[c(10,1,9,3,8,6)], name = "")
#dev.print(png, file=paste(PATH, "/Figures/Step_02/Density_Plot.png", sep = ""), width=12, height=4, units="in", res=200)


ggplot(compare, aes(x=Slingshot, y = 1, colour=factor(Sample, levels = CellTypeOrder))) +
  geom_point(shape=124) +
  theme_classic()+
  scale_color_manual(values = ColorBlind[c(11,11,11,3,8,6)], name = "")



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Saving CDS file after processing and object data after adding metadata
#DO NOT SAVE CDS FROM MONOCLE3 < 20GB
#saveRDS(data, paste0(PATH_SAVE, '/5_Pseudotime_ALL.rds'))
data <- readRDS(paste0(PATH_SAVE, '/5_Pseudotime_ALL.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## Subset for overLIPs line ----------------------------------------------------

overALL  <- subset(data, orig.ident %in% c('overLIP', 'NF-KO', 'TOX-KO'))
max_time <- max(overALL$Slingshot_2[!is.na(overALL$Slingshot_2)]) 

data <- subset(data, Slingshot_2 <= max_time)

## TEST LINEAR REPRESENTATION

compare <- data.frame(Slingshot = data$Slingshot_2, Sample = data$orig.ident)
compare <- compare[order(compare$Slingshot),]
compare <- compare[!is.na(compare$Slingshot),]
compare$Slingshot_order <- c(1:nrow(compare))


ggplot(compare, aes(x=Slingshot, fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(alpha=0.5, linewidth = 1) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,8,6)], name = "")

ggplot(data=compare, aes(x=Slingshot, fill=factor(Sample, levels = CellTypeOrder))) +
  geom_density(adjust=1.5, position="fill", linewidth = 0.5) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind[c(10,1,9,3,8,6)], name = "")
#dev.print(png, file=paste(PATH, "/Figures/Step_02/Density_Plot.png", sep = ""), width=12, height=4, units="in", res=200)


ggplot(compare, aes(x=Slingshot, y = 1, colour=factor(Sample, levels = CellTypeOrder))) +
  geom_point(shape=124) +
  theme_classic()+
  scale_color_manual(values = ColorBlind[c(11,11,11,3,6,8)], name = "")


#===============================================================================
## DC TRAJECTORIES -------------------------------------------------------------
#===============================================================================
# cDC2 are no longer present

DC <- subset(data, Slingshot_1 %!in% NA)
DimPlot(DC, reduction = red)

DC <- RunUMAP(DC, reduction = "cssz", dims = 1:ncol(Embeddings(data, "cssz")), n.neighbors = 50, reduction.name="DC", reduction.key="DC_")
DimPlot(DC, group.by = 'orig.ident', reduction = "DC", pt.size = 1, cols = ColorBlind[c(10,1,9,3,8,6)]) + scale_x_reverse() + scale_y_reverse()

FeaturePlot(DC, "Id2", reduction = "DC", pt.size = 1, order = T) + scale_x_reverse() + scale_y_reverse()

DC <- AddModuleScore(DC, DC_genelist, name = 'DC_score')

# Score spatial repartition in cells
plot_density(DC, 'DC_score1', pal = 'magma', size = 2) + scale_x_reverse() + scale_y_reverse()


################################################################################
################################ STOP HERE #####################################
################################################################################




#===============================================================================
## MONOCLE3 --------------------------------------------------------------------
#===============================================================================
#Based on : https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(monocle3))


## SETTING UP A CDS OBJECT -----------------------------------------------------
#Based on zihengxuwu answer : https://github.com/satijalab/seurat/issues/1658

# Gene annotations
gene_annotation <- as.data.frame(row.names(data@reductions[["pca"]]@feature.loadings), row.names = row.names(data@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# Cell metadata
cell_metadata <- as.data.frame(data@assays[["RNA"]]@data@Dimnames[[2]], row.names = data@assays[["RNA"]]@data@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# Counts sparse matrix
new_matrix        <- data@assays[["RNA"]]@data
new_matrix        <- new_matrix[rownames(data@reductions[["pca"]]@feature.loadings),]
expression_matrix <- new_matrix

# Construct the basic cds object
cds_from_seurat   <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Construct and assign the made up partition
recreate.partition        <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition        <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

# Assign the cluster information
list_cluster              <- data@meta.data[[sprintf("seurat_clusters")]]
names(list_cluster)       <- data@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

# Could be a space-holder, but essentially fills out louvain parameters
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

# Assign UMAP coordinate
cds_from_seurat@int_colData@listData[["reducedDims"]][["UMAP"]] <- data@reductions[["rev_umap_cssz2"]]@cell.embeddings

# Assign feature loading for downstream module analysis
cds_from_seurat[["preprocess_aux"]][["gene_loadings"]] <- data@reductions[["pca"]]@feature.loadings


## PSEUDO-TIME TRAJECTORY ------------------------------------------------------

#Use Monocle3 clustering function
cds   <- cluster_cells(cds_from_seurat, reduction_method = 'UMAP', cluster_method = 'leiden', resolution = 0.0002)

plot5 <- plot_cells(cds, color_cells_by = "cluster", group_label_size = 4)
plot5

#Learning graph from clusters
cds   <- learn_graph(cds)

#Plot trajectory and clustering by color
plot6 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE) +
  labs(x = "UMAP1", y = "UMAP2")
plot6

pseudotime <- order_cells(cds, reduction_method = "UMAP")


# Trajectory 1 subsetting
pseudotime1 <- choose_graph_segments(cds, clear_cds = FALSE)
plot7 <- plot_cells(pseudotime1, cell_size = 2)
plot7
pseudotime1 <- order_cells(pseudotime1, reduction_method = "UMAP")
plot8 <- plot_cells(cds = pseudotime1, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, cell_size = 1.5) + 
  labs(title = "", x = "UMAP1", y = "UMAP2") +
  scale_color_gradientn(name = "Pseudotime", colours = gradient)
plot8


# Trajectory 2 subsetting
pseudotime2 <- choose_graph_segments(cds, clear_cds = FALSE)
plot9 <- plot_cells(pseudotime2, cell_size = 2)
plot9
pseudotime2 <- order_cells(pseudotime2, reduction_method = "UMAP")
plot10 <- plot_cells(cds = pseudotime2, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, cell_size = 1.5) + 
  labs(title = "", x = "UMAP1", y = "UMAP2") +
  scale_color_gradientn(name = "Pseudotime", colours = gradient)
plot10


## ADDING PSEUDOTIMES TO DATA --------------------------------------------------

data$barcode <- row.names(data[[]])

times_list <- list(pseudotime1@principal_graph_aux@listData[["UMAP"]][["pseudotime"]],
                   pseudotime2@principal_graph_aux@listData[["UMAP"]][["pseudotime"]])

meta1 <- c()
meta2 <- c()

for(b in data$barcode){
  if(b %in% names(unlist(times_list[1]))){
    meta1 <- c(meta1, as.vector(unlist(times_list[1])[b]))
  } else{
    meta1 <- c(meta1, NA)
  }
  if(b %in% names(unlist(times_list[2]))){
    meta2 <- c(meta2, as.vector(unlist(times_list[2])[b]))
  } else{
    meta2 <- c(meta2, NA)
  }
}

data@meta.data[["Monocle_1"]] <- meta1
data@meta.data[["Monocle_2"]] <- meta2


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Saving CDS file after processing and object data after adding metadata
#DO NOT SAVE CDS FROM MONOCLE3 < 20GB
#saveRDS(data, paste0(PATH_SAVE, '/5_Pseudotime.rds'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## PLOTTING GLOBAL PSEUDOTIME --------------------------------------------------

gradient <- hcl.colors(100)

# Monocle plot cells with pseudotime
plot11 <- plot_cells(cds = pseudotime, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, cell_size = 1.5) + 
  labs(title = "", x = "UMAP1", y = "UMAP2") +
  scale_color_gradientn(name = "Pseudotime", colours = gradient)
plot11

plot12 <- plot_cells(pseudotime, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, cell_size = 1.5, trajectory_graph_segment_size = 1.5) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_color_gradientn(name = "Pseudotime", colours = gradient)
plot12



















################################################################################
### PSEUDOTIME ANALYSIS --------------------------------------------------------
################################################################################


data <- readRDS(paste(PATH, '/Saves/Step_02/CSS_4_pseudotime.rds', sep = ''))

compare <- data.frame(Slingshot = data$Slingshot_1, Sample = data$orig.ident)
compare <- compare[order(compare$Slingshot),]
compare <- compare[!is.na(compare$Slingshot),]
compare$Slingshot_order <- c(1:nrow(compare))


stats <- LinearStats(data$Monocle_1, data$Slingshot_1)

Posistats <- data.frame()
for(s in names(table(data$orig.ident))){
  mean <- mean(data$Slingshot_1[data$orig.ident == s])
  std  <- sd(data$Slingshot_1[data$orig.ident == s])
  Posistats <- rbind(Posistats,c(s, mean, std))
}
colnames(Posistats) <- c("sample", "mean", "std")



ggplot(compare, aes(y = 1, x = Slingshot_order)) +
  scale_shape_identity() +
  geom_point(mapping = aes(colour = factor(Sample, levels = c("ALP", "overLIP", "TULIP", "ILCpro")), shape = 73), size = 5) +
  scale_color_manual(values = ColorBlind12[c(6, 5, 1, 4)], name = "") + 
  scale_x_continuous(name = "Position au sein du Pseudotemps Slinghshot", limits = c(0, nrow(compare))) +
  scale_y_continuous(name = "", limits = c(0.5, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(0.5, 'cm')) +
  annotate("text", x = 300, y = 3700, label = stats[1], fontface = 2, colour = "black") +
  annotate("text", x = 500, y = 3550, label = stats[2], fontface = 2, colour = "black") +
  labs(title = "") 

ggplot() +
  geom_point(data = compare, aes(y = factor(Sample, levels = rev(c("ALP", "overLIP", "TULIP", "ILCpro"))), x = Slingshot_order, colour = factor(Sample, levels = c("ALP", "overLIP", "TULIP", "ILCpro"))), size = 1) +
  geom_point(data = Posistats, aes(y = factor(sample, levels = rev(c("ALP", "overLIP", "TULIP", "ILCpro"))),  x = round(as.numeric(mean)))) +
  geom_errorbarh(data = Posistats, aes(y = factor(sample, levels = rev(c("ALP", "overLIP", "TULIP", "ILCpro"))), xmax = as.numeric(mean) + as.numeric(std), xmin = as.numeric(mean) - as.numeric(std))) +
  scale_color_manual(values = ColorBlind12[c(6, 5, 1, 4)], name = "") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "", x = "", y ="") +
  NoLegend()










################################################################################

######################### NOT USED ANYMORE #####################################

################################################################################
### COMPARISON OF PSEUDOTIMES --------------------------------------------------
################################################################################


  
# Monocle3
pseudotime                       <- data.monocle@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
names                            <- data.monocle@colData@listData[["barcode"]]
monocle_pseudotime               <- data.frame(Cell = names, Monocle = pseudotime)
monocle_pseudotime               <- monocle_pseudotime[order(monocle_pseudotime$Monocle),]
monocle_pseudotime$Monocle_order <- c(1:length(names))

# Slingshot
pseudotime       <- data.slingshot@colData@listData[["slingPseudotime_1"]]
names            <- data.slingshot@colData@listData[["barcode"]]
sample           <- data.slingshot@colData@listData[["orig.ident"]]
sling_pseudotime <- data.frame(Cell = names, Slingshot = pseudotime, Sample = sample)
sling_pseudotime <- sling_pseudotime[order(sling_pseudotime$Slingshot),]
sling_pseudotime$Slingshot_order <- c(1:length(names))

# Merge the two dataframes based on barcode to compare pseudotimes and positions
compare <- merge(monocle_pseudotime, sling_pseudotime, by = "Cell")


## PLOTTING STRICT PSEUDOTIMES

# Linear correlation and p-value with Pearson's test
stats <- LinearStats(compare$Monocle, compare$Slingshot)
# R² = 0.9653 ans p-value < 2.2e-16

cor.test(compare$Monocle_order,compare$Slingshot_order, method="spearman")
#http://www.sthda.com/french/wiki/test-de-correlation-entre-deux-variables

plot8 <- ggplot(compare, aes(y = Monocle, x = Slingshot)) +
  geom_point(aes(colour = factor(Sample, levels = c("ALP", "overLIP", "Tulip", "ILCpro"))), size = 1.5) +
  geom_smooth(method="lm", colour = "black", fill = "red") +
  scale_color_manual(values = ColorBlind12[c(6, 5, 1, 4)], name = "") + 
  scale_x_continuous(name = "Pseudotemps Slingshot", limits = c(0, 22), breaks = seq(0,25,5)) +
  scale_y_continuous(name = "Pseudotemps Monocle3", limits = c(-1.5, 33), breaks = seq(0,30,5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(0.5, 'cm')) +
  labs(title = "", fill = " ") +
  annotate("text", x = 1.5, y = 32.5, label = stats[1], fontface = 2, colour = "black") +
  annotate("text", x = 2.55, y = 31, label = stats[2], fontface = 2, colour = "black")
plot8


## PLOTTING POSITIONS ACCROSS PSEUDOTIMES

# Linear correlation and p-value with Pearson's test
stats <- LinearStats(compare$Monocle_order, compare$Slingshot_order)
# R² = 0.9606 and p-value < 2.2e-16

plot9 <- ggplot(compare, aes(y = Monocle_order, x = Slingshot_order)) +
  geom_point(aes(colour = factor(Sample, levels = c("ALP", "overLIP", "Tulip", "ILCpro"))), size = 1.5) +
  geom_smooth(method="lm", colour = "black", fill = "red") +
  scale_color_manual(values = ColorBlind12[c(6, 5, 1, 4)], name = "") + 
  scale_x_continuous(name = "Position au sein du Pseudotemps Slinghshot", limits = c(0, 3750), breaks = seq(0,3800,500)) +
  scale_y_continuous(name = "Position au sein du Pseudotemps Monocle3", limits = c(0, 3750), breaks = seq(0,3800,500)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(0.5, 'cm')) +
  labs(title = "") +
  annotate("text", x = 300, y = 3700, label = stats[1], fontface = 2, colour = "black") +
  annotate("text", x = 500, y = 3550, label = stats[2], fontface = 2, colour = "black")
plot9



################################################################################
### LEAP MATRIX PREPARATION ----------------------------------------------------
################################################################################

setwd(paste(PATH, "/Saves", sep=""))
data <- readRDS(paste(PATH, "/Saves/Integration_Pseudotime", sep = ""))

# Generation of leap_matrix.csv (LEAP input) and gene_index.csv (for gene annotation after LEAP)
matrix_leap <- LeapMatrix(data, pseudotime = 'Slingshot', min.cells = 127, write = FALSE) 
# /!\ The .csv file written does not contain row names contrary to the variable matrix_leap



################################################################################
### VALIDATION WITH KNOWN GENES ------------------------------------------------
################################################################################
#Based on : Lineage specification in innate lymphocytes, Das et al. (Elsevier) 2018

gene_list <- c('Il7r', 'Nfil3', 'Tcf7', 'Gata3')

ordered <- OrderMatrix(data, order.by = 'Slingshot', features = gene_list)
ordered <- PseudotimeRepartition(ordered, 21, order.by = 'Slingshot')

plot10 <- DrawExpr(ordered, gene_list[1], std = TRUE, col = ColorBlind12[1], xlab = "Pseudotemps", ylab = "Expression moyenne")
plot10 <- recordPlot()
plot11 <- DrawExpr(ordered, gene_list[2], std = TRUE, col = ColorBlind12[1], xlab = "Pseudotemps", ylab = "Expression moyenne")
plot11 <- recordPlot()
plot12 <- DrawExpr(ordered, gene_list[3], std = TRUE, col = ColorBlind12[1])
plot12 <- recordPlot()
plot13 <- DrawExpr(ordered, gene_list[4], std = TRUE, col = ColorBlind12[1])
plot13 <- recordPlot()


## LINEAR REPARTITION OF CELL TYPE ACROSS PSEUDOTIME

ggplot(compare, aes(y = 1, x = Slingshot_order*3)) +
  geom_point(aes(colour = factor(Sample, levels = c("ALP", "overLIP", "Tulip", "ILCpro"))), size = 5) +
  scale_color_manual(values = ColorBlind12[c(6, 5, 1, 4)], name = "") + 
  scale_x_continuous(name = "Position au sein du Pseudotemps Slinghshot", limits = c(0, 3800*3)) +
  scale_y_continuous(name = "Position au sein du Pseudotemps Monocle3", limits = c(0.5, 1.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(0.5, 'cm')) +
  annotate("text", x = 300, y = 3700, label = stats[1], fontface = 2, colour = "black") +
  annotate("text", x = 500, y = 3550, label = stats[2], fontface = 2, colour = "black") +
  labs(title = "") 

SamplePos <- list() 
for(s in compare$Sample){
  SamplePos[[s]] <- compare$Slingshot_order[compare$Sample == s]
}


Posistats <- data.frame()
for(s in names(table(compare$Sample))){
  mean <- mean(compare$Slingshot_order[compare$Sample == s])
  std <- sd(compare$Slingshot_order[compare$Sample == s])
  Posistats <- rbind(Posistats,c(s, mean, std))
}
colnames(Posistats) <- c("sample", "mean", "std")


ggplot() +
  geom_point(data = compare, aes(y = factor(Sample, levels = rev(c("ALP", "overLIP", "Tulip", "ILCpro"))), x = Slingshot_order, colour = factor(Sample, levels = c("ALP", "overLIP", "Tulip", "ILCpro"))), size = 1) +
  geom_point(data = Posistats, aes(y = factor(sample, levels = rev(c("ALP", "overLIP", "Tulip", "ILCpro"))),  x = round(as.numeric(mean)))) +
  geom_errorbarh(data = Posistats, aes(y = factor(sample, levels = rev(c("ALP", "overLIP", "Tulip", "ILCpro"))), xmax = as.numeric(mean) + as.numeric(std), xmin = as.numeric(mean) - as.numeric(std))) +
  scale_color_manual(values = ColorBlind12[c(6, 5, 1, 4)], name = "") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "", x = "", y ="") +
  NoLegend()


ggplot(compare, aes(x=Slingshot_order, fill=factor(Sample, levels = c("ALP", "overLIP", "TULIP", "ILCpro")))) +
  geom_density(alpha=0.4) +
  theme_classic() +
  scale_fill_manual(values = ColorBlind12[c(6, 5, 1, 4)], name = "")
  
  













