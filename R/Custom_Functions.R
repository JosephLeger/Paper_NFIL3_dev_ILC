#!/usr/bin/env Rscript

#===============================================================================
# DESCRIPTION ------------------------------------------------------------------
#===============================================================================
# File containing useful custom functions employed for DNase/RNA-seq analyzes



#===============================================================================
# BASICS -----------------------------------------------------------------------
#===============================================================================

`%!in%` <- Negate(`%in%`)

# Colorblind friendly pack
Blues <- c('#016BA0', '#5F9ED1', '#A1C8EB', '#3399FF',  '#006ddb', '#330066')
YReds <- c('#FF0000', '#FF9933', '#F0E442', '#ffff6d', '#FF7070')
Greys <- c('#131313','#555555', '#809099', '#CFCFCF', '#F2F4F4')
Others <- c('#FFBBAB')

ColorBlind <- c(Blues, Others, YReds, Greys)
#pie(rep(1,17), col=ColorBlind)



#===============================================================================
# FILE WRITING -----------------------------------------------------------------
#===============================================================================

writePlot <- function(plot, path, filename=F, width=9, height=9, res=100){
  
  # Save plot as PNG formated file named by the variable containing it
  #
  # plot     = variable containing plot to save
  # path     = path of the directory to save the plot
  # filename = could be used to define filename other than the variable name
  # width    = width of saved plot (inches)
  # height   = height of saved plot (inches)
  # res      = resolution to save
  
  if(filename == F){
    file <- paste0(path, '/', deparse(substitute(plot)), '.png')
  }else{
    file <- paste0(path, '/', filename, '.png')
  }
  show(plot)
  dev.print(png, file, width=width, height=height, units='in', res=res)
}



#===============================================================================
# BULK RNA-SEQ FUNCTIONS -------------------------------------------------------
#===============================================================================

Quantile_Normalization <- function(df){
  
  # Perform quantile-normalization on given dataframe
  # from : https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
  #
  # df = dataframe object to normalize 
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_rank            <- apply(df, 2, rank, ties.method='min')
  df_sorted          <- data.frame(apply(df, 2, sort))
  df_mean            <- apply(df_sorted, 1, mean)
  df_final           <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  
  return(as.data.frame(df_final))
}

#-------------------------------------------------------------------------------

DrawHeatmap <- function(counts, sampleSheet, genes, 
                        by.group = F, 
                        groups = NA,
                        logT = F,
                        scale = 'row',
                        cluster_cols = F, 
                        cluster_rows = T,
                        treeheight_row = 50,
                        treeheight_col = 50,
                        color = colorRampPalette(c("blue","white","red"))(50),
                        border_color = NA,
                        gaps_col = 0,
                        gaps_row = 0,
                        fontsize = 10,
                        fontsize_row = fontsize,
                        fontsize_col = fontsize,
                        title = ''){
  
  # Draw a heatmap of gene expression from a table
  #
  # counts         = Table with gene expressions and names as Symbol column
  # sampleSheet    = metadata containing a column called 'Group'
  # genes          = gene names to draw
  # by.group       = whether draw heatmap for average expression by Group
  # groups         = groups to consider from sampleSheet$Group column
  # logT           = whether log transform expression values
  # scale          = see pheatmap parameters
  # cluster_cols   = see pheatmap parameters
  # cluster_rows   = see pheatmap parameters
  # treeheight_row = see pheatmap parameters
  # treeheight_col = see pheatmap parameters
  # color          = see pheatmap parameters
  # border_color   = see pheatmap parameters
  # gaps_col       = see pheatmap parameters
  # gaps_row       = see pheatmap parameters
  # fontsize       = see pheatmap parameters
  # fontsize_row   = see pheatmap parameters
  # fontsize_col   = see pheatmap parameters
  # title          = main title of the generated plot
  
  suppressPackageStartupMessages(library(pheatmap))
  
  # Abort function if 'Group' column not found in sampleSheet
  if('Group' %!in% colnames(sampleSheet)){
    stop(paste0('Provided sampleSheet does not contain \'Group\' column. ',
                'Please rename the column to consider to group samples.'))
  }
  
  # Subset table - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Subset sampleSheet based on defined or default groups
  if(length(groups) == 1 && is.na(groups)){
    groups        <- unique(sampleSheet$Group)
  }
  sampleSheet     <- sampleSheet[sampleSheet$Group %in% groups,]
  
  # Subset counts based on groups and expressed genes
  y               <- counts[counts$Symbol %in% genes, 
                            c('Symbol', sampleSheet$Sample)]
  
  # Remove zero expressed genes in the corresponding subset of samples
  y               <- y[rowSums(y[,colnames(y) %!in% 'Symbol']) != 0,]
  
  # Calculate mean expression if required or just filter count table
  if(by.group){
    y.mean        <- list(Symbol = y$Symbol)
    for(g in groups){
      samples     <- sampleSheet$Sample[sampleSheet$Group %in% g]
      y.sub       <- y[,colnames(y) %in% samples]
      y.mean[[g]] <- rowSums(y.sub)/ncol(y.sub)
    }
    y <- data.frame(y.mean)
    
  }else{
    samples       <- c()
    for(g in groups){
      samples     <- c(samples, sampleSheet$Sample[sampleSheet$Group %in% g])
    }  
    y             <- y[,c('Symbol', samples)]
  }
  
  # Manage logT if required
  if(logT){
    y.log <- log(y[,colnames(y) %!in% 'Symbol']+1)
    y     <- cbind(Symbol = y$Symbol, y.log)
    scale <- 'none'
  }
  
  # Draw Heatmap - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  plot <- pheatmap(y[,colnames(y) %!in% 'Symbol'],
                   cluster_cols = cluster_cols, 
                   cluster_rows = cluster_rows, 
                   treeheight_row = treeheight_row,
                   treeheight_col = treeheight_col,
                   scale = scale,
                   color = color,
                   border_color = border_color,
                   gaps_col = gaps_col,
                   gaps_row = gaps_row,
                   labels_row = y$Symbol, 
                   fontsize = fontsize,
                   fontsize_row = fontsize_row,
                   fontsize_col = fontsize_col,
                   main = title)
}

#-------------------------------------------------------------------------------

DrawDotplot <- function(counts, sampleSheet, gene, 
                        groups=NA, 
                        repeated=NA,
                        colstat=NA,
                        show.fc = TRUE,
                        main = gene){
  
  # Draw dotplot figure for a single gene using ggplot2
  #
  # counts      = Worked_Table with gene expressions and names as Symbol column
  # sampleSheet = metadata
  # gene        = gene name to draw
  # groups      = groups to consider from sampleSheet$Group column
  # repeated    = whether samples are linked and add a line for sample evolution
  # colstat     = define stat column from sampleSheet to use for plot title *
  # show.fc     = whether color the line for sample evolution
  # main        = Plot title
  
  suppressPackageStartupMessages(library(ggplot2))
  
  # Abort function if 'Group' column not found in sampleSheet
  if('Group' %!in% colnames(sampleSheet)){
    stop(paste0('Provided sampleSheet does not contain \'Group\' column. ',
                'Please rename the column to consider to group samples.'))
  }
  
  # Subset counts based on expressed genes and groups - - - - - - - - - - - - -
  if(length(groups) == 1 && is.na(groups)){
    groups    <- unique(sampleSheet$Group)
  }
  sampleSheet <- sampleSheet[sampleSheet$Group %in% groups,]
  
  # Extract information about provided gene
  geneline <- counts[counts$Symbol %in% gene, sampleSheet$Sample]
  
  if(nrow(geneline) > 1){
    warning(paste0('Multiple rows found for ', gene, 
                   '. Only first line will be considered.'))
  }else if(nrow(geneline) == 0){
    stop(paste(gene, 'gene not found.'))
  }
  
  # Organize matrix for ggplot2 with gene expression and corresponding group
  mat    <- data.frame(t(geneline))
  g_id   <- c()
  l_id   <- c()
  for(r in rownames(mat)){
    g_id <- c(g_id, sampleSheet$Group[sampleSheet$Sample %in% r])
    
    # Optionally add repeated parameter
    if(!is.na(repeated)){
      l_id <- c(l_id, sampleSheet[,repeated][sampleSheet$Sample %in% r])
    }else{
      l_id <- c(1:nrow(mat))
    }
  }
  mat           <- cbind(mat, g_id, l_id)
  colnames(mat) <- c('Expression', 'Group', 'Repeat')
  
  # Calculate FC by sample if needed
  if(!is.na(repeated) && show.fc){
    mat$log2FC <- rep(NA, nrow(mat))
    for(s in sampleSheet$Sample[sampleSheet$Group %in% groups[1]]){
      repnum <- sampleSheet[,repeated][sampleSheet$Sample %in% s]
      s2     <- sampleSheet$Sample[sampleSheet$Group %in% groups[2] &
                                     sampleSheet[,repeated] %in% repnum]
      mat[s,'log2FC']  <- log2(mat[s2,1]/mat[s,1])
      mat[s2,'log2FC'] <- log2(mat[s2,1]/mat[s,1])
    }
  }else{
    mat$log2FC <- runif(nrow(mat),0,1)
  }
  
  # Define plot title adding or not statistical relevance
  if(!is.na(colstat)){
    stat <- counts[,colstat][counts$Symbol %in% gene]
    if(as.numeric(stat) < 0.001){
      title <- paste(main, '****')
    }else if(as.numeric(stat) < 0.005){
      title <- paste(main, '***')
    }else if(as.numeric(stat) < 0.01){
      title <- paste(main, '**')
    }else if(as.numeric(stat) < 0.05){
      title <- paste(main, '*')
    }else{
      title <- main
    }
  }else{
    title <- main
  }
  
  # Draw plot depending on selected options
  if(!is.na(repeated) && show.fc){
    plot <- ggplot(mat, aes(x=factor(Group, levels=groups), y=Expression)) +
      geom_line(aes(group=Repeat, color=log2FC), alpha=1, linewidth=1.2, 
                show.legend = TRUE) +
      scale_color_gradient2(low='blue', mid = 'beige', high='red') +
      geom_point(size = 5)
  }else if(!is.na(repeated)){
    plot <- ggplot(mat, aes(x=factor(Group, levels=groups), y=Expression)) +
      geom_line(aes(group=Repeat, color=log2FC), alpha=1, linewidth=1.2, 
                show.legend = FALSE) +
      scale_color_gradient(low='grey', high='grey') +
      geom_point(size = 5)
    print('choice2')
  }else{
    plot <- ggplot(mat, aes(x=factor(Group, levels=groups), y=Expression)) +
      geom_point(size = 5)
  }
  plot <- plot + 
    stat_summary(fun.data=mean_sdl, fun.args=list(mult=1), geom='errorbar', 
                 color='black', width=0.1, linewidth=0.5) +
    stat_summary(fun='mean', geom='point', color='black', group=1, shape=3, 
                 size=10) +
    theme_classic() +
    ggtitle(title)
  
  return(plot)
}


  
#===============================================================================
# LEAP -------------------------------------------------------------------------
#===============================================================================

OrderMatrix <- function(x, trajectory, min.cells = 0, assay = x@active.assay, 
                        slot = "data"){
  
  # Return an ordered matrix by pseudotime for chosen trajectory
  #
  # x          = Seurat object with pseudotime in metadata
  # trajectory = name of metadata containing selected pseudotime
  # min.cells  = filters genes expressed in at least min.number cells
  # assay      = assay to consider for analysis
  # slot       = name of slot to use (data, scaledata or counts)
  
  
  # V2 - Get the matrix corresponding to slot chosen
  if(slot == "data"){
    tablex <- as.data.frame(t(as.matrix(x@assays[[x@active.assay]]@data)))
  }else if(slot == "scaledata"){
    tablex <- as.data.frame(t(as.matrix(x@assays[[x@active.assay]]@scale.data)))
  }else if(slot == "counts"){
    tablex <- as.data.frame(t(as.matrix(x@assays[[x@active.assay]]@counts)))
  }else{
    stop("slot must be in 'data', 'scaledata' or 'counts'")
  }
  
  
  pseudotime <- x@meta.data[[trajectory]]
  matrix_leap <- cbind(pseudotime, tablex)
  colnames(matrix_leap) <- c('Pseudotime', colnames(tablex))
  
  # Keeping only cells belonging to chosen trajectory and sort them
  matrix_leap  <- subset(matrix_leap, subset = !is.na(matrix_leap$Pseudotime))
  matrix_leap  <- matrix_leap[order(matrix_leap$Pseudotime),]
  
  
  # V2 - Eliminate genes expressed in less than min.cells
  cell_num <- colSums(matrix_leap != 0)
  features_to_keep <- unique(c('Pseudotime', names(cell_num[cell_num >= 
                                                              min.cells])))
  matrix_leap <- matrix_leap[,features_to_keep] 
  
  # V2 -  Warning if datascale slot is used
  if(min.cells != 0 & slot == 'scaledata'){
    warning('min.cells not efficient for scaledata slot')
  }
  
  return(matrix_leap)
}

#-------------------------------------------------------------------------------

CalculateBin <- function(x, bin.number){
  
  # Establish bin cutoffs for PseudotimeRepartition() and representation 
  #
  # x          = ordered matrix
  # bin.number = number of bin for representation
  
  
  p <- as.vector(x$Pseudotime)
  # Establish the min and max value of pseudotime to generate the bin cutoff
  min.val    <- min(as.numeric(p))
  max.val    <- max(as.numeric(p))
  # calculate the range of a single bin
  bin.range  <- ((max.val - min.val)/bin.number)
  result     <- c(min.val)
  for(i in 1:(bin.number-1)){
    previous <- result[length(result)]
    result   <- c(result, (previous + bin.range))
  }
  # Manually add max.val to ensure issues due to decimale
  result     <- c(result, max.val)
  return(result)
}

#-------------------------------------------------------------------------------

PseudotimeRepartition <- function(x, bin.number = 21, bin.level = NULL){
  
  # Attribute a bin value to each cell and return the ordered matrix
  #
  # x          = ordered matrix 
  # bin.number = number of bin for representation
  # bin.level  = result of CalculateBin()
  
  
  # Calculate bin.level if not already done and provided
  if(is.null(bin.level)){
    bin.level <- CalculateBin(x, bin.number)
  }
  
  bin.centers <- c()
  bin.range   <- max(x$Pseudotime)/bin.number
  for(b in 1:(length(bin.level)-1)){
    bin.centers <- c(bin.centers, bin.level[b]+(bin.range/2))
  } 
  
  i.bin       <- 1
  attribution <- c()
  
  # Attribute a bin to each cell
  for(i in 1:nrow(x)){
    if(as.numeric(x$Pseudotime[i]) <= bin.level[i.bin+1]){
      attribution <- c(attribution, bin.centers[i.bin])
    } else {
      attribution <- c(attribution, bin.centers[i.bin+1])
      i.bin <- i.bin+1
    }
  }
  x <- cbind(attribution, x)
  colnames(x) <- c("Bin", colnames(x)[2:ncol(x)])
  
  return(x)
}

#-------------------------------------------------------------------------------

LeapMatrix <- function(x, trajectory, min.cells = 0, assay = x@active.assay, 
                       slot = 'data', write = FALSE, dir = getwd(), 
                       suffix = NULL){
  
  # Generates matrix_leap and gene_index for LEAP analysis
  #
  # x = Seurat object with pseudotime is metadata
  # trajectory = name of metadata containing selected pseudotime
  # min.cells  = filters genes expressed in at least min.number cells
  # assay      = assay to consider for analysis
  # slot       = default uses @data, but could be @scaled 
  # write      = either save matrix_leap.txt and gene_index.txt files
  # dir        = directory in which save files
  # suffix     = suffix to use to associate and recognize files
  
  
  # If suffix is not defined, uses first and last character of trajectory
  if(is.null(suffix)){
    suffix <- paste0(substr(trajectory, 1, 1), 
                     substr(trajectory, nchar(trajectory),nchar(trajectory)))
  }
  
  # Run OrderMatrix and eliminates pseudotime row to get LEAP matrix format
  leap_matrix <- OrderMatrix(x, trajectory, min.cells, assay, slot)
  leap_matrix <- leap_matrix[,colnames(leap_matrix) %!in% "Pseudotime"]
  leap_matrix <- as.data.frame(t(leap_matrix))
  
  # Generates correspong gene_index
  gene_index <- data.frame(id = c(1:nrow(leap_matrix)), 
                           gene = rownames(leap_matrix))
  
  if(write){
    write.table(gene_index, paste0(dir, '/gene_index_', suffix, '.txt'), 
                row.names = FALSE)
    write.table(leap_matrix, paste0(dir, '/matrix_leap_', suffix, '.txt'), 
                row.names = FALSE)
  }
  
  result <- list(leap_matrix, gene_index)
  return(result)
}

#-------------------------------------------------------------------------------

AnnotateLeap <- function(MAC, index, write = FALSE, dir = getwd(), 
                         filename = 'results_indexed'){
  
  # Annotate LEAP output with gene names from an index previously generated
  #
  # MAC      = LEAP output of MAC_counter function (results.txt)
  # index    = name of variable which contains gene names and corresponding id
  # write    = whether write a results_indexed.csv file in working directory
  # dir      = directory in which save files
  # filename = output filename to save
  
  
  results    <- MAC
  gene_index <- index
  # Replacing gene id with gene names
  colnames(results)          <- c('Correlation','Lag','Row gene index', 
                                  'Column gene index')
  colnames(gene_index)       <- c('Row gene index', 'gene_row')
  merged                     <- merge(results, gene_index, 
                                      by = 'Row gene index')
  colnames(gene_index)       <- c('Column gene index', 'gene_col')
  merged                     <- merge(merged, gene_index, 
                                      by = 'Column gene index')
  merged$`Column gene index` <- NULL
  merged$`Row gene index`    <- NULL
  if(write == TRUE){
    write.table(merged, paste0(dir, '/', filename, '.txt'), row.names = FALSE, 
                col.names =T, quote = F)
  }
  return(merged)
}



#===============================================================================
# PSEUDOTIME PLOTS -------------------------------------------------------------
#===============================================================================

DrawExpr <- function(x, bin.number = 21, feature.list, 
                     by.order = F, std = F, scale = F, superposed = F, 
                     centroid = F, color = ColorBlind, lwd = 2, 
                     ylim = c(-0.75, 2), xlab = 'Pseudotime', 
                     ylab = 'Expression', main = NA){
  
  # Draw single or multiple curves to show expression from ordered matrix
  #
  # x            = ordered matrix 
  # bin.number   = number of bin for representation
  # feature.list = features to draw
  # by.order     = if TRUE, consider cell order and not pseudotime to make bins
  # std          = whether draw std in plot (only for superposed = FALSE)
  # scale        = whether scale shown expression levels from 0 to 1
  # superposed   = whether draw all features superposed in a single plot
  # centroid     = whether draw centroid curve (superposed must be True)
  # color        = colors for plots (several can be specified for superposed)
  # lwd          = line thickness
  # ylim         = ylim for plots (can be manually adjusted)
  # xlab         = label on x axis for plots (default = "Pseudotime")
  # ylab         = label on y axis for plots (default = "Expression")
  # main         = specified title for superposed plot
  
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  
  # If x axis is setted to be cell ordered by rank and not by pseudotime value
  if(by.order){
    x$Pseudotime <- 1:nrow(x)
  }
  
  # Select only existing genes
  absent_genes <- feature.list[feature.list %!in% colnames(x)]
  feature.list <- feature.list[feature.list %in% colnames(x)]
  if(length(absent_genes) != 0){
    warning(paste0('Following genes not found in provided matrix :', 
                   paste(absent_genes, collapse = ', ')))
  }
  
  # Subset table for more efficiency during following steps and order it
  x <- x[,colnames(x) %in% c('Pseudotime', feature.list)]
  x <- x[,order(match(colnames(x), c('Pseudotime', feature.list))), 
         drop = FALSE]
  
  # Attributes Bin to each cell
  x <- PseudotimeRepartition(x, bin.number)
  
  # Initiate loop
  x.dim         <- as.numeric(unique(x$Bin))
  expr_table    <- data.frame()

  # Generate expr_table containing gene, bin, mean expr, std - - - - - - - - - -
  for(feature in feature.list){
    for(b in x.dim){
      # Subset Bin
      values     <- x[x$Bin %in% b, feature]
      expr_table <- rbind(expr_table, c(feature, b, mean(values), sd(values)))
    }
    colnames(expr_table) <- c('Gene', 'Bin', 'Expr', 'Std')
    expr_table$Gene <- factor(expr_table$Gene, levels = feature.list)
    
    # Scale current gene mean expression if required - - - - - - - - - - - - -
    if(scale){
      scaled_expr <- as.numeric(expr_table$Expr[expr_table$Gene %in% feature])
      scaled_expr <- (scaled_expr-min(scaled_expr))/
        (max(scaled_expr)-min(scaled_expr))
      expr_table$Expr[expr_table$Gene %in% feature] <- scaled_expr
      # Modify parameters to correspond correctly to scaled plot
      ylim        <- c(0,1)
      std         <- F
    }
    
    # Draw individual curves if not superposed - - - - - - - - - - - - - - - -
    if(!superposed){
      p <- ggplot(expr_table[expr_table$Gene %in% feature,], 
                  aes(x = as.numeric(Bin), group = Gene)) +
        geom_line(linewidth = lwd, aes(y = as.numeric(Expr), color = Gene)) +
        scale_color_manual(values = color) +
        theme_classic() +
        ylim(ylim) 
      if(std){
        p <- p + geom_ribbon(
          aes(y = as.numeric(Expr), ymin = as.numeric(Expr) - as.numeric(Std), 
              ymax = as.numeric(Expr) + as.numeric(Std), fill = Gene), 
          alpha = 0.1,) + scale_fill_manual(values = color)
      }
      p <- p + labs(x=xlab, y=ylab, title=ifelse(is.na(main), feature, main))
      print(p)
    }
  }
  
  # Draw merged plot if required - - - - - - - - - - - - - - - - - - - - - - - - 
  if(superposed){
    p <- ggplot(expr_table, 
                aes(x = as.numeric(Bin), group = Gene)) +
      geom_line(linewidth = lwd, aes(y = as.numeric(Expr), color = Gene)) +
      scale_color_manual(values = color) +
      theme_classic() +
      ylim(ylim) 
    if(std){
      p <- p + geom_ribbon(aes(y = as.numeric(Expr), 
                               ymin = as.numeric(Expr) - as.numeric(Std), 
                               ymax = as.numeric(Expr) + as.numeric(Std),
                               fill = Gene), alpha = 0.1,) +
        scale_fill_manual(values = color)
    }
    if(centroid){
      centroid <- data.frame()
      for(b in x.dim){
        centroid <- rbind(
          centroid, c('centroid', b, 
                      mean(as.numeric(expr_table$Expr[expr_table$Bin %in% b]))))
      }
      colnames(centroid) <- c('Gene', 'Bin', 'Expr')
      p <- p + geom_line(data = centroid, linewidth = lwd,
                         aes(x = as.numeric(Bin), y = as.numeric(Expr)))
    }
    p <- p + labs(x=xlab, y=ylab, title=ifelse(is.na(main), '', main))
    print(p)
  }
  return(expr_table)
}

#-------------------------------------------------------------------------------

CompareExpr <- function(x, bin.number = 11, feature, 
                        by.order = F, scale = T, min.cells = 5, stat = 't', 
                        color = c("royalblue","red3"), lwd = 2, ylim = NULL, 
                        xlab = 'Pseudotime', ylab = 'Expression', 
                        main = feature, show = T){
  
  # Compare feature expression between samples along one pseudotime trajectory
  #
  # x          = ordered matrix 
  # bin.number = number of bin for representation
  # feature    = feature to draw
  # by.order   = if TRUE, consider cell order and not pseudotime to make bins
  # scale      = whether scale shown expression levels from 0 to 1
  # min.cells  = minimal number of cells to consider valid bin to plot mean
  # stat       = statistical test used, could be w (Wilcox) or t (Student T)
  # color      = colors for plots (several can be specified for superposed)
  # lwd        = line thickness
  # ylim       = ylim for plots (can be manually adjusted)
  # xlab       = label on x axis for plots (default = "Pseudotime")
  # ylab       = label on y axis for plots (default = "Expression")
  # main       = specified title for superposed plot
  # show       = whether show the produced plot (use false to avoid memory use)
  
  suppressPackageStartupMessages(library(ggplot2))
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  get_significance <- function(p) {
    if(is.na(p)){return(NA)
    }else if(p < 0.0001){return('****')
    }else if(p < 0.001){return('***')
    }else if(p < 0.01){return('**')
    }else if(p < 0.05){return('*')
    }else{return(' ')}
  }
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  # If feature is not present in the table, stop immediatly the process
  if(feature %!in% colnames(x)){
    stop(paste(feature, 'not found as a feature column in the provided matrix'))
  }
  
  # If x axis is set by cell rank and not by pseudotime value
  if (by.order) {
    x           <- x[rownames(x) %!in% 'Pseudotime',]
    x           <- rbind(c(1:ncol(x)), x)
    rownames(x) <- c('Pseudotime', rownames(x)[2:nrow(x)])
  }
  
  # Attribute Bin to each cell
  x             <- PseudotimeRepartition(x, bin.number)
  
  # Initiate loop with result_table and a list to store values by samples
  result_table  <- data.frame()
  expr_list     <- list()
  
  # Loop over each unique bin value
  for(b in unique(x$Bin)){
    bin_means <- c()
    bin_stats <- c()
    
    # Loop over each sample to get values for current bin and calculate means
    for(s in unique(x$Sample)){
      expr_list[[s]] <- x[feature][x$Bin %in% b & x$Sample %in% s,]
      # Only consider groups above a defined cellnumber threshold (5 cells)
      if(length(expr_list[[s]]) >= min.cells){
        bin_means    <- c(bin_means, mean(expr_list[[s]]))
      }else{
        bin_means    <- c(bin_means, NA)
      }
    }
    names(bin_means) <- unique(x$Sample)
    
    # Loop over all possible comparisions to perform statistical test
    for(i in 1:(length(expr_list)-1)){
      for(j in (i+1):length(expr_list)){
        
        # Only consider groups above a defined cellnumber threshold (5 cells)
        if(length(expr_list[[i]]) >= min.cells & 
           length(expr_list[[j]]) >= min.cells) {
          if (stat %in% 'w') {
            test_result <- wilcox.test(expr_list[[i]], expr_list[[j]])
          }else if(stat %in% 't'){
            test_result <- t.test(expr_list[[i]], expr_list[[j]])
          }
          bin_stats     <- c(bin_stats, test_result$p.value)
        }else{
          bin_stats     <- c(bin_stats, NA)
        }
        names(bin_stats)[length(bin_stats)] <- paste0(
          names(expr_list)[i], '_vs_', names(expr_list)[j])
      }
    }
    result_table           <- rbind(result_table, c(b, bin_means, bin_stats))
    colnames(result_table) <- c('Bin', names(bin_means), names(bin_stats))
  }
  
  # Replace mean expression by scaled ones if required
  if(scale){
    scaled <- result_table[,2:(length(expr_list)+1)]
    scaled <- (scaled-min(scaled, na.rm = T))/
      (max(scaled, na.rm = T)-min(scaled, na.rm = T))
    result_table[,2:(length(expr_list)+1)] <- scaled
  }
  
  # Define ylim if not already done
  if(length(ylim) %in% 0){
    ylim   <- c(min(result_table[,2:(length(expr_list)+1)], na.rm = T),
                max(result_table[,2:(length(expr_list)+1)], na.rm = T))
  }
  
  # Plot results - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  gg_results <- data.frame(
    Bin        = rep(result_table$Bin, length(expr_list)),
    Expression = stack(result_table[,2:(length(expr_list)+1)])[,1],
    Group      = stack(result_table[,2:(length(expr_list)+1)])[,2])
  
  p <- ggplot(gg_results, aes(x = Bin, y = Expression, color = Group)) +
    geom_line(linewidth = lwd) +
    scale_color_manual(values = color) + 
    theme_classic() +
    ylim(ylim) + labs(x = xlab, y = ylab, title = main)
  
  if(show){print(p)}
  
  return(result_table)
}


  
#===============================================================================
# ANNOTATION -------------------------------------------------------------------
#===============================================================================

EnsemblToGeneSymbol <- function(x, refindex, duplicate = T, sum.duplicate = T){
  
  # Convert rownames from Ensembl to Gene Symbol classification
  #
  # x             = count table with Ensembl rownames
  # refindex      = downloaded list for conversion (need to be updated)
  # duplicate     = whether duplicate rows sharing Gene Symbol (gene.1, gene.2)
  # sum.duplicate = whether add a row gene.ALL summing duplicated genes
  
  
  # Add GeneSymbol annotated information
  y <- cbind(Gene.stable.ID = rownames(x), x)
  y <- merge(refindex, y, by = "Gene.stable.ID", all.y = T)
  
  # Complete missing GeneSymbol
  missing <- which(y$Gene.name == "" | is.na(y$Gene.name))
  for(i in missing){
    y[i, "Gene.name"] <- y[i, "Gene.stable.ID"]
  }
  
  # Deal with non-unique gene names
  non_unique_num <- which(!isUnique(y$Gene.name))
  non_unique_names <- unique(y$Gene.name[!isUnique(y$Gene.name)])
  
  
  if(!duplicate){
    for(i in non_unique_num){
      y[i, "Gene.name"] <- y[i, "Gene.stable.ID"]
    }
    rownames(y) <- y$Gene.name
    y <- y[, colnames(y) %!in% c("Gene.stable.ID", "Gene.stable.ID.version", 
                                 "Gene.name")]
  } else {
    y$Gene.name <- make.names(y$Gene.name, unique = TRUE)
    
    if(sum.duplicate){
      z   <- data.frame()
      for(g in non_unique_names){
        pattern          <- paste("^", g, "|", g, "/.[123456789]", sep = "")
        sub_table        <- y[grepl(pattern, y$Gene.name),]
        sub_table        <- sub_table[, colnames(sub_table) %!in% c(
          "Gene.stable.ID", "Gene.stable.ID.version", "Gene.name")]
        new_line <- c(paste(g, ".ALL", sep = ""), apply(sub_table, 2, sum))
        z <- rbind(z, new_line)
      }
      rownames(z) <- z[,1]
      z <- z[,2:ncol(z)]
      colnames(z) <- colnames(y)[colnames(y) %!in% c(
        "Gene.stable.ID", "Gene.stable.ID.version", "Gene.name")]
    }
    rownames(y) <- y$Gene.name
    y <- y[, colnames(y) %!in% c(
      "Gene.stable.ID", "Gene.stable.ID.version", "Gene.name")]
    y <- rbind(y, z)
  }
  return(y)
}  

#-------------------------------------------------------------------------------
