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
ColorBlind12 <- c('#016BA0', '#CA5100', '#5F9ED1', '#FF9933', '#A1C8EB', 
                  '#F0E442', '#FFBBAB', '#555555', '#809099', '#CFCFCF', 
                  '#3399FF', '#330066')

# Custom 
Blues <- c('#016BA0', '#5F9ED1', '#A1C8EB', '#3399FF',  '#006ddb', '#330066' )
YReds <- c('#FF0000', '#FF9933', '#F0E442', '#ffff6d', "#FF7070")
Greys <- c('#131313','#555555', '#809099', '#CFCFCF', '#F2F4F4')
Others <- c('#FFBBAB')

ColorBlind <- c(Blues, Others, YReds, Greys)
#pie(rep(1,17), col=ColorBlind)



#===============================================================================
# PLOTS ------------------------------------------------------------------------
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



#===============================================================================
# LEAP -------------------------------------------------------------------------
#===============================================================================

OrderMatrix <- function(x, trajectory, min.cells = 0, assay = x@active.assay, slot = "data"){
  
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
  features_to_keep <- unique(c('Pseudotime', names(cell_num[cell_num >= min.cells])))
  matrix_leap <- matrix_leap[,features_to_keep] 
  
  # V2 -  Warning if datascale slot is used
  if(min.cells != 0 & slot == 'scaledata'){
    warning('min.cells not efficient for scaledata slot')
  }
  
  return(matrix_leap)
}

#------------------------------------------------------------------------------#

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

#------------------------------------------------------------------------------#

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

#------------------------------------------------------------------------------#

LeapMatrix <- function(x, trajectory, min.cells = 0, assay = x@active.assay, slot = "data",
                       write = FALSE, dir = getwd(), suffix = NULL){
  
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
    suffix <- paste(substr(trajectory, 1, 1), substr(trajectory, nchar(trajectory),nchar(trajectory)), sep = "")
  }
  
  # Run OrderMatrix and eliminates pseudotime row to get LEAP matrix format
  leap_matrix <- OrderMatrix(x, trajectory, min.cells, assay, slot)
  leap_matrix <- leap_matrix[,colnames(leap_matrix) %!in% "Pseudotime"]
  leap_matrix <- as.data.frame(t(leap_matrix))
  
  # Generates correspong gene_index
  gene_index <- data.frame(id = c(1:nrow(leap_matrix)), gene = rownames(leap_matrix))
  
  if(write){
    write.table(gene_index, paste(dir, "/gene_index_", suffix, ".txt", sep = ""), row.names = FALSE)
    write.table(leap_matrix, paste(dir, "/matrix_leap_", suffix, ".txt", sep = ""), row.names = FALSE)
  }
  
  result <- list(leap_matrix, gene_index)
  return(result)
}

#------------------------------------------------------------------------------#

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
  colnames(results)          <- c('Correlation','Lag','Row gene index', 'Column gene index')
  colnames(gene_index)       <- c('Row gene index', 'gene_row')
  merged                     <- merge(results, gene_index, by = 'Row gene index')
  colnames(gene_index)       <- c('Column gene index', 'gene_col')
  merged                     <- merge(merged, gene_index, by = 'Column gene index')
  merged$`Column gene index` <- NULL
  merged$`Row gene index`    <- NULL
  if(write == TRUE){
    write.table(merged, paste0(dir, '/', filename, '.txt'), row.names = FALSE, col.names =T, quote = F)
  }
  return(merged)
}

#------------------------------------------------------------------------------#

DrawExpr <- function(x, bin.number = 21, feature.list, 
                     by.order = F, std = F, scale = F, scale.method = 0,
                     superposed = F, write = F, dir = getwd(), 
                     col = "black", lwd = 2, ylim = c(-0.75, 2), 
                     xlab = "Pseudotime", ylab = "", main = "Untitled", 
                     width = 900, height = 800){
  
  # Faster way to plot multiple expression curves or superposed curves
  #
  # x            = ordered matrix 
  # bin.number   = number of bin for representation
  # feature_list = features to draw
  # by.order     = if set to TRUE, consider cell order rather than pseudotimeto make bins
  # std          = whether draw std in plot (only for superposed = FALSE)
  # scale        = whether scale shwon expression levels from 0 to 1
  # scale.method = 0 : absolute [0,1] / 1 : MinMax / 2 : amplitude
  # superposed   = whether draw all features superposed in a single plot
  # write        = whether save resulting plots
  # dir          = output directory to save plots
  # col          = colors for plots (several can be specified for superposed plot)
  # ylim         = ylim for plots (can be manually adjusted)
  # xlab         = label on x axis for plots (default = "Pseudoime")
  # ylab         = label on y axis for plots (default = "")
  # main         = specified title for superposed plot, if not superposed, plot title will be feature name
  # width        = width of png file if saving plot
  # height       = height of png file if saving plot
  
  suppressPackageStartupMessages(library(dplyr))
  
  # If x axis is setted to be cell ordered by rank and not by pseudotime value
  if(by.order){
    x$Pseudotime <- 1:nrow(x)
  }
  
  # Subset table for more efficiency during following steps
  x <- x[colnames(x) %in% c('Pseudotime', feature.list),]
  
  
  # Attributes Bin to each cell
  x <- PseudotimeRepartition(x, bin.number)
  
  # Initiate loop
  x.dim         <- as.numeric(names(table(x$Bin)))
  expr_table    <- data.frame()
  stdup_table   <- data.frame()
  stddown_table <- data.frame()
  absent_genes  <- c()
  
  # Generate new_table containing bin mean expression for each provided gene
  for(feature in feature.list){
    # Reset vectors for current feature
    mean.expr       <- c()
    stdup           <- c()
    stddown         <- c()
    if(feature %in% colnames(x)){
      for(b in x.dim){
        # Subset Bin
        values      <- x[,feature][as.character(x$Bin) == as.character(b)]
        # Scaling data if method 1 is precised
        if(scale & scale.method %in% 1){
          values <- (values-min(values))/(max(values)-min(values))
          values <- values %>% replace(is.na(.), 0)
          ylim      <- c(0,1)
        }
        
        # Calculate vectors
        mean.expr   <- c(mean.expr, mean(values))
        stdup       <- c(stdup, mean(values)+sd(values))
        stddown     <- c(stddown, mean(values)-sd(values))
      }
      
      # Scaling resulting means if scale method 0 is precised (default)
      if(scale & scale.method %in% 0){
        mean.expr <- (mean.expr-min(mean.expr))/(max(mean.expr)-min(mean.expr))
        mean.expr <- mean.expr %>% replace(is.na(.), 0)
        ylim      <- c(0,1)
      }else if(scale & scale.method %in% 2){
        mean.expr <- mean.expr/(abs(max(mean.expr))+abs(min(mean.expr)))
        mean.expr <- mean.expr %>% replace(is.na(.), 0)
        ylim      <- c(-1,1)
      }
      
      # Complete tables with current vectors
      expr_table    <- rbind(expr_table, c(feature, mean.expr))
      stdup_table   <- rbind(stdup_table, c(feature, stdup))
      stddown_table <- rbind(stddown_table, c(feature, stddown))
    }else{
      absent_genes <- c(absent_genes, feature)
    }
  }
  
  # Formation final new_table
  rownames(expr_table)    <- expr_table[,1]
  rownames(stdup_table)   <- stdup_table[,1]
  rownames(stddown_table) <- stddown_table[,1]
  expr_table              <- expr_table[,2:ncol(expr_table)]
  stdup_table             <- stdup_table[,2:ncol(stdup_table)]
  stddown_table           <- stddown_table[,2:ncol(stddown_table)]
  colnames(expr_table)    <- x.dim
  colnames(stdup_table)   <- x.dim
  colnames(stddown_table) <- x.dim
  
  # Define function to select colors if multiple curves are drawn
  choose_color <- function(defcol, n){
    if(length(defcol) == 1){
      return(defcol)
    }else{
      if(n <= length(defcol)){
        return(defcol[n])
      }else{
        return("lightgrey")
      }
    }
  }
  
  # Plotting depending on selected options
  if(write == FALSE){
    if(superposed == FALSE && compare.with == FALSE){
      for(feature in rownames(expr_table)){
        # Draw expression plot anyway
        plot(x.dim, expr_table[feature,], type = "l", lwd = lwd, ylim = ylim, col = col[1],  
             xlab = xlab, ylab = ylab, main = feature)
        if(std == TRUE){
          # Calculate and draw STD surfaces
          par(new=TRUE)
          x.surface <- c(x.dim, rev(x.dim))
          y.surface <- c(stdup_table[feature,], rev(stddown_table[feature,]))
          polygon(x.surface, y.surface, col = "#e9e9e9", border = NA)
          par(new=TRUE)
          # Draw expression curve again (upper than std surfaces)
          plot(x.dim, expr_table[feature,], type = "l", lwd = lwd, ylim = ylim, col = col[1],  
               xlab = xlab, ylab = ylab, main = feature)
          par(new=FALSE)
        }
      }
    }else if(superposed == TRUE && compare.with == FALSE){
      n <- 1
      for(feature in rownames(expr_table)){
        # Draw all expression plot
        plot(x.dim, expr_table[feature,], type = "l", lwd = lwd, ylim = ylim, col = choose_color(col, n),  
             xlab = xlab, ylab = ylab, main = main)
        par(new=TRUE)
        n <- n+1
      }
      par(new=FALSE)
    }else if(compare.with != FALSE){
      for(feature in rownames(expr_table)){
        n <- 1
        plot(x.dim, expr_table[compare.with,], type = "l", lwd = lwd, ylim = ylim, col = choose_color(col, n),  
             xlab = xlab, ylab = ylab, main = "")
        par(new=TRUE)
        n <- n+1
        plot(x.dim, expr_table[feature,], type = "l", lwd = lwd, ylim = ylim, col = choose_color(col, n),  
             xlab = xlab, ylab = ylab, main = feature)
        par(new=FALSE)
      }
    }
  }else if(write == TRUE){
    if(superposed == FALSE && compare.with == FALSE){
      for(feature in rownames(expr_table)){
        png(paste(dir, "/", feature, ".png", sep = ""), width = width, height = height)
        # Draw expression plot anyway
        plot(x.dim, expr_table[feature,], type = "l", lwd = lwd, ylim = ylim, col = col[1],  
             xlab = xlab, ylab = ylab, main = feature)
        if(std == TRUE){
          # Calculate and draw STD surfaces
          par(new=TRUE)
          x.surface <- c(x.dim, rev(x.dim))
          y.surface <- c(stdup_table[feature,], rev(stddown_table[feature,]))
          polygon(x.surface, y.surface, col = "#e9e9e9", border = NA)
          par(new=TRUE)
          # Draw expression curve again (upper than std surfaces)
          plot(x.dim, expr_table[feature,], type = "l", lwd = lwd, ylim = ylim, col = col[1],  
               xlab = xlab, ylab = ylab, main = feature)
          par(new=FALSE)
        }
        dev.off()
      }
    }else if(superposed == TRUE && compare.with == FALSE){
      n <- 1
      png(paste(dir, "/", main, ".png", sep = ""),width = 900, height = 800)
      for(feature in rownames(expr_table)){
        # Draw all expression plot
        plot(x.dim, expr_table[feature,], type = "l", lwd = lwd, ylim = ylim, col = choose_color(col, n),  
             xlab = xlab, ylab = ylab, main = main)
        par(new=TRUE)
        n <- n+1
      }
      par(new=FALSE)
      dev.off()
    }else if(compare.with != FALSE){
      for(feature in rownames(expr_table)){
        n <- 1
        png(paste(dir, "/", feature, ".png", sep = ""), width = width, height = height)
        plot(x.dim, expr_table[compare.with,], type = "l", lwd = lwd, ylim = ylim, col = choose_color(col, n),  
             xlab = xlab, ylab = ylab, main = "")
        par(new=TRUE)
        n <- n+1
        plot(x.dim, expr_table[feature,], type = "l", lwd = lwd, ylim = ylim, col = choose_color(col, n),  
             xlab = xlab, ylab = ylab, main = feature)
        par(new=FALSE)
        dev.off()
      }
    }  
  }    
  return(expr_table)
}     

#------------------------------------------------------------------------------#

CompareExpr <- function(x, bin.number = 21, feature = 'Nfil3', by.order = F,
                        scale = T, col = "black", lwd = 2, ylim = c(-0.75, 2), 
                        xlab = "Pseudotime", ylab = "", main = ''){
  
  # Draw feature expression by separating provided celltypes along pseudotime
  #
  # x            = ordered matrix with celltypes described in "Sample" row
  # bin.number   = number of bin for representation
  # feature      = feature to draw
  # by.order     = if set to TRUE, consider cell order rather than pseudotimeto make bins
  # scale        = whether scale shwon expression levels from 0 to 1
  # col          = colors for plots (several can be specified for superposed plot)
  # ylim         = ylim for plots (can be manually adjusted)
  # xlab         = label on x axis for plots (default = "Pseudoime")
  # ylab         = label on y axis for plots (default = "")
  # main         = specified title for superposed plot, if not superposed, plot title will be feature name
  
  
  # Get sample names by default
  sample.names <- names(table(as.character(as.vector(x$Sample))))
  
  # If x axis is setted to be cell ordered by rank and not by pseudotime value
  if(by.order){
    x <- x[,colnames(x) %!in% 'Pseudotime']
    x <- cbind(c(1:nrow(x)), x)
    colnames(x) <- c('Pseudotime', colnames(x)[2:ncol(x)])
  }
  
  # Attribute Bin to each cell
  x <- PseudotimeRepartition(x, bin.number)
  
  # Initiate loop
  x.dim         <- as.numeric(names(table(x$Bin)))
  expr_table    <- data.frame()
  
  
  if(feature %in% colnames(x)){
    # For each celltype it creates a new mean binned vector
    for(s in sample.names){
      subx          <- x[x$Sample %in% s, c('Bin', feature)]
      mean.expr     <- c()
      NA.list       <- c()
      for(b in x.dim){
        # Remove current point if number of values < 5
        values      <- as.numeric(subx[,feature][as.character(subx$Bin) %in% as.character(b)])
        if(length(values) > 5){
          mean.expr <- c(mean.expr, mean(values))
        }else{
          mean.expr <- c(mean.expr, NA)
        }
      }
      # Complete tables with current vectors
      expr_table    <- rbind(expr_table, c(s, mean.expr))
    }
    
    # Formating final table
    res_table <- as.data.frame(sapply(expr_table[,2:ncol(expr_table)],as.numeric))
    
    # If scale is precised
    if(scale){
      res_table <- (res_table-min(res_table, na.rm = T))/(max(res_table, na.rm = T)-min(res_table, na.rm = T))
      ylim      <- c(0,1)
    }
    
    rownames(res_table)    <- expr_table[,1]
    colnames(res_table)    <- x.dim
    res_table              <- res_table[,order(x.dim)]
    
    ### DRAW
    
    n <- 1
    # Add extra space to right of plot area; change clipping to figure
    par(mar=c(5,4,2,7))
    for(s in rownames(res_table)){
      # Draw all expression plot
      plot(x.dim, as.numeric(res_table[s,]), type = "l", lwd = lwd, ylim = ylim, 
           xlim=c(min(x$Pseudotime), max(x$Pseudotime)), col = col[n],  
           xlab = xlab, ylab = ylab, main = main)
      par(new=TRUE)
      n <- n+1
    }
    legend(x="topleft", inset=c(1,0), xpd=T, bty="n", legend = sample.names, col = col, lwd = lwd, text.font = 0.5)
    par(new=F)
  }else{
    warnings(paste(feature, 'gene not found.' ))
  }
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
    y <- y[, colnames(y) %!in% c("Gene.stable.ID", "Gene.stable.ID.version", "Gene.name")]
  } else {
    y$Gene.name <- make.names(y$Gene.name, unique = TRUE)
    
    if(sum.duplicate){
      z   <- data.frame()
      for(g in non_unique_names){
        pattern          <- paste("^", g, "|", g, "/.[123456789]", sep = "")
        sub_table        <- y[grepl(pattern, y$Gene.name),]
        sub_table        <- sub_table[, colnames(sub_table) %!in% c("Gene.stable.ID", "Gene.stable.ID.version", "Gene.name")]
        new_line <- c(paste(g, ".ALL", sep = ""), apply(sub_table, 2, sum))
        z <- rbind(z, new_line)
      }
      rownames(z) <- z[,1]
      z <- z[,2:ncol(z)]
      colnames(z) <- colnames(y)[colnames(y) %!in% c("Gene.stable.ID", "Gene.stable.ID.version", "Gene.name")]
    }
    rownames(y) <- y$Gene.name
    y <- y[, colnames(y) %!in% c("Gene.stable.ID", "Gene.stable.ID.version", "Gene.name")]
    y <- rbind(y, z)
  }
  return(y)
}  

#------------------------------------------------------------------------------#
