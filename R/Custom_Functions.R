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
