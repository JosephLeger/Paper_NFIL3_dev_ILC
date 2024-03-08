#!/usr/bin/env Rscript

#===============================================================================
## DESCRIPTION -----------------------------------------------------------------
#===============================================================================
# Perform DiffBind using data pre-filtered with Clumpify
# Peak Calling was performed using MACS2
# DiffBind FRD cutoff is setted to 0.1 (default = 0.05)



#===============================================================================
## SETUP -----------------------------------------------------------------------
#===============================================================================

rm(list=ls(all.names=TRUE))

PATH         <- 'C:/Users/E15639P/Doctorat/NFIL3_dev_ILC'
BAM_folder   <- 'C:/Users/E15639P/Data/Epigenetic/DNAse-seq/mm39/BAM'
BED_folder   <- 'C:/Users/E15639P/Data/Epigenetic/DNAse-seq/mm39/Peaks/MACS2'
SAMPLE_SHEET <- 'C:/Users/E15639P/Data/Epigenetic/DNAse-seq/mm39/SampleSheet_DNAse.csv'
setwd(PATH)

# Load required packages and functions
source('C:/Users/E15639P/Doctorat/NFIL3_dev_ILC/Custom_functions.R')
library(DiffBind)
library(tidyverse)
library(stats)

# Create required directories
dir.create(file.path(PATH, 'Figures'))
dir.create(file.path(paste0(PATH, '/Figures'), 'DNAse-seq'))
dir.create(file.path(PATH, 'Saves'))
dir.create(file.path(paste0(PATH, '/Saves'), 'DNAse-seq'))
dir.create(file.path(paste0(PATH, '/Saves/DNAse-seq'), 'Peaks'))

PATH_FIG  <- paste0(PATH, '/Figures/DNAse-seq')
PATH_SAVE <- paste0(PATH, '/Saves/DNAse-seq')



#===============================================================================
## READ FILES ------------------------------------------------------------------
#===============================================================================

# Read files
SHEET          <- read.csv(SAMPLE_SHEET)
SHEET$bamReads <- paste(BAM_folder, SHEET$bamReads, sep = '/')
SHEET$Peaks    <- paste(BED_folder, SHEET$Peaks, sep = '/')
SHEET$ScoreCol <- rep(5, nrow(SHEET))



#===============================================================================
## DATA PROCESSING -------------------------------------------------------------
#===============================================================================

## Object generation -----------------------------------------------------------
# Generate DBA object
data <- dba(sampleSheet=SHEET)
plot(data)
dev.print(png, file = paste0(PATH_FIG, "/data_0.png"), 
          width = 9, height = 9, units = 'in', res = 100)

# Calculate Binding Matrix and explore it
summary(data$binding[,3]-data$binding[,2]) 
data <- dba.count(data,
                  bUseSummarizeOverlaps = TRUE, 
                  summits = 100, 
                  score = DBA_SCORE_RPKM)

# Saving dba.count object
saveRDS(data, paste0(PATH_SAVE, '/data_DBAcount.RDS'))


plot(data)
dev.print(png, file=paste0(PATH_FIG, '/data_1.png'), 
          width=9, height=9, units='in', res=100)

info               <- dba.show(data)
libsizes           <- cbind(LibReads = info$Reads, 
                            FRiP = info$FRiP, 
                            PeakReads = round(info$Reads * info$FRiP))
rownames(libsizes) <- paste0(info$Condition, info$Replicate)


## Normalization and contrast analysis -----------------------------------------
COMP <- data.frame(Ctrl = c('ALP', 'EILP'), 
                   Tst = c('EILP', 'ILCP'))

# Normalization
data <- dba.normalize(data)

# Set auto threshold to 0.25 
FDR  <- 0.25                  
data$config$th <- FDR        

# Model design and contrast
for(i in 1:nrow(COMP)){
  data <- dba.contrast(data, 
                       categories = DBA_CONDITION,
                       minMembers = 3,
                       contrast = c('Condition', COMP$Ctrl[i], COMP$Tst[i]))
}

# Analysis
data <- dba.analyze(data)
dba.show(data, bContrasts = TRUE)

# PCA plot using all regions
dba.plotPCA(data,
            DBA_FACTOR,
            label = DBA_CONDITION)
dev.print(png, paste0(PATH_FIG, '/PCA_plot.png'),
          width = 9, height = 9, units = 'in', res = 100)

# Saving dba.analyze object
saveRDS(data, paste0(PATH_SAVE, '/data_DBAanalyze.RDS'))


## Full Table Preparation ------------------------------------------------------
# Generates big table with all merged peaks
TABLE <- data.frame(chr = data$peaks[[1]]$Chr, 
                    start = data$peaks[[1]]$Start,
                    end = data$peaks[[1]]$End,
                    ID = paste0('Peak_',1:nrow(data$peaks[[1]])))

# Add scores of each sample
for(i in 1:length(data$peaks)){
  TABLE <- cbind(TABLE, data$peaks[[i]]$Score)
  colnames(TABLE)[ncol(TABLE)] <- paste0('score_', i)
}

Merged.bed <- cbind(TABLE[,1:4], rep(1, nrow(TABLE)), rep('+', nrow(TABLE)))

# Saving native TABLE object and all peaks in a bedfile
write.table(TABLE, paste0(PATH_SAVE, '/Peaks_Table.txt'), 
            quote = F, row.names = F)
write.table(Merged.bed, paste0(PATH_SAVE, "/Peaks/All_Overlap_Peaks.bed"),
            sep = "\t", quote = F, row.names = F, col.names = F)



#===============================================================================
## SAVING ALL COMPARISONS RESULTS ----------------------------------------------
#===============================================================================

for(i in 1:nrow(COMP)){
  curr <- paste0(COMP$Ctrl[i], '_vs_', COMP$Tst[i])
  
  dir.create(file.path(paste0(PATH_SAVE, '/Peaks'), curr))
  dir.create(file.path(PATH_FIG, curr))
  
  ## SAVING PLOTS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Volcano plot
  dba.plotVolcano(data, 
                  contrast = i)
  
  dev.print(png, paste(PATH_FIG, curr, 'VolcanoPlot.png', sep = '/'),
            width = 9, height = 9, units = 'in', res = 100)
  
  
  # Venn Diagram
  dba.plotVenn(data, 
               contrast = i, 
               bDB = TRUE, 
               bGain = TRUE, 
               bLoss = TRUE, 
               bAll=FALSE)
  dev.print(png, paste(PATH_FIG, curr, 'VennDiagram.png', sep = '/'),
            width = 9, height = 9, units = 'in', res = 100)

  ## WRITE PEAKSETS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Table of results generation
  data.DB         <- dba.report(data, contrast = i)
  scores          <- -10*(log10(data.DB$FDR))
  chr             <- as.character(data.DB@seqnames)
  start           <- data.DB@ranges@start
  end             <- data.DB@ranges@start+(data.DB@ranges@width-1)
  # Add common Peak_ID to trace back all peaks
  peak_id         <- c()
  for(j in 1:length(data.DB)){
    current_chr   <- chr[j]
    current_start <- start[j]
    line          <- TABLE[TABLE$chr %in% current_chr,]
    line          <- line[line$start %in% current_start,]
    peak_id       <- c(peak_id, line$ID)
  }
  data.DB@elementMetadata@listData[['ID']] <- peak_id
  
  # BED files
  data.DB.BED  <- data.frame(chr = chr,
                             start = start,
                             end = end, 
                             ID = peak_id, 
                             score = scores, 
                             strand = rep('+', length(data.DB)), 
                             fold = data.DB$Fold) 
  data.DB.loss <- data.DB.BED[data.DB.BED$fold > 0,]  
  data.DB.gain <- data.DB.BED[data.DB.BED$fold < 0,]
  
  # DESeq2 report
  data.DB.DE  <- data[['contrasts']][[i]][['DESeq2']][['de']]
  
  # Update full TABLE
  new_cols <- data.frame(chr = as.character(data.DB@seqnames),
                         start = data.DB@ranges@start, 
                         end = data.DB@ranges@start+(data.DB@ranges@width-1),
                         fold = data.DB@elementMetadata@listData[["Fold"]], 
                         score = scores)
  colnames(new_cols) <- c('chr', 'start', 'end', 
                          paste0(curr, '_fold'), 
                          paste0(curr, '_score'))
  
  TABLE <- merge(TABLE, new_cols, by = c('chr', 'start', 'end'), all = T)
  
  # Write files
  write.table(data.DB, 
              paste(PATH_SAVE, 'Peaks', curr, 'Report.txt', sep = '/'),
              quote = F)
  write.table(data.DB.BED, 
              paste0(PATH_SAVE, '/Peaks/', curr, '/', curr, '.bed'), 
              quote = F, col.names = F, row.names = F, sep = "\t")
  write.table(data.DB.DE, 
              paste(PATH_SAVE, 'Peaks', curr, 'Contrast.txt', sep = '/'), 
              quote = F)
  
  write.table(data.DB.loss, 
              paste0(PATH_SAVE, '/Peaks/', curr, '/', curr, '_loss.bed'), 
              quote = F, col.names = F, row.names = F, sep = "\t")
  write.table(data.DB.gain,
                paste0(PATH_SAVE, '/Peaks/', curr, '/', curr, '_gain.bed'), 
              quote = F, col.names = F, row.names = F, sep = "\t")
  
}


# Saving final TABLE object
write.table(TABLE, paste0(PATH_SAVE, '/Peaks_Table_Contrasted.txt'), 
            quote = F, row.names = F)
 

## Saving files for all non DB sites -------------------------------------------
# Create new folder
dir.create(file.path(paste0(PATH_SAVE, '/Peaks'), 'non_DB'))

# Filtering to keep lines with full NA for contrasts
all_folds  <- grep('fold', colnames(TABLE)) 
non_DB     <- TABLE %>% 
  filter_at(vars(colnames(TABLE[all_folds])),all_vars(is.na(.))) 

non_DB.BED <- cbind(non_DB[,1:3], 
                    ID = non_DB$ID,  
                    score = rep(1, nrow(non_DB)),
                    strand = rep('+', nrow(non_DB)))

# Save non DB BED file
write.table(non_DB.BED, 
            paste0(PATH_SAVE, '/Peaks/non_DB/non_DB.bed'), 
            quote = F, row.names = F, col.names = F, sep = "\t")



#===============================================================================
## GENERATE BED FILES ----------------------------------------------------------
#===============================================================================

FDR   <- 0.25
data  <- readRDS(paste0(PATH_SAVE, '/data_DBAanalyze.RDS'))
TABLE <- read.table(paste0(PATH_SAVE, '/Peaks_Table_Contrasted.txt'), 
                    header = T)

## Precise peak profiles -------------------------------------------------------
# Vonlontary inversed ! ALP_vs_EILP < 0 means UP in EILP
patterns <- list('[0][0]' = c(0,0), 
                 '[+][0]' = c(-1,0),
                 '[-][+]' = c(1,-1),
                 '[0][-]' = c(0,1),
                 '[0][+]' = c(0,-1),
                 '[+][-]' = c(-1,1),
                 '[-][0]' = c(1,0),
                 '[+][+]' = c(-1,-1),
                 '[-][-]' = c(1,1))

for(n in names(patterns)){
  dir.create(file.path(paste0(PATH_SAVE, '//Peaks'), n))
}

# Generate bed file for each matching pattern
peak.groups <- list()

for(i in 1:length(patterns)){
  pat <- patterns[[i]]
  # first comp
  if(pat[1] == 0){
    subset <- TABLE[is.na(TABLE$ALP_vs_EILP_fold),]
  }else if(pat[1] == 1){
    subset <- TABLE[!is.na(TABLE$ALP_vs_EILP_fold),]
    subset <- subset[subset$ALP_vs_EILP_fold > 0,]
  }else if(pat[1] == -1){
    subset <- TABLE[!is.na(TABLE$ALP_vs_EILP_fold),]
    subset <- subset[subset$ALP_vs_EILP_fold < 0,]
  }
  # second comp
  if(pat[2] == 0){
    subset <- subset[is.na(subset$EILP_vs_ILCP_fold),]
  }else if(pat[2] == 1){
    subset <- subset[!is.na(subset$EILP_vs_ILCP_fold),]
    subset <- subset[subset$EILP_vs_ILCP_fold > 0,]
  }else if(pat[2] == -1){
    subset <- subset[!is.na(subset$EILP_vs_ILCP_fold),]
    subset <- subset[subset$EILP_vs_ILCP_fold < 0,]
  }
  
  peak.groups[[names(patterns)[i]]] <- subset
  BED <- subset[,c('chr', 'start', 'end', 'ID')]
  write.table(BED, 
              paste0(PATH_SAVE, '/Peaks/', names(patterns)[i], '/', 
                     names(patterns)[i], '.bed'), 
              quote = F, col.names = F, row.names = F, sep = "\t")
  print(paste(names(patterns)[i], ' : ', nrow(subset), ' peaks', sep = ''))
  
}


## All patterns of interest merged ---------------------------------------------

#Create new folder
dir.create(file.path(paste0(PATH_SAVE, '/Peaks'), 'POI_Merged'))

patterns_of_interest  <- c('[+][0]',
                           '[-][+]',
                           '[+][-]',
                           '[-][0]',
                           '[+][+]',
                           '[-][-]')

POI_peaks <- c()
for(p in patterns_of_interest){
  x <- read.table(paste0(PATH_SAVE, '/Peaks/', p, '/', p, '.bed'))
  POI_peaks <- c(POI_peaks, x[,4])
}

POI_merged     <- TABLE[TABLE$ID %in% POI_peaks,]

POI_merged.BED <- cbind(POI_merged[,1:3], 
                        ID = POI_merged$ID,  
                        score = rep(1, nrow(POI_merged)),
                        strand = rep('+', nrow(POI_merged)))

# Save ALP open BED file
write.table(POI_merged.BED, 
            paste0(PATH_SAVE, '/Peaks/POI_merged/POI_merged_sorted.bed'), 
            quote = F, row.names = F, col.names = F, sep = "\t")


## ALP open regions ------------------------------------------------------------
# Generates a file for all regions open in ALP as control for Motif enrichment

#Create new folder
dir.create(file.path(paste0(PATH_SAVE, '/Peaks'), 'ALP_open'))

# Read bed files containing ALP open regions
ALP_open_patterns <- c('[0][0]', 
                       '[0][-]', 
                       '[-][0]',
                       '[-][+]',
                       '[-][-]')

ALP_open_peaks    <- c()
for(p in ALP_open_patterns){
  x <- read.table(paste0(PATH_SAVE, '/Peaks/', p, '/', p, '.bed'))
  ALP_open_peaks <- c(ALP_open_peaks, x[,4])
}

ALP_open <- TABLE[TABLE$ID %in% ALP_open_peaks,]

ALP_open.BED <- cbind(ALP_open[,1:3], 
                      ID = ALP_open$ID,  
                      score = rep(1, nrow(ALP_open)),
                      strand = rep('+', nrow(ALP_open)))

# Save ALP open BED file
write.table(ALP_open.BED, 
            paste0(PATH_SAVE, '/Peaks/ALP_open/ALP_open_sorted.bed'), 
            quote = F, row.names = F, col.names = F, sep = "\t")



################################################################################
#   READY TO USE ALL GENERATED BED FILES FOR ENRICHMENT AND MOTIFS ANALYSIS    #
################################################################################



#===============================================================================
## TF MOTIFS IN OPEN CHROMATIN -------------------------------------------------
#===============================================================================
# Explore HOMER motif enrichment results
# Here are written only the first member of top 5 TF families identified, and
# then we add all members of its Subfamily

TF_fam <- list()

## MOTIFS OF [-][+] PATTERN ----------------------------------------------------

# 1) ETS1 (ETS)
# 2) RUNX (Runt)
# 3) AtGRF6 (GRF) /!\ vegetals
# 4) CTCF (Zf)
# 5) Ascl1 (bHLH)

TF_fam[["ETS"]]  <- c("Cets1","Cets2","Etv2","Gabpa","Fli1","Erg","Fev","Etv3",
                      "Erf","Ets1","Ets2")
TF_fam[["RUNT"]] <- c("Runx1","Runx2","Runx3")
TF_fam[["GRF"]]  <- c()
TF_fam[["ZF"]]   <- c("Ctcf","Ctcfl")
TF_fam[["BHLH"]] <- c("Ash1","Ash2","Ash3","Ash4","Ash5","Ash1l","Ash2l")


## MOTIFS OF [+][-] PATTERN ----------------------------------------------------

# 1) SpiB (ETS)
# 2) RUNX (Runt)
# 3) IRF8 (IRF)
# 4) CEBP (bZIP - CEBP related factors)
# 5) NFIL3 (PAR)

TF_fam[["IRF"]]  <- c("Spdef","Irf1","Irf2","Irf3","Irf4","Irf5","Irf6","Irf7",
                      "Irf8","Irf9","Tef1","Tef3","Tef4","Tef5")
TF_fam[["BZIP"]] <- c("Cebpa","Cebpb","Cebpg","Cebpd","Cebpe","Ddit3")
TF_fam[["PAR"]]  <- c("Dbp","Hlf","Tef","Tefa","Tefb","Nfil3")


## MOTIFS OF ALP OPEN PATTERN --------------------------------------------------

# 1) CTCF (Zf)
# 2) ETS1 (ETS)
# 3) RUNX (Runt)
# 4) IRF8 (IRF)
# 5) Ronin (THAP)

TF_fam[["THAP"]] <- c("Thap1","Thap2","Thap3","Thap4","Thap7","Thap11","Prkrir")


## BZIP LARGE FAMILY -----------------------------------------------------------
# Because BZIP factors group is a huge class, we add all BZIP-derived TF

TF_fam[["BZIP_large"]] <- c(
  # Jun factors
  "Cjun","Jun","Junb","Jund",
  # NF-E2-like factors
  "Nfe2","Nfe2l1","Nfe2l2","Nfe2l3","Bach1","Bach2", 
  # ATF-2-like factors
  "Atf2","Atf7","Creb5",
  # Fos factors
  "Cfos","Fos","Fosb","Fosl1", "Fosl2","Fra1","Fra2", 
  # ATF-3-like factors
  "Atf3","Jdp2", 
  # Large Maf factors
  "Cmaf", "Maf","Mafa","Mafb","Nrl",
  # Small Maf factors
  "Maff","Mafg","Mafk","Batf","Batf2","Batf3","Xbp1","Xbp1u","Xbp1s","Atf4",
  "Atf5",
  # CREB-like factors
  "Creb","Atf1","Crem","Cremt","Crema","Cremb","Cremg","Cremtag","Cremt2a",
  "Cremt1a","Cremt1ag","Cremt1g","Cremag",
  # CREB-3-like factors
  "Creb3","Lzip","Creb3l1","Creb3l2","Creb3l3","Creb3l4",
  # ATF-6
  "Atf6a","Atf6b",
  # CREBZF-like factors
  "Crebzf",
  # CREBL2-like factors
  "Crebl2",
  # C/EBP
  "Cebpa","Cebpb","Cebpg","Cebpd","Cebpe","Ddit3",
  # PAR factors
  "Dbp","Hlf","Tef","Tefa","Tefb","Nfil3",
  # TSC22
  "Tsc22d1","Tsc22d2","Tsc22d3","Tsc22d4",
  # UTF
  "Utf1","Spib"
)

# Saving R object containing family genes list
saveRDS(TF_fam, paste0(PATH_SAVE, '/TF_fam.RDS'))


