# NFIL3_dev_ILC
Source code used in the paper [The transcription factor NFIL3 drives ILC specification from lymphoid progenitors]()   
J. Léger, E. Artano, D. Coulais, N. Belletoise, R. Fadhloun, D. Kenney, A. Bhandoola, C. Harly  
  
## Description 
This repository contains scripts and resources used to perform DNase-seq, Bulk RNA-seq, scRNA-seq and CUT&RUN post-processing analyzes.  
Information about pre-processing steps are detailed in [PRE-PROCESSING.md](https://github.com/JosephLeger/NFIL3_dev_ILC/blob/main/pre-processing/README.md).  

Complete explanation is avaible ([https://github.com/JosephLeger/NFIL3_dev_ILC/Full_Tutorial/](https://www.youtube.com/watch?v=dQw4w9WgXcQ)).

## Session Info
```R
> sessionInfo()
R version 4.1.3 (2022-03-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=French_France.1252  LC_CTYPE=French_France.1252    LC_MONETARY=French_France.1252 LC_NUMERIC=C                  
[5] LC_TIME=French_France.1252    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] EnhancedVolcano_1.12.0      ggrepel_0.9.2               org.Mm.eg.db_3.14.0         AnnotationDbi_1.56.2       
 [5] edgeR_3.36.0                limma_3.50.3                sva_3.35.2                  BiocParallel_1.28.3        
 [9] genefilter_1.76.0           mgcv_1.8-39                 nlme_3.1-155                pheatmap_1.0.12            
[13] tximport_1.22.0             forcats_0.5.2               stringr_1.5.0               dplyr_1.0.10               
[17] purrr_0.3.5                 readr_2.1.3                 tidyr_1.2.1                 tibble_3.1.8               
[21] ggplot2_3.4.0               tidyverse_1.3.2             DiffBind_3.4.11             SummarizedExperiment_1.24.0
[25] Biobase_2.54.0              MatrixGenerics_1.6.0        matrixStats_0.63.0          GenomicRanges_1.46.1       
[29] GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4            BiocGenerics_0.40.0        

loaded via a namespace (and not attached):
  [1] readxl_1.4.1             backports_1.4.1          plyr_1.8.8               splines_4.1.3            amap_0.8-19             
  [6] digest_0.6.30            invgamma_1.1             htmltools_0.5.3          SQUAREM_2021.1           fansi_1.0.3             
 [11] magrittr_2.0.3           memoise_2.0.1            BSgenome_1.62.0          googlesheets4_1.0.1      tzdb_0.3.0              
 [16] Biostrings_2.62.0        annotate_1.72.0          extrafont_0.18           modelr_0.1.10            systemPipeR_2.0.8       
 [21] extrafontdb_1.0          bdsmatrix_1.3-6          timechange_0.1.1         jpeg_0.1-10              colorspace_2.0-3        
 [26] blob_1.2.3               rvest_1.0.3              apeglm_1.16.0            haven_2.5.1              crayon_1.5.2            
 [31] RCurl_1.98-1.9           jsonlite_1.8.3           survival_3.2-13          glue_1.6.2               gtable_0.3.1            
 [36] gargle_1.2.1             zlibbioc_1.40.0          XVector_0.34.0           proj4_1.0-12             DelayedArray_0.20.0     
 [41] Rttf2pt1_1.3.11          maps_3.4.1               scales_1.2.1             mvtnorm_1.1-3            DBI_1.1.3               
 [46] Rcpp_1.0.9               xtable_1.8-4             emdbook_1.3.13           bit_4.0.5                truncnorm_1.0-9         
 [51] htmlwidgets_1.5.4        httr_1.4.4               gplots_3.1.3             RColorBrewer_1.1-3       ellipsis_0.3.2          
 [56] pkgconfig_2.0.3          XML_3.99-0.13            dbplyr_2.2.1             deldir_1.0-6             locfit_1.5-9.6          
 [61] utf8_1.2.2               tidyselect_1.2.0         rlang_1.0.6              munsell_0.5.0            cellranger_1.1.0        
 [66] tools_4.1.3              cachem_1.0.6             cli_3.4.1                generics_0.1.3           RSQLite_2.2.18          
 [71] broom_1.0.1              fastmap_1.1.0            yaml_2.3.6               bit64_4.0.5              fs_1.5.2                
 [76] caTools_1.18.2           KEGGREST_1.34.0          ash_1.0-15               ggrastr_1.0.1            xml2_1.3.3              
 [81] compiler_4.1.3           rstudioapi_0.14          beeswarm_0.4.0           png_0.1-8                reprex_2.0.2            
 [86] stringi_1.7.8            ggalt_0.4.0              lattice_0.20-45          Matrix_1.5-3             vctrs_0.5.1             
 [91] pillar_1.8.1             lifecycle_1.0.3          bitops_1.0-7             irlba_2.3.5.1            rtracklayer_1.54.0      
 [96] R6_2.5.1                 BiocIO_1.4.0             latticeExtra_0.6-30      hwriter_1.3.2.1          ShortRead_1.52.0        
[101] KernSmooth_2.23-20       vipor_0.4.5              MASS_7.3-55              gtools_3.9.4             assertthat_0.2.1        
[106] rjson_0.2.21             withr_2.5.0              GenomicAlignments_1.30.0 Rsamtools_2.10.0         GenomeInfoDbData_1.2.7  
[111] parallel_4.1.3           hms_1.1.2                grid_4.1.3               coda_0.19-4              GreyListChIP_1.26.0     
[116] ashr_2.2-54              googledrive_2.0.0        mixsqp_0.3-48            bbmle_1.0.25             numDeriv_2016.8-1.1     
[121] lubridate_1.9.0          ggbeeswarm_0.6.0         interp_1.1-3             restfulr_0.0.15   
```
