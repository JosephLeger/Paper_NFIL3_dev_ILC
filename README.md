# Paper_NFIL3_dev_ILC
Source code used in the paper : [The transcription factor NFIL3 drives ILC specification from lymphoid progenitors](https://pubmed.ncbi.nlm.nih.gov/41172988/)   
Léger J, Artano E, Coulais D, Belletoise N, Fadhloun R, Kenney D, Bhandoola A, Harly C.  
  
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
[1] LC_COLLATE=French_France.1252  LC_CTYPE=French_France.1252    LC_MONETARY=French_France.1252
[4] LC_NUMERIC=C                   LC_TIME=French_France.1252    

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] VennDiagram_1.7.3           futile.logger_1.4.3         igraph_1.3.5               
 [4] slingshot_2.2.1             TrajectoryUtils_1.2.0       SingleCellExperiment_1.16.0
 [7] princurve_2.1.6             Nebulosa_1.4.0              patchwork_1.1.2            
[10] simspec_0.0.0.9000          SeuratObject_4.1.3          Seurat_4.3.0               
[13] EnhancedVolcano_1.12.0      ggrepel_0.9.2               org.Mm.eg.db_3.14.0        
[16] AnnotationDbi_1.56.2        edgeR_3.36.0                limma_3.50.3               
[19] sva_3.35.2                  BiocParallel_1.28.3         genefilter_1.76.0          
[22] mgcv_1.8-39                 nlme_3.1-155                pheatmap_1.0.12            
[25] ggfortify_0.4.17            tximport_1.22.0             Hmisc_4.7-2                
[28] Formula_1.2-4               survival_3.2-13             lattice_0.20-45            
[31] forcats_0.5.2               stringr_1.5.0               dplyr_1.0.10               
[34] purrr_0.3.5                 readr_2.1.3                 tidyr_1.2.1                
[37] tibble_3.3.0                ggplot2_3.4.0               tidyverse_1.3.2            
[40] DiffBind_3.4.11             SummarizedExperiment_1.24.0 Biobase_2.54.0             
[43] MatrixGenerics_1.6.0        matrixStats_0.63.0          GenomicRanges_1.46.1       
[46] GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4           
[49] BiocGenerics_0.40.0        

loaded via a namespace (and not attached):
  [1] rtracklayer_1.54.0       scattermore_0.8          coda_0.19-4              bit64_4.0.5             
  [5] knitr_1.41               irlba_2.3.5.1            DelayedArray_0.20.0      data.table_1.14.6       
  [9] rpart_4.1.16             hwriter_1.3.2.1          KEGGREST_1.34.0          RCurl_1.98-1.9          
 [13] generics_0.1.3           lambda.r_1.2.4           cowplot_1.1.1            RSQLite_2.2.18          
 [17] RANN_2.6.1               future_1.29.0            bit_4.0.5                tzdb_0.3.0              
 [21] spatstat.data_3.0-0      xml2_1.3.3               lubridate_1.9.0          httpuv_1.6.6            
 [25] assertthat_0.2.1         gargle_1.2.1             amap_0.8-19              apeglm_1.16.0           
 [29] xfun_0.35                hms_1.1.2                promises_1.2.0.1         restfulr_0.0.15         
 [33] caTools_1.18.2           dbplyr_2.2.1             readxl_1.4.1             DBI_1.1.3               
 [37] htmlwidgets_1.5.4        spatstat.geom_3.0-3      googledrive_2.0.0        ellipsis_0.3.2          
 [41] ks_1.14.0                backports_1.4.1          annotate_1.72.0          deldir_1.0-6            
 [45] vctrs_0.5.1              ROCR_1.0-11              abind_1.4-5              cachem_1.0.6            
 [49] withr_2.5.0              BSgenome_1.62.0          progressr_0.11.0         bdsmatrix_1.3-6         
 [53] checkmate_2.1.0          sctransform_0.3.5        GenomicAlignments_1.30.0 mclust_6.0.0            
 [57] goftest_1.2-3            cluster_2.1.2            lazyeval_0.2.2           crayon_1.5.2            
 [61] spatstat.explore_3.0-5   pkgconfig_2.0.3          vipor_0.4.5              nnet_7.3-17             
 [65] rlang_1.1.6              globals_0.16.2           lifecycle_1.0.3          miniUI_0.1.1.1          
 [69] extrafontdb_1.0          modelr_0.1.10            invgamma_1.1             polyclip_1.10-4         
 [73] ggrastr_1.0.1            cellranger_1.1.0         lmtest_0.9-40            Matrix_1.5-3            
 [77] ashr_2.2-54              zoo_1.8-11               reprex_2.0.2             base64enc_0.1-3         
 [81] beeswarm_0.4.0           ggridges_0.5.4           googlesheets4_1.0.1      png_0.1-8               
 [85] viridisLite_0.4.1        rjson_0.2.21             bitops_1.0-7             KernSmooth_2.23-20      
 [89] Biostrings_2.62.0        blob_1.2.3               mixsqp_0.3-48            SQUAREM_2021.1          
 [93] spatstat.random_3.0-1    ShortRead_1.52.0         parallelly_1.32.1        jpeg_0.1-10             
 [97] scales_1.2.1             memoise_2.0.1            magrittr_2.0.3           plyr_1.8.8              
[101] ica_1.0-3                gplots_3.1.3             zlibbioc_1.40.0          compiler_4.1.3          
[105] BiocIO_1.4.0             bbmle_1.0.25             RColorBrewer_1.1-3       ash_1.0-15              
[109] fitdistrplus_1.1-8       Rsamtools_2.10.0         cli_3.4.1                systemPipeR_2.0.8       
[113] XVector_0.34.0           listenv_0.8.0            pbapply_1.6-0            formatR_1.12            
[117] htmlTable_2.4.1          MASS_7.3-55              tidyselect_1.2.0         stringi_1.7.8           
[121] proj4_1.0-12             emdbook_1.3.13           yaml_2.3.6               locfit_1.5-9.6          
[125] latticeExtra_0.6-30      tools_4.1.3              timechange_0.1.1         future.apply_1.10.0     
[129] parallel_4.1.3           rstudioapi_0.14          foreign_0.8-82           gridExtra_2.3           
[133] Rtsne_0.16               digest_0.6.30            pracma_2.4.2             shiny_1.7.3             
[137] Rcpp_1.0.9               broom_1.0.1              ggalt_0.4.0              later_1.3.0             
[141] RcppAnnoy_0.0.20         httr_1.4.4               colorspace_2.0-3         tensor_1.5              
[145] rvest_1.0.3              XML_3.99-0.13            fs_1.5.2                 reticulate_1.26         
[149] truncnorm_1.0-9          splines_4.1.3            uwot_0.1.14              spatstat.utils_3.0-1    
[153] sp_1.5-1                 plotly_4.10.1            xtable_1.8-4             futile.options_1.0.1    
[157] jsonlite_1.8.3           R6_2.5.1                 pillar_1.10.2            htmltools_0.5.3         
[161] mime_0.12                glue_1.8.0               fastmap_1.1.0            codetools_0.2-18        
[165] maps_3.4.1               GreyListChIP_1.26.0      mvtnorm_1.1-3            spatstat.sparse_3.0-0   
[169] numDeriv_2016.8-1.1      ggbeeswarm_0.6.0         leiden_0.4.3             gtools_3.9.4            
[173] Rttf2pt1_1.3.11          interp_1.1-3             munsell_0.5.0            GenomeInfoDbData_1.2.7  
[177] haven_2.5.1              reshape2_1.4.4           gtable_0.3.1             extrafont_0.18
```
