1. system requirements

R Session information: 
*********************************************************************
R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.9 (Final)

Matrix products: default
BLAS: /usr/lib64/libblas.so.3.2.1
LAPACK: /net/shendure/vol1/home/cao1025/anaconda3/lib/libmkl_rt.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] splines   parallel  stats4    stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] bindrcpp_0.2.2        pheatmap_1.0.10       rgl_0.99.16          
 [4] dplyr_0.7.5           monocle_2.9.0         L1Graph_0.1.0        
 [7] proxy_0.4-22          lpSolveAPI_5.5.2.0-17 igraph_1.2.1         
[10] sse_0.1.0             DDRTree_0.1.5         VGAM_1.0-5           
[13] ggplot2_2.2.1         Biobase_2.40.0        DelayedArray_0.6.0   
[16] BiocParallel_1.14.1   IRanges_2.14.10       S4Vectors_0.18.3     
[19] BiocGenerics_0.26.0   matrixStats_0.53.1    irlba_2.3.2          
[22] Matrix_1.2-14        

loaded via a namespace (and not attached):
 [1] nlme_3.1-137             gmodels_2.16.2           RColorBrewer_1.1-2      
 [4] repr_0.15.0              docopt_0.4.5             tools_3.5.0             
 [7] R6_2.2.2                 spData_0.2.8.3           lazyeval_0.2.1          
[10] colorspace_1.3-2         manipulateWidget_0.9.0   sp_1.3-1                
[13] tidyselect_0.2.4         gridExtra_2.3            compiler_3.5.0          
[16] glmnet_2.0-16            expm_0.999-2             labeling_0.3            
[19] slam_0.1-43              scales_0.5.0             pbdZMQ_0.3-3            
[22] stringr_1.3.1            digest_0.6.15            sparsesvd_0.1-4         
[25] base64enc_0.1-3          pkgconfig_2.0.1          htmltools_0.3.6         
[28] limma_3.36.1             htmlwidgets_1.2          rlang_0.2.1             
[31] DelayedMatrixStats_1.2.0 FNN_1.1                  shiny_1.1.0             
[34] bindr_0.1.1              jsonlite_1.5             gtools_3.5.0            
[37] crosstalk_1.0.0          spdep_0.7-7              magrittr_1.5            
[40] Rcpp_0.12.17             IRkernel_0.8.12.9000     munsell_0.4.3           
[43] viridis_0.5.1            reticulate_1.7           stringi_1.2.2           
[46] MASS_7.3-50              Rtsne_0.13               plyr_1.8.4              
[49] grid_3.5.0               gdata_2.18.0             promises_1.0.1          
[52] ggrepel_0.8.0            crayon_1.3.4             miniUI_0.1.1.1          
[55] deldir_0.1-15            lattice_0.20-35          IRdisplay_0.5.0         
[58] knitr_1.20               pillar_1.2.3             uuid_0.1-2              
[61] boot_1.3-20              reshape2_1.4.3           codetools_0.2-15        
[64] LearnBayes_2.15.1        glue_1.2.0               evaluate_0.10.1         
[67] httpuv_1.4.3             foreach_1.4.4            gtable_0.2.0            
[70] RANN_2.5.1               purrr_0.2.5              assertthat_0.2.0        
[73] mime_0.5                 xtable_1.8-2             coda_0.19-1             
[76] later_0.7.2              viridisLite_0.3.0        HSMMSingleCell_0.114.0  
[79] qlcMatrix_0.9.7          tibble_1.4.2             iterators_1.0.9         
[82] cluster_2.0.7-1          fastICA_1.2-1            densityClust_0.3  
*********************************************************************

2. Installation guide (Install time: 4 hours)
(1) Install R.3.5.0 from online source (https://www.r-project.org/)
(2) Install Anaconda (https://conda.io/docs/user-guide/install/index.html)
(3) install monocle dependent R packages:
For sse, L1-graph, and monocle, source code are included in the folder (Package_source) and can be installed by "R CMD INSTALL (package name)" Other dependent packages can be installed from CRAN  (The Comprehensive R Archive Network) or Bioconductor.
(4) Install UMAP package
UMAP is installed from the github source (https://github.com/lmcinnes/umap)
conda install -c conda-forge umap-learn