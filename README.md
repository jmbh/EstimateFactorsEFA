## Reproducibility Archive

This repository allows you to reproduce the simulation study, and all results and figures in the preprint "Estimating the Number of Factors in Exploratory Factor Analysis via out-of-sample Prediction Errors" (https://psyarxiv.com/qktsd)

### Simulation study

- `Simulation.R` contains the simulation script for one iteration of the simulation
- `aux_functions.R` contains a function to generate data from different factor models
- `submit_jobs.sh` is a batch script that runs a single iteration of `Simulation.R` on the LISA custer of UvA
- `submit_all.sh` is a batch script that submits `submit_jobs.sh` with seeds `1:200` to the LISA cluster

The output of the simulation study is in the folder /output. In principle the simulation can also be run locally, bt running `Simulation.R` sequentially with seeds `1:200`. Each iteration took around 3h, when running parallel on 12 cores. Note that some runs failed because of a rare convergence issue of an eigen-decomposition. We worked around this issue by adding additional runs.

### Results
- `Evaluation.R` preprocesses the simulation output, computes the numerical results shown in the paper and the results Figures 1 and 2, and the figures in Appendices H and I
- `Additional_Analyses.R` contains additional analyses reported in the paper that are not computed from the simulation results

### Tutorial
- `Tutorial.R` contains the code of the tutorial in Appendix A and creates the figure in that appendix


### sessionInfo()

R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 10 (buster)

Matrix products: default
BLAS/LAPACK: /sara/eb/AVX2/Debian10/EB_production/2021/software/FlexiBLAS/3.0.4-GCC-10.3.0/lib/libflexiblas.so.3.0

locale:
 [1] LC_CTYPE=en_US       LC_NUMERIC=C         LC_TIME=en_US       
 [4] LC_COLLATE=en_US     LC_MONETARY=en_US    LC_MESSAGES=en_US   
 [7] LC_PAPER=en_US       LC_NAME=C            LC_ADDRESS=C        
[10] LC_TELEPHONE=C       LC_MEASUREMENT=en_US LC_IDENTIFICATION=C 

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] doParallel_1.0.16    iterators_1.0.13     foreach_1.5.1       
 [4] GAIPE_1.0            corpcor_1.6.10       mvtnorm_1.1-1       
 [7] GPArotation_2022.4-1 psych_2.2.5          fspe_0.1.0          
[10] EGAnet_1.0.0         lavaan_0.6-8        

loaded via a namespace (and not attached):
 [1] codetools_0.2-18 lattice_0.20-44  grid_4.1.0       nlme_3.1-152    
 [5] stats4_4.1.0     pbivnorm_0.6.0   tools_4.1.0      compiler_4.1.0  
 [9] mnormt_2.0.2     tmvnsim_1.0-2   


### Package Version of EGAnet

We followed the author's advice and used the latest version of their package from Github. We added the source of that version in this reproducibility archive (EGAnet.zip)

