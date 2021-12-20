## Reproducibility Archive for preprint https://psyarxiv.com/qktsd

This repository allows you to reproduce the simulation study, and all results and figures in the preprint "Estimating the Number of Factors in Exploratory Factor Analysis via out-of-sample Prediction Errors" (https://psyarxiv.com/qktsd)

### Simulation study

- `Simulation.R` contains the simulation script for one iteration of the simulation
- `aux_functions.R` contains a function to generate data from different factor models
- `submit_jobs.sh` is a batch script that runs a single iteration of `Simulation.R` on the LISA custer of UvA
- `submit_all.sh` is a batch script that submits `submit_jobs.sh` with seeds `1:200` to the LISA cluster

The output of the simulation study is in the folder /output. In principle the simulation can also be run locally, bt running `Simulation.R` sequentially with seeds `1:200`. Each iteration took around 30 minutes, when running parallel on 12 cores. Note that five runs failed because of a rare convergence issue of an eigen-decomposition. We worked around this issue by running five extra runs.

### Results
- `Evaluation.R` preprocesses the simulation output, computes the numerical results shown in the paper and th results figures 1 and 2
- `Additional_Analyses.R` contains additional analyses reported in the paper that are not computed from the simulation results
- `mini_Sim_MajorMinor.R` contains code to run the additional simulation study on minor/minor factors in the Appendix

### Tutorial
- `Tutorial.R` contains the code of the tutorial in Appendix A and creates the figure in that appendix


### sessionInfo()

R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 10 (buster)

Matrix products: default
BLAS:   /sara/eb/AVX2/Debian10/EB_production/2020/software/R/4.0.2-intel-2020a/lib/R/lib/libR.so
LAPACK: /sara/eb/AVX2/Debian10/EB_production/2020/software/R/4.0.2-intel-2020a/lib/R/modules/lapack.so

locale:
 [1] LC_CTYPE=en_US       LC_NUMERIC=C         LC_TIME=en_US       
 [4] LC_COLLATE=en_US     LC_MONETARY=en_US    LC_MESSAGES=en_US   
 [7] LC_PAPER=en_US       LC_NAME=C            LC_ADDRESS=C        
[10] LC_TELEPHONE=C       LC_MEASUREMENT=en_US LC_IDENTIFICATION=C 

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] doParallel_1.0.15     iterators_1.0.12      foreach_1.5.0        
 [4] corpcor_1.6.9         mvtnorm_1.1-0         GPArotation_2014.11-1
 [7] psych_2.0.9           fspe_0.1.0            EGAnet_0.9.8         
[10] lavaan_0.6-5         

loaded via a namespace (and not attached):
[1] codetools_0.2-16 lattice_0.20-41  grid_4.0.2       nlme_3.1-147    
[5] stats4_4.0.2     pbivnorm_0.6.0   tools_4.0.2      compiler_4.0.2  
[9] mnormt_1.5-6    


### Package Version of EGAnet

We followed the author's advice and used the latest version of their package from Github. We added the source of that version in this reproducibility archive (EGAnet.zip)

