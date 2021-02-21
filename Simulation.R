# jonashaslbeck@gmail.com; Feb 21, 2021

# --------------------------------------------------------------
# ---------- Get Iteration Number ------------------------------
# --------------------------------------------------------------

# !/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
print(iter)
iter <- as.numeric(iter)


# -----------------------------------------------------
# -------- Load Packages & Aux functions --------------
# -----------------------------------------------------

# -------- Data generation ----------------------------
library(lavaan)

# -------- Algorithms ----------------------------
library(EGAnet) # for EGA approach
library(fspe)
library(psych)
library(GPArotation)
library(mvtnorm)
library(corpcor)

# Parallel
library(foreach)
library(parallel)
library(doParallel)


# -------- Source aux functions -----------------------
source("aux_functions.R")


# -----------------------------------------------------
# -------- Simulation Specs ---------------------------
# -----------------------------------------------------

# Specifications: Data generation
v_nfactor <- c(1, 2, 4, 6)
v_items_pf <- c(3, 6, 12)
v_psi <- c(0, .4, .8)

# check whether I have enough n in lowest n-condition
k <- max(v_nfactor)
p <- max(v_items_pf)
k*p + k*(k-1)/2

v_n <- round(exp(seq(4.61, 8.5171, length=12)))
v_n

# -----------------------------------------------------
# -------- Algorithm Specs ----------------------------
# -----------------------------------------------------

# Specifications: Algorithm
maxK_global <- 10

# -----------------------------------------------------
# -------- Simulate -----------------------------------
# -----------------------------------------------------

# Setup parallelization
cluster <- 12
cl <- makeCluster(cluster, outfile="")
registerDoParallel(cl)


timer_total <- proc.time()[3]

out <- foreach(ni = 1:12,
               .packages = c("psych", "lavaan", "corpcor", "RGenData",
                             "mvtnorm", "GPArotation", "EGAnet", "fspe"),
               .export = c("f_datagen2", "f_impcor_fa", "FactorEst", "maxK_global"),
               .verbose = TRUE) %dopar% {

                 # Reproducibility
                 the_seed <- iter * 1000 + ni
                 set.seed(the_seed)

                 # Storage
                 a_out <- array(NA, dim=c(4,3,3,10))

                 for(nf in 1:4) {
                   for(it in 1:3) {
                     for(ps in 1:3) {

                       maxK <- maxK_global

                       # ----- get cell specs ------
                       nfactor <- v_nfactor[nf]
                       items_pf <- v_items_pf[it]
                       psi <- v_psi[ps]
                       n <- v_n[ni]

                       # number of variables
                       p <- items_pf * nfactor

                       # Calculate max number of factors in cell
                       maxK_cell <- NA
                       for(q in 1:maxK_global) {
                         k <- q
                         df <- p*(p-1)/2 - (k*p - k*(k-1)/2)
                         if(df>0) maxK_cell <- q else break
                       }

                       maxK <- maxK_cell

                       # Do not run the cell with factor=1, item/factor=3, because no compariosn is permitted
                       if(!(nf==1 & it == 1)) {

                         # ----- simulate data ------
                         timer_spec <- proc.time()[3]
                         dg_out <- f_datagen2(nfactor = nfactor,
                                              items_pf = items_pf,
                                              psi = psi,
                                              n = n,
                                              ra = c(0.3, 1))

                         data <- dg_out$data

                         timer_spec2 <- round(proc.time()[3] - timer_spec)

                         print(paste("iter = ", iter, " ni = ", ni, " nf = ", nf, " it = ", it, " ps = ", ps,
                                     "datagen: time = ", timer_spec2))


                         # ------ 1+2+3) OoS PE & COV w CV  ------

                         timer_spec <- proc.time()[3]

                         # PE + 10fold
                         k_PE_10 <- fspe(data = data,
                                         maxK = maxK,
                                         nfold = 10,
                                         rep = 1,
                                         method = "PE")

                         # COV + 10fold
                         k_COV_10 <- fspe(data = data,
                                          maxK = maxK,
                                          nfold = 10,
                                          rep = 1,
                                          method = "Cov")

                         # COV + 2fold
                         # For 2 fold, for small n, one often gets a not-pos-def matrix
                         # in the training set, which creates an error;
                         # Since we only ran this for comparison of the nfold-
                         # hyperparameter, we limit these runs to ni>4
                         if(ni > 4) {


                           k_COV_2 <- fspe(data = data,
                                           maxK = maxK,
                                           nfold = 2,
                                           rep = 1,
                                           method = "Cov")

                         } else {

                           k_COV_2 <- list()
                           k_COV_2$nfactor <- NA

                         }


                         timer_spec2 <- round(proc.time()[3] - timer_spec)
                         print(paste("iter = ", iter, " ni = ", ni, " nf = ", nf, " it = ", it, " ps = ", ps,
                                     "OoS PE: time = ", timer_spec2))


                         # ------ 4) Parallel  ------

                         timer_spec <- proc.time()[3]

                         options(mc.cores = 1) # work around weird psych-package bug
                         fa_fits_out2 <- fa.parallel(x = data,
                                                     n.iter = 20,
                                                     fa = "fa",
                                                     fm = "minres",
                                                     plot = FALSE,
                                                     quant = 0.95)

                         k_PA <- fa_fits_out2$nfact

                         timer_spec2 <- round(proc.time()[3] - timer_spec)

                         print(paste("iter = ", iter, " ni = ", ni, " nf = ", nf, " it = ", it, " ps = ", ps,
                                     "parallel: time = ", timer_spec2))

                         # ------ 5) EGA  ------

                         timer_spec <- proc.time()[3]
                         out_EGA <- EGA(data = data, plot.EGA = FALSE, verbose = FALSE)
                         if(is.na(out_EGA$n.dim)) out_EGA$n.dim <- 0
                         k_EGA <- out_EGA$n.dim



                         # ------ 6) Kaiser-Guttman ------

                         k_KG <- sum(fa_fits_out2$fa.values > 1)

                         # ------ 7) MAP ------

                         timer_spec <- proc.time()[3]
                         options(mc.cores = 1) # work around weird psych-package bug
                         fa_fits_out <- vss(x = data,
                                            rotate = "oblimin",
                                            n = maxK,
                                            plot = FALSE)
                         timer_spec2 <- round(proc.time()[3] - timer_spec)

                         print(paste("iter = ", iter, " ni = ", ni, " nf = ", nf, " it = ", it, " ps = ", ps,
                                     "psych VVS: time = ", timer_spec2))


                         k_MAP <- which.min(fa_fits_out$map)

                         # ------ 8) BIC  ------

                         k_BIC <- which.min(fa_fits_out$vss.stats$BIC)


                         # ------ 9) VSS ------

                         k_VSS <- which.max(fa_fits_out$vss.stats$cfit.1)


                         # ------ 10) RMSEA ------

                         k_RMSEA <- which.min(fa_fits_out$vss.stats$RMSEA)


                         # ------ Save  ------

                         a_out[nf, it, ps, 1]  <- k_PE_10$nfactor
                         a_out[nf, it, ps, 2]  <- k_COV_10$nfactor
                         a_out[nf, it, ps, 3]  <- k_COV_2$nfactor
                         a_out[nf, it, ps, 4]  <- k_PA
                         a_out[nf, it, ps, 5]  <- k_EGA
                         a_out[nf, it, ps, 6]  <- k_KG
                         a_out[nf, it, ps, 7]  <- k_MAP
                         a_out[nf, it, ps, 8]  <- k_BIC
                         a_out[nf, it, ps, 9]  <- k_VSS
                         a_out[nf, it, ps, 10]  <- k_RMSEA

                       } # end if: permitted cells

                     } # end for: nf
                   } # end for: it
                 }# end for: psu

                 return(a_out)

               } # end:foreach

# print total time of nodes
print(paste0("Full Timing Iteration ", iter, ":"))
proc.time()[3] - timer_total

stopCluster(cl)

# -----------------------------------------------------
# -------- Postprocess & Save -------------------------
# -----------------------------------------------------

# Combine n-variations
a_out_all <- array(NA, dim=c(12, 4, 3, 3, 10)) # Storage
for(i in 1:12) a_out_all[i, , , , ] <- out[[i]]

# Output file
saveRDS(a_out_all, file = paste0("Simres_Iter_lowFL3_", iter, ".RDS"))



