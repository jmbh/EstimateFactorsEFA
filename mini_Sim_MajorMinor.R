# jonashaslbeck@gmail.com; Dec 17, 2021

# --------------------------------------------------------------
# ---------- What are we doing here? ---------------------------
# --------------------------------------------------------------

# Investigate drop in performance in methods when drastically
# decreasing the size of the factor loadings of one factor

# --------------------------------------------------------------
# ---------- Load Packages -------------------------------------
# --------------------------------------------------------------

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



# --------------------------------------------------------------
# ---------- Aux Functions -------------------------------------
# --------------------------------------------------------------

datagenMM_sim <- function(exp=0, n) {

  # Draw factor loadings
  nfactor <- 6
  items_pf <- 6 # items / factor
  psi <- 0.40
  l_loadings <- list()

  for(i in 1:nfactor) l_loadings[[i]] <- rep(0.65, items_pf)
  if(exp==1) l_loadings[[nfactor]] <- rep(0.65/2, items_pf)

  # Each item loads on one factor
  spec_items <- rep(NA, nfactor)
  for(i in 1:nfactor) {
    base <- (0:(nfactor-1))[i] * items_pf
    spec_items[i] <- paste0("F", i, " =~ ", paste0(l_loadings[[i]], "*x", (base+1):(base + items_pf), collapse = " + "))
  }
  spec_items <- paste(spec_items, collapse = " \n ")
  full_model <- spec_items

  # Create factor correlations
  if(nfactor>1) {
    nfac_combn <- combn(1:nfactor, 2)
    n_combn <- ncol(nfac_combn)
    spec_psi <- rep(NA, n_combn)
    for(i in 1:n_combn) spec_psi[i] <- paste0("F", nfac_combn[1, i], " ~~", psi," * F", nfac_combn[2, i], " \n")
    spec_psi <- paste(spec_psi, collapse = "")

    full_model <- paste(spec_items, "\n", spec_psi)
  }

  # Sample data
  data <- simulateData(full_model,
                       sample.nobs = n)

  return(data)

} # eoF

# --------------------------------------------------------------
# ---------- Simulation ----------------------------------------
# --------------------------------------------------------------

# ----- Simulation Setup ------
v_n <- round(exp(seq(4.61, 8.5171, length=12)))
n_v_n <- length(v_n)
nIter <- 100
a_res <- array(NA, dim = c(nIter, n_v_n, 2, 6)) # iterations, n-variations, exp cases, methods
maxK <- 10

set.seed(1)

for(i in 1:nIter) {
  for(n in 1:n_v_n) {
    for(e in 1:2) {

      ## Generate Data
      data <- datagenMM_sim(exp = c(0,1)[e], n = v_n[n])

      ## Estimate using top 5 methods

      # PE
      k_PE_10 <- fspe(data = data,
                      maxK = maxK,
                      nfold = 10,
                      rep = 1,
                      method = "PE",
                      pbar = FALSE)
      a_res[i, n, e, 1] <- k_PE_10$nfactor


      # Parallel
      fa_fits_out2 <- fa.parallel(x = data,
                                  n.iter = 20,
                                  fa = "fa",
                                  fm = "minres",
                                  plot = FALSE,
                                  quant = 0.95)
      a_res[i, n, e, 2] <- fa_fits_out2$nfact


      # EGA
      out_EGA <- EGA(data = data, plot.EGA = FALSE, verbose = FALSE)
      if(is.na(out_EGA$n.dim)) out_EGA$n.dim <- 0
      a_res[i, n, e, 3] <- out_EGA$n.dim

      # CovE
      k_COV_10 <- fspe(data = data,
                       maxK = maxK,
                       nfold = 10,
                       rep = 1,
                       method = "Cov", pbar = FALSE)
      a_res[i, n, e, 4] <- k_COV_10$nfactor

      # BIC
      fa_fits_out <- vss(x = data,
                         rotate = "oblimin",
                         n = maxK,
                         plot = FALSE)
      a_res[i, n, e, 6] <- which.min(fa_fits_out$vss.stats$BIC)

      # AIC
      AIC_seq <- fa_fits_out$vss.stats$chisq - 2*fa_fits_out$vss.stats$dof
      k_AIC <- which.min(AIC_seq)
      a_res[i, n, e, 5] <- k_AIC

      print(paste0("i = ", i, " n = ", n, " exp = ", c(0,1)[e]))

    } # end for: 2

  } # end for: n

} # end for: i


# saveRDS(a_res, "files/SimRes_MM_withAIC.RDS")
a_res <- readRDS("files/SimRes_MM_withAIC.RDS")


# --------------------------------------------------------------
# ---------- Make nice Figure for paper ------------------------
# --------------------------------------------------------------


# Color scheme
library(RColorBrewer)
cols <- brewer.pal(5, "Set1")
v_methods <- c("PE", "Parallel", "EGA", "CovE", "AIC", "BIC")

# ----- Compute accuracy -----
a_res_ind <- a_res == 6
a_res_ind_agg <- apply(a_res_ind, 2:4, function(x) mean(x, na.rm=TRUE))

# ----- Plotting -----

pdf("figures/Fig_ExtraSim_MM.pdf", width=7, height=5)

par(mfrow=c(2,3))
for(i in 1:5) {

  # Canvas
  plot.new()
  plot.window(xlim=c(1,12), ylim=c(0,1))
  title(main = v_methods[i], font.main=1)
  axis(2, las=2)
  axis(1, labels = v_n, at=1:12, las=2, cex=.5)
  if(i %in% c(1,4)) title(ylab="Accuracy")
  if(i %in% c(3:5)) title(xlab="Sample size")

  # Data
  lines(a_res_ind_agg[, 1, i], col=cols[i], lwd=2) # Standard
  lines(a_res_ind_agg[, 2, i], col=cols[i], lwd=2, lty=2) # With minor factor
}

Legend
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))

legend("center",
       legend=c("6 Major", "5 Major, 1 Minor"),
       lty=1:2,
       lwd=c(2,2), bty="n", cex=1.2)

dev.off()





