# jonashaslbeck@gmail.com; May 30, 2022

# Data generation
library(lavaan)
source("aux_functions.R")
library(MASS)

# Estimation
library(fspe)
library(psych)
library(EGAnet)

# Visualization
library(RColorBrewer)


# --------------------------------------------------------------
# ---------- Get Average correlation in 1-Factor model ---------
# --------------------------------------------------------------

# Reported Mean/SD in Section 3.1

nIter <- 1000
m_save <- matrix(NA, nIter, 6*5/2)

set.seed(1)

for(i in 1:nIter) {

  dg_out <- f_datagen2(nfactor = 1,
                       items_pf = 6,
                       psi = 0,
                       n = 5000, ra = c(.3, 1))

  corm <- cor(dg_out$data)
  m_save[i, ] <- corm[upper.tri(corm)]

  print(i)

}

round(mean(corm), 2)
round(sd(corm), 2)


# --------------------------------------------------------------
# ---------- Appendix B: CV folds ------------------------------
# --------------------------------------------------------------

# Here we check whether the fold size leads to an increase
# We pick the cell with 4 factors, 3 items, and moderate correlation
# because we would like to pick a cell in which the performance is not
# too high already so the potential impact of fold size would be quite
# small

# We compare fold sizes of 5, 10, 20, 50 and run 100 repetitions of the design

# Settings
nIter <- 1000
v_folds <- c(2, 5, 10, 20)
n_v_folds <- length(v_folds)

# Storage
a_perf <- array(NA, dim=c(nIter, 4, 2))

# Reproducibility
set.seed(1)

# Loop through
for(i in 1:nIter) {

  dg_out <- f_datagen2(nfactor = 4,
                       items_pf = 3,
                       psi = 0.4,
                       n = 593)
  data <- dg_out$data

  for(fold in 1:n_v_folds) {

    # PE method
    out <- fspe(data = data,
                maxK = 7,
                nfold = v_folds[fold],
                method = "PE",
                pbar=FALSE)
    a_perf[i, fold, 1] <- out$nfactor

    # PE method
    out <- fspe(data = data,
                maxK = 7,
                nfold = v_folds[fold],
                method = "Cov",
                pbar=FALSE)
    a_perf[i, fold, 2] <- out$nfactor

  } # end: folds

  print(i)

} # end: iter


# Save
saveRDS(a_perf, file="files/FoldSim_1000iter.RDS")
# a_perf <- readRDS(file="files/FoldSim_1000iter.RDS")

# Accuracy
acc <- apply(a_perf==4, 2:3, mean)

# ----- Make Figure -----

# Same colors as Figure 2/3 in main text
cols <- brewer.pal(5, "Set1")[c(1,4)]

pdf("figures/Fig_FoldSim.pdf", width = 5, height = 4)
par(mar=c(4,4,3,1))
plot(acc[,1], ylim=c(.4,.75), type="l", axes=F, xlab="", ylab="", col=cols[1], lwd=2, lty=2)
points(acc[,1], col=cols[1], lwd=2, pch=20, cex=1.5)
points(acc[,2], col=cols[2], lwd=2, pch=20, cex=1.5)
lines(acc[,2], col=cols[2], lwd=2, lty=2)
axis(1, labels=v_folds, at=1:4)
axis(2, las=2)
title(ylab="Mean Accuracy")
title(xlab="Folds in Cross Validation")
legend("bottomright", legend=c("PE", "CovE"), lwd=c(2,2), col=cols, bty="n")
dev.off()


# --------------------------------------------------------------
# ---------- Appendix C: CV scheme repetitions -----------------
# --------------------------------------------------------------

# Question: can we bump up the performance by running CVs several
# times and choosing the most chosen factor number?

nIter <- 100
m_reps <- matrix(NA, nIter, 11)
set.seed(2)

# With this model: factor model with max 7 factors is identified

for(i in 1:nIter) {

  # Generate Data
  dg_out <- f_datagen2(nfactor = 4,
                       items_pf = 3,
                       psi = .4,
                       n = 593)

  # OoS PE scheme
  time <- proc.time()[3]
  out <- fspe(data = dg_out$data,
              maxK = 6,
              nfold = 10,
              rep = 10,
              method = "PE",
              pbar = FALSE)
  print(proc.time()[3] - time)

  # Save performances
  for(r in 1:10) {
    v_PEs <- apply(out$PE_array[, , , r], 1, mean)
    m_reps[i, r] <- which.min(v_PEs)
  }

  m_reps[i, 11] <- out$nfactor # based on majority rule

  print(i)
} # end iter

# saveRDS(m_reps, file="files/m_reps_seed2.RDS")
m_reps <- readRDS(file="files/m_reps_seed1.RDS")

# Accuracy improvement
m_acc <- m_reps == 4
v_acc <- colMeans(m_acc)
mean(v_acc[1:10]) # mean accuracy of individual runs
v_acc[11] # accuracy of majority rule

v_acc[11] / mean(v_acc[1:10])


# Checking free parameters

# ----- get cell specs ------
nfactor <- 4
items_pf <- 3
psi <- 0.4
n <- v_n[6]

# number of variables
p <- items_pf * nfactor

# Correlation between re-runs
corm <- cor(m_reps[, 1:10])
mean(corm[lower.tri(corm)])



# --------------------------------------------------------------
# ---------- Figures in Appendix D (Parallel Analysis) ---------
# --------------------------------------------------------------
# Code Adapted from rietvanbork@hotmail.com

set.seed(7)

# ---------- Figure 4 ---------------------------------------

dg_out <- f_datagen2(nfactor = 4,
                     items_pf = 3,
                     psi = 0,
                     n = 1000,
                     ra = c(0.3, 1))
p <- 4*3

fa_fits_out2 <- fa.parallel(x = dg_out$data,
                            n.iter = 20,
                            fa = "fa",
                            fm = "minres",
                            plot = T,
                            quant = 0.95)

sim_factor <- fa_fits_out2$values[ ,(p+1):(p+1+p-1)]
qntls <- apply(sim_factor, 2, function(x) quantile(x, probs = .95))
means <- apply(sim_factor, 2, function(x) mean(x))


pdf("figures/APP_paralleltest_4factors.pdf",width=6,height=5)
plot(qntls, xaxt="n",type="l",cex.axis = 1.2, col="red", cex.lab=1.2, lwd=2, ylab="eigen values",las=1,lty=2, ylim=c(0, 2),xlab="eigen value number ordered by absolute value", bty="n",main="Parallel test: 4 factors")
lines(fa_fits_out2$fa.values, type="b", lwd=2)
lines(means, col="red", lwd=2)
axis(side=1, at=seq(1,12),cex.axis = 1.2)
legend(8,2,
       legend = c("Actual data", "Mean", "95 quantile"),
       col = c("black", "red","red"),
       bty = "n",
       lty = c(1,1,2),
       pt.cex = 1,
       cex = 1.2,
       lwd = 2,
       text.col = "black")
dev.off()


# ---------- Figure 5 ---------------------------------------


v_n <- c(10000, 50000, 100000)

for(i in 1:3) {

  corr <- 0.3 #with .3 it always selects 0 factors and with .4 it selects 2 factors
  sig<-matrix(0,6,6)
  diag(sig)<-1
  sig[1,2]<-corr
  sig[2,1]<-corr
  sig[1,3]<-corr
  sig[3,1]<-corr
  sig[3,2]<-corr
  sig[2,3]<-corr
  sig[5,4]<-corr
  sig[4,5]<-corr
  sig[4,6]<-corr
  sig[6,4]<-corr
  sig[6,5]<-corr
  sig[5,6]<-corr
  sig
  p<-6
  data <- mvrnorm(n = v_n[i],
                  mu=rep(0,6),
                  Sigma=sig,
                  empirical=T)

  fa_fits_out2 <- fa.parallel(x = data,
                              n.iter = 20,
                              fa = "fa",
                              fm = "minres",
                              plot = F,
                              quant = 0.95)



  #simulated F1:12 instead of F1:12.
  sim_factor <- fa_fits_out2$values[ ,(3*p+1):(3*p+1+p-1)]
  qntls <- apply(sim_factor, 2, function(x) quantile(x, probs = .95))
  means <- apply(sim_factor, 2, function(x) mean(x))




  pdf(paste0("figures/2factors_cor5_n",v_n[i],".pdf"),width=4.5,height=5)

  plot(qntls, xaxt="n",type="l",cex.axis = 1.2, col="red", cex.lab=1.2, lwd=2, ylab="eigen values",las=1,lty=2, ylim=c(0, 2),xlab="eigen value number", bty="n",main=" ")
  lines(fa_fits_out2$fa.values, type="b", lwd=2)
  lines(means, col="red", lwd=2)
  axis(side=1, at=seq(1,12),cex.axis = 1.2)
  if(i==3) legend(3,2,
                  legend = c("Actual data", "Mean", "95 quantile"),
                  col = c("black", "red","red"),
                  bty = "n",
                  lty = c(1,1,2),
                  pt.cex = 1,
                  cex = 1.2,
                  lwd = 2,
                  text.col = "black")
  dev.off()



}

# ---------- Figure 6 ---------------------------------------

v_cor <- c(.3, .4, .5)

for(i in 1:3) {

  corr <- v_cor[i]
  sig<-matrix(0,6,6)
  diag(sig)<-1
  sig[1,2]<-corr
  sig[2,1]<-corr
  sig[1,3]<-corr
  sig[3,1]<-corr
  sig[3,2]<-corr
  sig[2,3]<-corr
  sig[5,4]<-corr
  sig[4,5]<-corr
  sig[4,6]<-corr
  sig[6,4]<-corr
  sig[6,5]<-corr
  sig[5,6]<-corr
  sig
  p<-6
  data <- mvrnorm(n = 1000,
                  mu=rep(0,6),
                  Sigma=sig,
                  empirical=T)

  fa_fits_out2 <- fa.parallel(x = data,
                              n.iter = 20,
                              fa = "fa",
                              fm = "minres",
                              plot = F,
                              quant = 0.95)



  #simulated F1:12 instead of F1:12.
  sim_factor <- fa_fits_out2$values[ ,(3*p+1):(3*p+1+p-1)]
  qntls <- apply(sim_factor, 2, function(x) quantile(x, probs = .95))
  means <- apply(sim_factor, 2, function(x) mean(x))




  pdf(paste0("figures/2factors_cor",v_cor[i],"_n1000.pdf"),width=4.5,height=5)

  plot(qntls, xaxt="n",type="l",cex.axis = 1.2, col="red", cex.lab=1.2, lwd=2, ylab="eigen values",las=1,lty=2, ylim=c(0, 2),xlab="eigen value number", bty="n",main=" ")
  lines(fa_fits_out2$fa.values, type="b", lwd=2)
  lines(means, col="red", lwd=2)
  axis(side=1, at=seq(1,12),cex.axis = 1.2)
  if(i==3) legend(3,2,
                  legend = c("Actual data", "Mean", "95 quantile"),
                  col = c("black", "red","red"),
                  bty = "n",
                  lty = c(1,1,2),
                  pt.cex = 1,
                  cex = 1.2,
                  lwd = 2,
                  text.col = "black")
  dev.off()



}



# --------------------------------------------------------------
# ---------- Figure Appendix F (Polytomous responses) ----------
# --------------------------------------------------------------

# We investigate the performance drop when moving from continuous
# slowly to binary responses
# We pick a condition in which all methods are making some mistakes
# to get a good differentiation


nIter <- 200
v_cats <- 2:10

a_res <- array(NA, dim=c(nIter, 10, 3))

for(i in 1:nIter) {

  for(j in 1:10) {

    if(j<10) {
      cats_j <- v_cats[j]
    } else {
      cats_j <- NULL
    }

    set.seed(i)

    ## Generate Data
    dg_out <- f_datagen2(nfactor = 4,
                         items_pf = 6,
                         psi = .8,
                         n = 1000,
                         ra = c(0.3, 1),
                         cats = cats_j)

    data <- dg_out$data

    ## Estimate
    # PE
    k_PE_10 <- fspe(data = data,
                    maxK = 6,
                    nfold = 10,
                    rep = 10,
                    method = "PE",
                    pbar = FALSE)

    a_res[i, j, 1] <- k_PE_10$nfactor

    # EGA
    out_EGA <- EGA(data = data, plot.EGA = FALSE, verbose = FALSE)
    a_res[i, j, 2] <- out_EGA$n.dim

    # PA
    fa_fits_out2 <- fa.parallel(x = data,
                                n.iter = 20,
                                fa = "fa",
                                fm = "minres",
                                plot = FALSE,
                                quant = 0.95)

    a_res[i, j, 3] <- fa_fits_out2$nfact


  }
  print(i)
}

# saveRDS(a_res, file="files/PolytomousSim.RDS")

# Compute accuracy
a_acc <- a_res==4
m_acc <- apply(a_acc, 2:3, mean)

# -------- Make Figure for appendix ----------

sc <- 0.8
pdf("figures/App_Polytomous.pdf", width=7*sc, height=6*sc)

# Colors
cols <- brewer.pal(5, "Set1")

# Polytomous
plot.new()
plot.window(xlim=c(2, 11), ylim=c(0,1))
axis(1, labels=c(2:10, "Cont."), at=2:11, las=1)
axis(2, las=2)
lines(2:10, m_acc[-10, 1], col=cols[1], lty=1, lwd=2)
lines(2:10, m_acc[-10, 3], col=cols[2], lty=2, lwd=2)
lines(2:10, m_acc[-10, 2], col=cols[3], lty=3, lwd=2)
points(2:11, m_acc[, 1], col=cols[1], pch=20, lwd=2)
points(2:11, m_acc[, 3], col=cols[2], pch=20, lwd=2)
points(2:11, m_acc[, 2], col=cols[3], pch=20, lwd=2)

legend("bottomright", legend = c("PE", "Parallel", "EGA"),
       col=cols, lty=1:3, bty="n", lwd=rep(3,1))

title(xlab=c("Response categories"), line=2.5)
title(ylab=c("Accuracy"), line=2.5)

dev.off()



# --------------------------------------------------------------
# ---------- Assess Bias-Variance Trade-off [for reply letter] -
# --------------------------------------------------------------

# ----- True: 6-factor with varying psi -----

# Storage
nIter <- 25
n_seq <- seq(100, 1000, length=7)
n_n_seq <- length(n_seq)
psi_seq <- c(0, .8)
a_res <- array(NA, dim=c(nIter, n_n_seq, 2, 6))

# Loop
set.seed(1)
for(i in 1:nIter){
  print(i)
  for(psi in 1:2) {
    for(n in 1:n_n_seq) {

      print(paste0("i = ", i, "; psi = ", psi, "; n = ", n))

      # Generate Data
      dg_out <- f_datagen2(nfactor = 6,
                           items_pf = 6,
                           psi = psi_seq[psi],
                           n = n_seq[n],
                           ra = c(0.3, 1),
                           cats = NULL)

      # Fit factor models
      out <- fspe(dg_out$data, maxK = 6, pbar = FALSE)

      a_res[i, n, psi, ] <- out$PEs
    }
  }
}

saveRDS(a_res, "files/BiasVar_Sim.RDS")


# Aggregate
m_res <- apply(a_res, 2:4, mean)

m_ylim <- rbind(c(.8, 1), c(.76, .82))

# Plot Figure
cols <- brewer.pal(7, "Set1")[-6]

sc <- 1.1
pdf("figures/Fig_BiasVariance.pdf", width = 8*sc, height = 4.5*sc)

par(mfrow=c(1,2))
for(psi in 1:2) {

  # layout
  plot.new()
  plot.window(xlim=c(1,n_n_seq), ylim=m_ylim[psi, ])
  axis(1, labels=n_seq, at=1:n_n_seq, cex.axis=.9)
  axis(2, las=2)
  title(main=paste0("Psi = ", psi_seq[psi]), font.main=1, line=1.5)
  title(xlab="N")
  title(ylab="Prediction Error")

  # data
  for(i in 1:6) lines(m_res[, psi, i], col=cols[i])

  # legend
  if(psi==2) legend("topright", legend=paste0(1:6, " Factor"),
                    col=cols, lwd = rep(1,6), bty="n")

}

dev.off()















