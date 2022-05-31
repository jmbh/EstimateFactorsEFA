# jonashaslbeck@gmail.com; May 31, 2022

# --------------------------------------------------------------
# ---------- Load Packages -------------------------------------
# --------------------------------------------------------------

library(RColorBrewer)

# --------------------------------------------------------------
# ---------- Load Simulation Output ----------------------------
# --------------------------------------------------------------

simResDir <- "output/"

n_methods <- 21

v_files <- list.files(simResDir)
n_files <- length(v_files)
l_files <- list()
a_results <- array(NA, dim=c(n_files, 12, 4, 3, 3, n_methods))
for(i in 1:n_files) a_results[i, , , , , ] <- readRDS(paste0(simResDir, v_files[i]))


# --------------------------------------------------------------
# ---------- Preprocess ----------------------------------------
# --------------------------------------------------------------

# -------- Make Composite Prediction --------
# This is used for an Appendix, in which we show the performance of
# combining the predictions of various methods with a majority rule

n_methods_aug <- n_methods + 8
a_results_aug <- array(NA, dim=c(n_files, 12, 4, 3, 3, n_methods_aug))
a_results_aug[, , , , , 1:n_methods] <- a_results

# using top 3: PE, Parallel, EGA
nIter <- n_files
for(i in 1:nIter) { # iterations
  for(n in 1:12) { # n-var
    for(k in 1:4) { # factors
      for(j in 1:3) { # items
        for(s in 1:3) { # correlation

          # The predictions of  the top 10 performing methods ordered
          preds_top10 <- a_results[i, n, k, j, s, c(1,4,5,18,2,8,7,6,9,10)]

          # Create predictions for Top 3,4,...,10
          for(q in 1:8) {

            # MAJORITY RULE
            pred_q <- preds_top10[1:(2+q)]
            if(sum(is.na(pred_q)) != (2+q)) {
              # Majority rule; if draw: random across ties
              tb <- table(pred_q)
              set_max <- as.numeric(names(which(tb == max(tb) )))
              final_pred <- set_max[sample(1:length(set_max), size=1)] # if there is only 1, no randomness

              a_results_aug[i, n, k, j, s, n_methods + q] <- final_pred
            }

          } # for q: majority of Top 2+q


        } # for: psi
      }
    }
  }
} # for: iter


# -------- Calculate Accuracy [with Composite] --------

# Set the conditions with psi=2,3 for nfactor = 1 to NA
a_results_aug[, , 1, , 2:3, ] <- NA

# Get indicator: correct
v_nfactor <- c(1, 2, 4, 6) # true factors
a_indcorrect <- array(NA, dim=c(n_files, 12, 4, 3, 3, n_methods_aug))
for(nf in 1:4)  a_indcorrect[, , nf, , , ] <- a_results_aug[, , nf, , , ] == v_nfactor[nf]
# Average across iterations (compute accuracy)
a_acc_comp <- apply(a_indcorrect, 2:6, mean)

labels_comp <-  c("PE nfold=10", "CovE nfold=10", "OoS Cov nfold=2", "Parallel",
                  "EGA", "K-G",
                  "MAP_ml", "BIC_ml",  "VSS_ml",
                  "RMSEA_ml", "RMSEA_CI_ml", "AIC_ml",
                  "MAP_minres", "BIC_minres",  "VSS_minres",
                  "RMSEA_minres", "RMSEA_CI_minres", "AIC_minres",
                  "PE nfold=2", "PE nfold10, rep10", "CovE nfold10, rep10", paste0("Top ", 3:10))

# ----- Get marginal performance -------
marg_acc_comp <- apply(a_acc_comp, 5, function(x) mean(x, na.rm = TRUE))
names(marg_acc_comp) <- labels
ord <- order(marg_acc_comp, decreasing = TRUE)

par(mar=c(8,4,2,1))
round(marg_acc_comp[ord], 3)
barplot(marg_acc_comp[ord], las=2, ylim=c(0, 1))
abline(h=marg_acc_comp[20], lty=2)



# -------- Calculate Accuracy [without composite] --------
# These are the results we report in the main text

# Set the conditions with psi=2,3 for nfactor = 1 to NA
a_results[, , 1, , 2:3, ] <- NA

# Get indicator: correct
v_nfactor <- c(1, 2, 4, 6) # true factors
a_indcorrect <- array(NA, dim=c(n_files, 12, 4, 3, 3, n_methods))
for(nf in 1:4)  a_indcorrect[, , nf, , , ] <- a_results[, , nf, , , ] == v_nfactor[nf]
# Average across iterations (compute accuracy)
a_acc <- apply(a_indcorrect, 2:6, mean)

dim(a_acc)

a_acc[6, 3, 1, 2, 1]  # n-var, nfactor, nitems, Psi, methods [1=PE]


# -------- Calculate Bias --------

# Note this is not the actual bias, but the mean; I subtract below the corresponding numnber of true factors
a_bias <- apply(a_results, 2:6, function(x) mean(x, na.rm=T))

# --------------------------------------------------------------
# ---------- Overview & Sanity Checks --------------------------
# --------------------------------------------------------------

# ----- Get Labels/ordering from Simulation study ------

labels <- c("PE nfold=10", "CovE nfold=10", "OoS Cov nfold=2", "Parallel",
            "EGA", "K-G",
            "MAP_ml", "BIC_ml",  "VSS_ml",
            "RMSEA_ml", "RMSEA_CI_ml", "AIC_ml",
            "MAP_minres", "BIC_minres",  "VSS_minres",
            "RMSEA_minres", "RMSEA_CI_minres", "AIC_minres",
            "PE nfold=2", "PE nfold10, rep10", "CovE nfold10, rep10")

# ----- Get marginal performance -------
marg_acc <- apply(a_acc, 5, function(x) mean(x, na.rm = TRUE))
names(marg_acc) <- labels
ord <- order(marg_acc, decreasing = TRUE)

par(mar=c(8,4,2,1))
round(marg_acc, 3)
round(marg_acc[ord], 2)
BP <- barplot(marg_acc[ord][-c(1,2)], las=2, ylim=c(0, 1))
abline(h=max(marg_acc[ord][-c(1,2)]), lty=2)

vals <- marg_acc[ord][-c(1,2)]
text(BP, vals + .1, labels = round(vals, 3), srt=90)



# --------------------------------------------------------------
# ---------- Subset Methods ------------------------------------
# --------------------------------------------------------------

# Delete:
# (1) all "minres" methods
# (2) the worse nfold for PE and CovE
# (3) the worse version of RMSEA

methods_keep <- c(20,21,4,5,6,7,8,9,10,11,12)
length(methods_keep)
labels[methods_keep]

a_acc_ss <- a_acc[, , , , methods_keep]
labels_ss <- labels[methods_keep]
labels_ss[1:2] <- c("PE", "CovE")
labels_ss[6:11] <- c("MAP", "BIC", "VSS", "RMSEA", "RMSEA(CI)", "AIC")


# --------------------------------------------------------------
# ---------- Stability of those results ------------------------
# --------------------------------------------------------------

# Storage & batching
batch_end <- round(seq(0, n_files, length=21))
m_batchresults <- matrix(NA, nrow=20, ncol=11)

# Subset to Methods in Figure 1
a_indcorrect_ss <- a_indcorrect[, , , , , methods_keep]

# Make new ordering on full-data iteration data subset
marg_acc_ss <- apply(a_acc_ss, 5, function(x) mean(x, na.rm = TRUE))
names(marg_acc_ss) <- labels_ss
ord_ss <- order(marg_acc_ss, decreasing = TRUE)

# Loop over batches
for(i in 1:20) {
  a_acc_batch <- apply(a_indcorrect_ss[(batch_end[i]+1):batch_end[i+1], , , , , ], 2:6, mean)
  marg_acc_ss <- apply(a_acc_batch, 5, function(x) mean(x, na.rm = TRUE))
  names(marg_acc_ss) <- labels_ss
  m_batchresults[i, ] <- round(marg_acc_ss[ord_ss], 3)
}

m_batchresults

m_batchresults[, 1] > m_batchresults[, 2]
m_batchresults[, 2] > m_batchresults[, 3]
m_batchresults[, 3] > m_batchresults[, 4]
m_batchresults[, 4] > m_batchresults[, 5]
m_batchresults[, 5] > m_batchresults[, 6]
m_batchresults[, 6] > m_batchresults[, 7]
m_batchresults[, 7] > m_batchresults[, 8]
m_batchresults[, 8] > m_batchresults[, 9]
m_batchresults[, 9] > m_batchresults[, 10] # only the last two are exchanged 3/10



# --------------------------------------------------------------
# ---------- EGA Undecided -------------------------------------
# --------------------------------------------------------------

tb <- table(a_results[, , , , , 5] == 0)

tb / sum(tb)

# --------------------------------------------------------------
# ---------- Figure 1b: Marginal performance [split by N] ------
# --------------------------------------------------------------

# ----- Calculate Aggregates -----

# Overall mean and ordering
n_methods_ss <- 11
m_acc_cp <- matrix(NA, 432, n_methods_ss)
for(i in 1:n_methods_ss) m_acc_cp[, i] <-  as.numeric(a_acc_ss[, , , , i])
ord <- order(colMeans(m_acc_cp,  na.rm=T), decreasing=T)
colnames(m_acc_cp) <- labels_ss
m_acc_cp_ord <- m_acc_cp[, ord]

# Get numbers for best methods ordered by peformance
sort(colMeans(m_acc_cp,  na.rm=T), decreasing=T)

# N=100
m_acc_cp_100 <- matrix(NA, 432/12, n_methods_ss)
for(i in 1:n_methods_ss) m_acc_cp_100[, i] <-  as.numeric(a_acc_ss[1, , , , i])
m_acc_cp_100_ord <- m_acc_cp_100[, ord] # order based on grand-average

# N=847
m_acc_cp_847 <- matrix(NA, 432/12, n_methods_ss)
for(i in 1:n_methods_ss) m_acc_cp_847[, i] <-  as.numeric(a_acc_ss[7, , , , i])
m_acc_cp_847_ord <- m_acc_cp_847[, ord] # order based on grand-average

# N=5000
m_acc_cp_5000 <- matrix(NA, 432/12, n_methods_ss)
for(i in 1:n_methods_ss) m_acc_cp_5000[, i] <-  as.numeric(a_acc_ss[12, , , , i])
m_acc_cp_5000_ord <- m_acc_cp_5000[, ord] # order based on grand-average

# Make list
m_acc_PLOT <- list(m_acc_cp_100_ord,
                   m_acc_cp_847_ord,
                   m_acc_cp_5000_ord)


sc <- 1.1
pdf("figures/Figure1_Marginals_byN_May25.pdf", width = 8.5*sc, height = 3*sc)

# ----- Make Design -----
lmat <- matrix(1:2, 1, 2)
lo <- layout(mat = lmat, widths = c(1, .2))

# ----- Plot canvas 1 -----
par(mar=c(2, 4, 0.5, .5))
plot.new()
plot.window(xlim=c(0.65,n_methods_ss), ylim=c(0,1))
axis(2, las=2)
axis(1, labels=labels_ss[ord], at=1:n_methods_ss,
     tick = FALSE, cex.axis=0.70, line = -.75)
title(ylab = "Accuracy", cex.lab=1.1)


# ----- Some specs -----
sshift <- 0.25
ss <- c(-sshift, 0, sshift)
cols <- c("lightgrey", "darkgrey", "black")

# ----- Plot data -----

for(i in 1:3) {

  # Caclulate means and quantiles
  m_acc_PLOT_i <- m_acc_PLOT[[i]]
  means <- apply(m_acc_PLOT_i, 2, function(x) mean(x, na.rm=TRUE))
  qu <- apply(m_acc_PLOT_i, 2, function(x) quantile(x, na.rm=TRUE, prob=c(0.10, 0.90)))

  # Plot Means
  wid <- 0.09
  bar <- 0.0025
  rect((1:n_methods_ss)-wid+ss[i], means-bar, (1:n_methods_ss)+wid+ss[i], means+bar,
       col=cols[i],
       border=cols[i])

  lwd <- 1.5
  segments(1:n_methods_ss+ss[i], qu[1,], 1:n_methods_ss+ss[i], qu[2,], lwd=lwd, col=cols[i])
  wisk <- wid #0.075
  segments((1:n_methods_ss)-wisk+ss[i], qu[1,], (1:n_methods_ss)+wisk+ss[i], qu[1,], lwd=lwd, col=cols[i])
  segments((1:n_methods_ss)-wisk+ss[i], qu[2,], (1:n_methods_ss)+wisk+ss[i], qu[2,], lwd=lwd, col=cols[i])

} # end for: i

# ----- Add Legend -----

par(mar=c(2, 1, 0.5, .5))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))

# Legend 1
rect(0.2-wid, .9-bar ,0.2+wid, .9+bar,
     col=cols[1],
     border=cols[1])
text(0.35, .9, "Mean", adj=0, cex=0.8)

fakequ <- .05
segments(0.2, .75-fakequ, 0.2, .75+fakequ, lwd=lwd, col=cols[1])
wisk <- wid #0.075
segments(.2-wisk, .75-fakequ, .2+wisk, .75-fakequ, lwd=lwd, col=cols[1])
segments(.2-wisk, .75+fakequ, .2+wisk, .75+fakequ, lwd=lwd, col=cols[1])

text(0.35, .79, "10 & 90%", adj=0, cex=0.8)
text(0.35, .71, "Quantiles", adj=0, cex=0.8)


# Legend 2
legend(-.2, .6, c("N = 100", "N = 847", "N = 5000"),
       bty="n",
       text.col=cols,
       cex=.95)

dev.off()




# --------------------------------------------------------------
# ---------- Copy Simulation Setup -----------------------------
# --------------------------------------------------------------

# Specifications: Data generation
v_nfactor <- c(1, 2, 4, 6)
v_items_pf <- c(3, 6, 12)
v_psi <- c(0, .4, .8)
v_n <- round(exp(seq(4.61, 8.5171, length=12)))
v_n


# --------------------------------------------------------------
# ---------- Figure 2: Performance per Cell --------------------
# --------------------------------------------------------------

# meth_select <- c(1,2,4,5,8)

meth_select <- c(1,4,5,2,12) # present methods in order of Figure 1
# Color scheme
library(RColorBrewer)
cols <- brewer.pal(5, "Set1")

sc <- .05
pdf("figures/Figure2_PerfScenarios.pdf", width = 210*sc, height=297*sc)

# ----- Define Layout -------

lmat <- matrix(NA, nrow=13, ncol=6)
# top line
lmat[1, ] <- c(0,0,13:16)

# block 1
lmat[2, ] <- c(1,4,17:20)
lmat[3, ] <- c(1,5,21:24)
lmat[4, ] <- c(1,6,25:28)
lmat[5, ] <- rep(0, 6)

# block 2
lmat[6, ] <- c(2,7,29:32)
lmat[7, ] <- c(2,8,33:36)
lmat[8, ] <- c(2,9,37:40)
lmat[9, ] <- rep(0, 6)

# block 3
lmat[10, ] <- c(3,10,41:44)
lmat[11, ] <- c(3,11,45:48)
lmat[12, ] <- c(3,12,49:52)
lmat[13, ] <- rep(0, 6)

# Legend
# lmat[14, ] <- rep(53, 6)

lo <- layout(lmat,
             width = c(.2, .1, 1, 1, 1, 1),
             height = c(.15,
                        1,1,1,.1,
                        1,1,1,.1,
                        1,1,1,.1))
# layout.show(lo)

# ----- Plot all labels -------

# Label plot function
plotLab <- function(label, rot, cex) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
  text(0, 0, label, srt=rot, cex=cex)
}


# Side labels big
plotLab(label="Variables per Factor = 3", rot=90, cex=2)
plotLab(label="Variables per Factor = 6", rot=90, cex=2)
plotLab(label="Variables per Factor = 12", rot=90, cex=2)

# Side labels small
for(i in 1:3) {
  plotLab(label=expression(paste(Psi, " = 0")), rot=90, cex=2)
  plotLab(label=expression(paste(Psi, " = 0.40")), rot=90, cex=2)
  plotLab(label=expression(paste(Psi, " = 0.80")), rot=90, cex=2)
}

# Top labels
plotLab(label="  Factors = 1", rot=0, cex=2)
plotLab(label="  Factors = 2", rot=0, cex=2)
plotLab(label="  Factors = 4", rot=0, cex=2)
plotLab(label="  Factors = 6", rot=0, cex=2)


# ----- Plot data -------

for(it in 1:3) {
  for(ps in 1:3) {
    for(nf in 1:4) {

      # Do not show: nfactors 1 for item/fac=3; and: variations of psi for nfac=1
      data_cell <- !any(c((it==1 & nf == 1),
                          nf==1 & ps > 1))

      # Plot data when available
      if(data_cell) {

        # Plotting
        par(mar=c(2,3,1,0))
        plot.new()
        plot.window(xlim=c(1,12), ylim=c(0,1))

        # Conditional y-axes
        if(nf==1) {
          axis(2, las=2)
        } else if(nf==2 & !((it==2 & ps==1) | (it==3 & ps==1))) {
          axis(2, las=2)
        } else {
          axis(2, labels=rep("", 6), at=seq(0, 1, length=6), las=2)
        }

        # Conditional x-axes
        if(ps == 3) {
          axis(1, labels = v_n, at=1:12, las=2, cex=.5)
        } else {
          axis(1, labels = rep("", 12), at=1:12, las=2)
        }

        # Plot data
        for(meth in 1:length(meth_select)) {
          lines(a_acc[, nf, it, ps, meth_select[meth]],
                col=cols[meth],
                lwd=2,
                lty=meth)
        }

      } else {

        plot.new()
        plot.window(xlim = c(-1, 1), ylim=c(-1, 1))

        if(it==1 & ps == 2)

          legend("center", text.col=cols, lty=1:5, col=cols,
                 legend = c("PE", "Parallel", "EGA", "CovE", "AIC"),
                 ncol=1, bty="n", cex=1.8, lwd=rep(2, 5))

      }
    }
  }
} # end for: it


dev.off()



# --------------------------------------------------------------
# ---------- Figure 2b: Bias per Cell --------------------------
# --------------------------------------------------------------

meth_select <- c(1,4,5,2,12) # present methods in order of Figure 1

# Color scheme
library(RColorBrewer)
cols <- brewer.pal(5, "Set1")

sc <- .05
pdf("figures/Figure_AppE.pdf", width = 210*sc, height=297*sc)

# ----- Define Layout -------

lmat <- matrix(NA, nrow=13, ncol=6)
# top line
lmat[1, ] <- c(0,0,13:16)

# block 1
lmat[2, ] <- c(1,4,17:20)
lmat[3, ] <- c(1,5,21:24)
lmat[4, ] <- c(1,6,25:28)
lmat[5, ] <- rep(0, 6)

# block 2
lmat[6, ] <- c(2,7,29:32)
lmat[7, ] <- c(2,8,33:36)
lmat[8, ] <- c(2,9,37:40)
lmat[9, ] <- rep(0, 6)

# block 3
lmat[10, ] <- c(3,10,41:44)
lmat[11, ] <- c(3,11,45:48)
lmat[12, ] <- c(3,12,49:52)
lmat[13, ] <- rep(0, 6)

lo <- layout(lmat,
             width = c(.2, .1, 1, 1, 1, 1),
             height = c(.15,
                        1,1,1,.1,
                        1,1,1,.1,
                        1,1,1,.1))
# layout.show(lo)

# ----- Plot all labels -------

# Label plot function
plotLab <- function(label, rot, cex) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
  text(0, 0, label, srt=rot, cex=cex)
}


# Side labels big
plotLab(label="Variables per Factor = 3", rot=90, cex=2)
plotLab(label="Variables per Factor = 6", rot=90, cex=2)
plotLab(label="Variables per Factor = 12", rot=90, cex=2)

# Side labels small
for(i in 1:3) {
  plotLab(label=expression(paste(Psi, " = 0")), rot=90, cex=2)
  plotLab(label=expression(paste(Psi, " = 0.40")), rot=90, cex=2)
  plotLab(label=expression(paste(Psi, " = 0.80")), rot=90, cex=2)
}

# Top labels
plotLab(label="  Factors = 1", rot=0, cex=2)
plotLab(label="  Factors = 2", rot=0, cex=2)
plotLab(label="  Factors = 4", rot=0, cex=2)
plotLab(label="  Factors = 6", rot=0, cex=2)


# ----- Plot data -------

for(it in 1:3) {
  for(ps in 1:3) {
    for(nf in 1:4) {

      # Do not show: nfactors 1 for item/fac=3; and: variations of psi for nfac=1
      data_cell <- !any(c((it==1 & nf == 1),
                          nf==1 & ps > 1))

      # Plot data when available
      if(data_cell) {

        # Plotting
        par(mar=c(2,3,1,0))
        plot.new()
        plot.window(xlim=c(1,12), ylim=c(-6,6))

        # Conditional y-axes
        if(nf==1) {
          axis(2, las=2)
        } else if(nf==2 & !((it==2 & ps==1) | (it==3 & ps==1))) {
          axis(2, las=2)
        } else {
          axis(2, labels=rep("", 7), at=seq(-6, 6, length=7), las=2)
        }


        # Conditional x-axes
        if(ps == 3) {
          axis(1, labels = v_n, at=1:12, las=2, cex=.5)
        } else {
          axis(1, labels = rep("", 12), at=1:12, las=2)
        }

        # Correct line
        abline(h=0, col="grey")

        # Plot data
        for(meth in 1:length(meth_select)) {
          lines(a_bias[, nf, it, ps, meth_select[meth]] - c(1,2,4,6)[nf],
                col=cols[meth],
                lwd=2,
                lty=meth)
        }

      } else {

        plot.new()
        plot.window(xlim = c(-1, 1), ylim=c(-1, 1))

        if(it==1 & ps == 2)

          legend("center", text.col=cols, lty=1:5, col=cols,
                 legend = c("PE", "Parallel", "EGA", "CovE", "AIC"),
                 ncol=1, bty="n", cex=1.8, lwd=rep(2, 5))

      }
    }
  }
} # end for: it


dev.off()


# --------------------------------------------------------------
# ---------- Appendidix: Bias per Cell (incl. BIC, RMSEA) ------
# --------------------------------------------------------------

meth_select <- c(1,4,5,2,12, 10, 11, 8) # present methods in order of Figure 1

# Color scheme
cols <- brewer.pal(9, "Set1")[-6]

sc <- .05
pdf("figures/Figure_AppE_addMethods.pdf", width = 210*sc, height=297*sc)

# ----- Define Layout -------

lmat <- matrix(NA, nrow=13, ncol=6)
# top line
lmat[1, ] <- c(0,0,13:16)

# block 1
lmat[2, ] <- c(1,4,17:20)
lmat[3, ] <- c(1,5,21:24)
lmat[4, ] <- c(1,6,25:28)
lmat[5, ] <- rep(0, 6)

# block 2
lmat[6, ] <- c(2,7,29:32)
lmat[7, ] <- c(2,8,33:36)
lmat[8, ] <- c(2,9,37:40)
lmat[9, ] <- rep(0, 6)

# block 3
lmat[10, ] <- c(3,10,41:44)
lmat[11, ] <- c(3,11,45:48)
lmat[12, ] <- c(3,12,49:52)
lmat[13, ] <- rep(0, 6)

lo <- layout(lmat,
             width = c(.2, .1, 1, 1, 1, 1),
             height = c(.15,
                        1,1,1,.1,
                        1,1,1,.1,
                        1,1,1,.1))
# layout.show(lo)

# ----- Plot all labels -------

# Label plot function
plotLab <- function(label, rot, cex) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
  text(0, 0, label, srt=rot, cex=cex)
}


# Side labels big
plotLab(label="Variables per Factor = 3", rot=90, cex=2)
plotLab(label="Variables per Factor = 6", rot=90, cex=2)
plotLab(label="Variables per Factor = 12", rot=90, cex=2)

# Side labels small
for(i in 1:3) {
  plotLab(label=expression(paste(Psi, " = 0")), rot=90, cex=2)
  plotLab(label=expression(paste(Psi, " = 0.40")), rot=90, cex=2)
  plotLab(label=expression(paste(Psi, " = 0.80")), rot=90, cex=2)
}

# Top labels
plotLab(label="  Factors = 1", rot=0, cex=2)
plotLab(label="  Factors = 2", rot=0, cex=2)
plotLab(label="  Factors = 4", rot=0, cex=2)
plotLab(label="  Factors = 6", rot=0, cex=2)


# ----- Plot data -------

for(it in 1:3) {
  for(ps in 1:3) {
    for(nf in 1:4) {

      # Do not show: nfactors 1 for item/fac=3; and: variations of psi for nfac=1
      data_cell <- !any(c((it==1 & nf == 1),
                          nf==1 & ps > 1))

      # Plot data when available
      if(data_cell) {

        # Plotting
        par(mar=c(2,3,1,0))
        plot.new()
        plot.window(xlim=c(1,12), ylim=c(-6,6))

        # Conditional y-axes
        if(nf==1) {
          axis(2, las=2)
        } else if(nf==2 & !((it==2 & ps==1) | (it==3 & ps==1))) {
          axis(2, las=2)
        } else {
          axis(2, labels=rep("", 7), at=seq(-6, 6, length=7), las=2)
        }


        # Conditional x-axes
        if(ps == 3) {
          axis(1, labels = v_n, at=1:12, las=2, cex=.5)
        } else {
          axis(1, labels = rep("", 12), at=1:12, las=2)
        }

        # Correct line
        abline(h=0, col="grey")

        # Plot data
        for(meth in 1:length(meth_select)) {
          lines(a_bias[, nf, it, ps, meth_select[meth]] - c(1,2,4,6)[nf],
                col=cols[meth],
                lwd=2,
                lty=meth)
        }

      } else {

        plot.new()
        plot.window(xlim = c(-1, 1), ylim=c(-1, 1))

        if(it==1 & ps == 2)

          legend("center", text.col=cols, lty=1:8, col=cols,
                 legend = c("PE", "Parallel", "EGA", "CovE", "AIC", "RMSEA", "RMSEA_CI", "BIC"),
                 ncol=1, bty="n", cex=1.1, lwd=rep(2, 8))

      }
    }
  }
} # end for: it


dev.off()





# --------------------------------------------------------------
# ---------- Barplot: Average performance incl. Composite ------
# --------------------------------------------------------------


marg_acc_comp
labels_comp

# ----- Choose Methods Chosen in Figure 1 + Top1,...,10 -------
# Delete:
# (1) all "minres" methods
# (2) the worse nfold for PE and CovE
# (3) the worse version of RMSEA

methods_keep <- c(20,21,4,5,6,7,8,9,10,11,12, 22:29)
labels_comp_ss <- labels_comp[methods_keep]
labels_comp_ss[1:2] <- c("PE", "CovE")
labels_comp_ss[6:11] <- c("MAP", "BIC", "VSS", "RMSEA", "RMSEA(CI)", "AIC")
marg_acc_ss <- marg_acc_comp[methods_keep]
names(marg_acc_ss) <- labels_comp_ss
marg_acc_ss_ord <- sort(marg_acc_ss, decreasing = TRUE)

sc <- 1.1
pdf("figures/Fig_App_CompMeasures.pdf", width=6*sc, height=4.5*sc)
par(mar=c(6,4,2,1))
BP <- barplot(marg_acc_ss_ord, las=2, ylim=c(0, 1), ylab="Accuracy")
abline(h=marg_acc_ss_ord[6], lty=2)

vals <- marg_acc_ss_ord
text(BP, vals + .1, labels = round(vals, 3), srt=90)

dev.off()






