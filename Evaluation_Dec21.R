# jonashaslbeck@gmail.com; Dec 16, 2021

# --------------------------------------------------------------
# ---------- Load Simulation Output ----------------------------
# --------------------------------------------------------------

simResDir <- "output/"

n_methods <- 18

v_files <- list.files(simResDir)
n_files <- length(v_files)
l_files <- list()
a_results <- array(NA, dim=c(n_files, 12, 4, 3, 3, n_methods))
for(i in 1:n_files) a_results[i, , , , , ] <- readRDS(paste0(simResDir, v_files[i]))

# --------------------------------------------------------------
# ---------- Preprocess ----------------------------------------
# --------------------------------------------------------------

# -------- Calculate Accuracy --------

# Set the conditions with psi=2,3 for nfactor = 1 to NA
a_results[, , 1, , 2:3, ] <- NA

# Get indicator: correct
v_nfactor <- c(1, 2, 4, 6) # true factors
a_indcorrect <- array(NA, dim=c(n_files, 12, 4, 3, 3, n_methods))
for(nf in 1:4)  a_indcorrect[, , nf, , , ] <- a_results[, , nf, , , ] == v_nfactor[nf]
# Average across iterations (compute accuracy)
a_acc <- apply(a_indcorrect, 2:6, mean)


# -------- Calculate Bias --------
a_bias <- apply(a_results, 2:6, function(x) mean(x, na.rm=T))


# -------- Sanity Checks --------

# a_results[, , , , , 9]
#
# dim(a_bias)
# a_results[, 12, 3, 2, 2, n_methods]
# hist(a_results[, 12, 3, 2, 2, n_methods])



# --------------------------------------------------------------
# ---------- Overview & Sanity Checks --------------------------
# --------------------------------------------------------------

# ----- Get Labels/ordering from Simulation study ------

labels <- c("PE", "CovE", "OoS Cov 2", "Parallel",
            "EGA", "K-G",
            "MAP_ml", "BIC_ml",  "VSS_ml",
            "RMSEA_ml", "RMSEA_CI_ml", "AIC_ml",
            "MAP_minres", "BIC_minres",  "VSS_minres",
            "RMSEA_minres", "RMSEA_CI_minres", "AIC_minres")

# ----- Get marginal performance -------
marg_acc <- apply(a_acc, 5, function(x) mean(x, na.rm = TRUE))
names(marg_acc) <- labels
ord <- order(marg_acc, decreasing = TRUE)

par(mar=c(8,4,2,1))
round(marg_acc, 3)
round(marg_acc[ord], 2)
barplot(marg_acc[ord], las=2, ylim=c(0, 1))
abline(h=max(marg_acc), lty=2)


# --------------------------------------------------------------
# ---------- Subset Methods ------------------------------------
# --------------------------------------------------------------

# Delete:
# (1) all "minres" methods
# (2) the worse nfold for PE and CovE
# (3) the worse version of RMSEA

methods_keep <- c(1,2,4,5,6,7,8,9,10,12)
length(methods_keep)
labels[methods_keep]

a_acc_ss <- a_acc[, , , , methods_keep]
labels_ss <- labels[methods_keep]
labels_ss[6:10] <- c("MAP", "BIC", "VSS", "RMSEA", "AIC")

# --------------------------------------------------------------
# ---------- Stability of those results ------------------------
# --------------------------------------------------------------

# Storage & batching
batch_end <- round(seq(0, 200, length=21))
m_batchresults <- matrix(NA, nrow=20, ncol=10)

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
n_methods_ss <- 10
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


pdf("figures/Figure1_Marginals_byN.pdf", width = 8.5, height = 3)

n_methods_ss <- 10

# ----- Make Design -----
lmat <- matrix(1:2, 1, 2)
lo <- layout(mat = lmat, widths = c(1, .2))

# ----- Plot canvas 1 -----
par(mar=c(2, 4, 0.5, .5))
plot.new()
plot.window(xlim=c(0.65,n_methods_ss), ylim=c(0,1))
axis(2, las=2)
axis(1, labels=labels_ss[ord], at=1:n_methods_ss, tick = FALSE, cex.axis=0.95, line = -.75)
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

  # points(1:9, means, pch=20, cex=1.5)
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



