# jonashaslbeck@gmail.com; February 21, 2021

# --------------------------------------------------------------
# ---------- Load Simulation Output ----------------------------
# --------------------------------------------------------------

simResDir <- "output/"

v_files <- list.files(simResDir)
n_files <- length(v_files)
l_files <- list()
a_results <- array(NA, dim=c(n_files, 12, 4, 3, 3, 10))
for(i in 1:n_files) a_results[i, , , , , ] <- readRDS(paste0(simResDir, v_files[i]))


# --------------------------------------------------------------
# ---------- Preprocess ----------------------------------------
# --------------------------------------------------------------

# -------- Calculate Accuracy --------

# Set the conditions with psi=2,3 for nfactor = 1 to NA
a_results[, , 1, , 2:3, ] <- NA

# Get indicator: correct
v_nfactor <- c(1, 2, 4, 6) # true factors
a_indcorrect <- array(NA, dim=c(n_files, 12, 4, 3, 3, 10))
for(nf in 1:4)  a_indcorrect[, , nf, , , ] <- a_results[, , nf, , , ] == v_nfactor[nf]
# Average across iterations (compute accuracy)
a_acc <- apply(a_indcorrect, 2:6, mean)


# -------- Calculate Bias --------
a_bias <- apply(a_results, 2:6, function(x) mean(x, na.rm=T))



# --------------------------------------------------------------
# ---------- Overview & Sanity Checks --------------------------
# --------------------------------------------------------------

# ----- Get marginal performance -------
marg_acc <- apply(a_acc, 5, function(x) mean(x, na.rm = TRUE))
ord <- order(marg_acc, decreasing = TRUE)
labels <- c("PE", "CovE", "OoS Cov 2", "Parallel", "EGA", "K-G", "MAP", "BIC",  "VSS", "RMSEA")
names(marg_acc) <- labels

round(marg_acc, 3)
round(marg_acc[ord], 2)
barplot(marg_acc[ord], las=2, ylim=c(0, 1))
abline(h=max(marg_acc), lty=2)


# --------------------------------------------------------------
# ---------- Stability of those results ------------------------
# --------------------------------------------------------------

batch_end <- seq(0, 200, length=21)

m_batchresults <- matrix(NA, nrow=20, ncol=10)

for(i in 1:20) {

  a_acc_batch <- apply(a_indcorrect[(batch_end[i]+1):batch_end[i+1], , , , , ], 2:6, mean)
  marg_acc <- apply(a_acc_batch, 5, function(x) mean(x, na.rm = TRUE))
  ord <- order(marg_acc, decreasing = TRUE)
  names(marg_acc) <- labels
  m_batchresults[i, ] <- round(marg_acc[ord], 3)

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
m_batchresults[, 9] > m_batchresults[, 10]



# --------------------------------------------------------------
# ---------- EGA Undecided -------------------------------------
# --------------------------------------------------------------

tb <- table(a_results[, , , , , 5] == 0)

tb / sum(tb)


# --------------------------------------------------------------
# ---------- Figure 1: Marginal perforance ---------------------
# --------------------------------------------------------------

a_acc_m3 <- a_acc[, , , , -3]
labels2 <- labels[-3]
m_acc_cp <- matrix(NA, 432, 9)
for(i in 1:9) m_acc_cp[, i] <-  as.numeric(a_acc_m3[, , , , i])
ord <- order(colMeans(m_acc_cp,  na.rm=T), decreasing=T)
colnames(m_acc_cp) <- labels2
m_acc_cp_ord <- m_acc_cp[,ord]


pdf("figures/Figure1_Marginals.pdf", width = 7, height = 3)
# yet another alternative

par(mar=c(3, 4, 0.5, 1))
plot.new()
plot.window(xlim=c(0.65,9), ylim=c(0,1))
axis(2, las=2)
axis(1, labels=colnames(m_acc_cp_ord), at=1:9, tick = FALSE)
title(ylab="Accuracy", cex.lab=1.2)
# means + qu
means <- apply(m_acc_cp_ord, 2, function(x) mean(x, na.rm=TRUE))
qu <- apply(m_acc_cp_ord, 2, function(x) quantile(x, na.rm=TRUE, prob=c(0.10, 0.90)))

wid <- 0.3
# rect((1:9)-wid, rep(0, 9), (1:9)+wid, means, col="lightgrey")
bar <- 0.005
rect((1:9)-wid, means-bar, (1:9)+wid, means+bar, col="black")
# points(means, pch=20, cex=2)

# points(1:9, means, pch=20, cex=1.5)
lwd <- 1.5
segments(1:9, qu[1,], 1:9, qu[2,], lwd=lwd)
wisk <- 0.075
segments((1:9)-wisk, qu[1,], (1:9)+wisk, qu[1,], lwd=lwd)
segments((1:9)-wisk, qu[2,], (1:9)+wisk, qu[2,], lwd=lwd)

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

meth_select <- c(1,4,5,2,8) # present methods in order of Figure 1
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
                 legend = c("PE", "Parallel", "EGA", "CovE", "BIC"),
                 ncol=1, bty="n", cex=1.8, lwd=rep(2, 5))

      }
    }
  }
} # end for: it


dev.off()



# --------------------------------------------------------------
# ---------- Figure 2b: Bias per Cell --------------------------
# --------------------------------------------------------------

meth_select <- c(1,4,5,2,8) # present methods in order of Figure 1

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
                 legend = c("PE", "Parallel", "EGA", "CovE", "BIC"),
                 ncol=1, bty="n", cex=1.8, lwd=rep(2, 5))

      }
    }
  }
} # end for: it


dev.off()



