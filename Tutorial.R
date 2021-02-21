# jonashaslbeck@gmail.com; February 18, 2020

# -----------------------------------------------------
# -------- Tutorial -----------------------------------
# -----------------------------------------------------

# library(devtools)
# install_github("jmbh/fspe")
library(fspe)

data(holzinger19)

set.seed(1)
fspe_out <- fspe(data = holzinger19,
                 maxK = 10,
                 nfold = 10,
                 rep = 1,
                 method = "PE")

fspe_out$nfactor # estimated factors = 4



# -----------------------------------------------------
# -------- Make Figure for paper ----------------------
# -----------------------------------------------------

pdf("figures/Tutorial_PE_path.pdf", width = 6, height = 4)
par(mar=c(4,4,1,1))
plot.new()
plot.window(xlim=c(1, 10), ylim=c(.6, .8))
axis(1, 1:10)
axis(2, las=2)
lines(fspe_out$PEs)
points(fspe_out$PEs, pch=20, cex=1.5)
title(xlab="Number of Factors", ylab="Prediction Error")
dev.off()
