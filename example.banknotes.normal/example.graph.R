library(mclust)
library(PPcovMcd)

setwd("???")

asymptotic.quantiles.flag = FALSE

if (!asymptotic.quantiles.flag) {
  load("quantiles.rew.RData")
}

####

data("banknote")
key = "banknotes"

x = as.matrix(banknote[, -1])

####

n.in  = 100

I.in  = 1:n.in
I.out = (n.in + 1):(2*n.in)

dataset = x[c(I.in, I.out), ]

p     = ncol(dataset)
n.max = nrow(dataset)

####

purity.ppmcd   = rep(0.0, n.in)
purity.fastmcd = rep(0.0, n.in)
purity.detmcd  = rep(0.0, n.in)

TP.ppmcd   = rep(0.0, n.in)
TP.fastmcd = rep(0.0, n.in)
TP.detmcd  = rep(0.0, n.in)

FP.ppmcd   = rep(0.0, n.in)
FP.fastmcd = rep(0.0, n.in)
FP.detmcd  = rep(0.0, n.in)

for (n in (n.in + 1):n.max) {
  set.seed(n - n.in)
  
  cat("Status: ", 100*(n - n.in)/n.in, "%...\n", sep = "")

  x = as.matrix(dataset)
  rownames(x) = NULL
  colnames(x) = NULL

  p = ncol(x)

  x = x[1:n, ]

  mah.std = mah.standard(x)

  wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))

  ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            wgtFUN = wgtFUN)
  fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             wgtFUN = wgtFUN)
  detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", wgtFUN = wgtFUN)

  nout.ppmcd   = sum(ppmcd$best > n.in)
  nout.fastmcd = sum(fastmcd$best > n.in)
  nout.detmcd  = sum(detmcd$best > n.in)

  purity.ppmcd[n - n.in]   = 1 - nout.ppmcd/length(ppmcd$best)
  purity.fastmcd[n - n.in] = 1 - nout.fastmcd/length(fastmcd$best)
  purity.detmcd[n - n.in]  = 1 - nout.detmcd/length(detmcd$best)

  limit.max = if (asymptotic.quantiles.flag) qchisq(1 - (1 - 0.05)^n, df = p) else ppmcd.quantiles.max.rew.0.95[n - n.in + 1]

  if (max(ppmcd$mah) > limit.max) {
    limit.rnd = if (asymptotic.quantiles.flag) qchisq(1 - 0.05, df = p) else ppmcd.quantiles.rnd.rew.0.95[n - n.in + 1]
    
    TP.ppmcd[n - n.in] = sum(ppmcd$mah[(n.in + 1):n] > limit.rnd)/(n - n.in)
    FP.ppmcd[n - n.in] = sum(ppmcd$mah[1:n.in] > limit.rnd)/n.in
  }
  
  limit.max = if (asymptotic.quantiles.flag) qchisq(1 - (1 - 0.05)^n, df = p) else fastmcd.quantiles.max.rew.0.95[n - n.in + 1]

  if (max(fastmcd$mah) > limit.max) {
    limit.rnd = if (asymptotic.quantiles.flag) qchisq(1 - 0.05, df = p) else fastmcd.quantiles.rnd.rew.0.95[n - n.in + 1]
    
    TP.fastmcd[n - n.in] = sum(fastmcd$mah[(n.in + 1):n] > limit.rnd)/(n - n.in)
    FP.fastmcd[n - n.in] = sum(fastmcd$mah[1:n.in] > limit.rnd)/n.in
  }
  
  limit.max = if (asymptotic.quantiles.flag) qchisq(1 - (1 - 0.05)^n, df = p) else detmcd.quantiles.max.rew.0.95[n - n.in + 1]

  if (max(detmcd$mah) > limit.max) {
    limit.rnd = if (asymptotic.quantiles.flag) qchisq(1 - 0.05, df = p) else detmcd.quantiles.rnd.rew.0.95[n - n.in + 1]
    
    TP.detmcd[n - n.in] = sum(detmcd$mah[(n.in + 1):n] > limit.rnd)/(n - n.in)
    FP.detmcd[n - n.in] = sum(detmcd$mah[1:n.in] > limit.rnd)/n.in
  }
}

## TP and FP

  file.name = paste("TP.FP.graphs.", key, ".eps", sep = "")
  grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 10, height = 4, horizontal = TRUE)

  par(mfrow = c(1, 2))

  plot(TP.ppmcd*100, type = "l", lty = 1, lwd = 1, col = "black", xlab = "Number of outliers in the sample", ylab = "True positive rate (TPR) (in %)",
       xlim = c(0, n.in), ylim = c(0, 100)) #, xaxt = "n")
  xtick = c(0, 10, 20, 30, 40, 50, 77, 88, 100)
  axis(side = 1, at = xtick, labels = TRUE)
  
  lines(TP.fastmcd*100, col = "blue", lty = 2, lwd = 1)
  lines(TP.detmcd*100, col = "red", lty = 4, lwd = 1)
  
  points(88, TP.ppmcd[88]*100, pch = 1, col = "black")
  points(77, TP.fastmcd[77]*100, pch = 0, col = "blue")
  points(50, TP.detmcd[50]*100, pch = 2, col = "red")
  legend("bottomleft", legend = c("PP MCD", "FastMCD", "DetMCD"), col = c("black", "blue", "red"), lty = c(1, 2, 4), pch = c(1, 0, 2))

  plot(FP.ppmcd*100, type = "l", lty = 1, lwd = 1, col = "black", xlab = "Number of outliers in the sample", ylab = "False positive rate (FPR) (in %)",
       xlim = c(0, n.in), ylim = c(0, 100))
  lines(FP.fastmcd*100, col = "blue", lty = 2, lwd = 1)
  lines(FP.detmcd*100, col = "red", lty = 4, lwd = 1)
  legend("topleft", legend = c("PP MCD", "FastMCD", "DetMCD"), col = c("black", "blue", "red"), lty = c(1, 2, 4), lwd = c(2, 2, 2), pch = c(1, 0, 2))

  grDevices::dev.off()

# Purity

  file.name = paste("bulk.purity.graphs.", key, ".eps", sep = "")
  grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 5, height = 4, horizontal = FALSE)

  par(mfrow = c(1, 1))

  plot(purity.ppmcd*100, type = "l", lty = 1, lwd = 1, col = "black", xlab = "Number of outliers in the sample", ylab = "Bulk purity (in % of bulk size)",
       xlim = c(0, n.in), ylim = c(0, 100)) #, xaxt = "n")
  xtick = c(0, 10, 20, 30, 40, 50, 77, 88, 100)
  axis(side = 1, at = xtick, labels = TRUE)

  lines(purity.fastmcd*100, col = "blue", lty = 2, lwd = 1)
  lines(purity.detmcd*100, col = "red", lty = 4, lwd = 1)

  points(50, purity.detmcd[50]*100.0, col = "red", pch = 2)
  points(77, purity.fastmcd[77]*100.0, col = "blue", pch = 0)
  points(88, purity.ppmcd[88]*100, col = "black", pch = 1)

  legend("bottomleft", legend = c("PP MCD", "FastMCD", "DetMCD"), col = c("black", "blue", "red"), lty = c(1, 2, 4), lwd = c(2, 2, 2), pch = c(1, 0, 2))

  grDevices::dev.off()

# Purity, TP, FP

    file.name = paste("bulk.purity.TP.FP.graphs.", key, ".eps", sep = "")
    grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 12, height = 4, horizontal = TRUE)

    par(mfrow = c(1, 3))

    plot(purity.ppmcd*100, type = "l", lty = 1, lwd = 1, col = "black", 
         main = "(a) Bulk purity", xlab = "Number of outliers in the sample", ylab = "Bulk purity (in % of bulk size)",
         xlim = c(0, n.in), ylim = c(0, 100), xaxt = "n")
    xtick = c(0, 10, 20, 30, 40, 50, 77, 88, 100)
    axis(side = 1, at = xtick, labels = TRUE)

    lines(purity.fastmcd*100, col = "blue", lty = 2, lwd = 1)
    lines(purity.detmcd*100, col = "red", lty = 4, lwd = 1)

    points(50, purity.detmcd[50]*100.0, col = "red", pch = 2)
    points(77, purity.fastmcd[77]*100.0, col = "blue", pch = 0)
    points(88, purity.ppmcd[88]*100, col = "black", pch = 1)

    legend("bottomleft", legend = c("PP MCD", "FastMCD", "DetMCD"), col = c("black", "blue", "red"), lty = c(1, 2, 4), lwd = c(2, 2, 2), pch = c(1, 0, 2))

    plot(TP.ppmcd*100, type = "l", lty = 1, lwd = 1, col = "black", 
         main = "(b) True positive rate (TPR)", xlab = "Number of outliers in the sample", ylab = "True positive rate (TPR) (in %)",
         xlim = c(0, n.in), ylim = c(0, 100), xaxt = "n")
    xtick = c(0, 10, 20, 30, 40, 50, 77, 88, 100)
    axis(side = 1, at = xtick, labels = TRUE)
    
    lines(TP.fastmcd*100, col = "blue", lty = 2, lwd = 1)
    lines(TP.detmcd*100, col = "red", lty = 4, lwd = 1)
    
    points(88, TP.ppmcd[88]*100, pch = 1, col = "black")
    points(77, TP.fastmcd[77]*100, pch = 0, col = "blue")
    points(50, TP.detmcd[50]*100, pch = 2, col = "red")
    legend("bottomleft", legend = c("PP MCD", "FastMCD", "DetMCD"), col = c("black", "blue", "red"), lty = c(1, 2, 4), lwd = c(2, 2, 2), pch = c(1, 0, 2))

    plot(FP.ppmcd*100, type = "l", lty = 1, lwd = 1, col = "black", 
         main = "(c) False positive rate (FPR)", xlab = "Number of outliers in the sample", ylab = "False positive rate (FPR) (in %)",
         xlim = c(0, n.in), ylim = c(0, 100))
    lines(FP.fastmcd*100, col = "blue", lty = 2, lwd = 1)
    lines(FP.detmcd*100, col = "red", lty = 4, lwd = 1)
    legend("topleft", legend = c("PP MCD", "FastMCD", "DetMCD"), col = c("black", "blue", "red"), lty = c(1, 2, 4), lwd = c(2, 2, 2), pch = c(1, 0, 2))

    grDevices::dev.off()

# Tables

table = cbind(1:n.in, zapsmall(TP.detmcd, 4), zapsmall(TP.fastmcd, 4), zapsmall(TP.ppmcd, 4),
                    zapsmall(FP.detmcd, 4), zapsmall(FP.fastmcd, 4), zapsmall(FP.ppmcd, 4))
table = as.data.frame(table)
colnames(table) = c("Bad pts. cnt.", "TP DetMCD", "TP FastMCD", "TP PP MCD", "FP DetMCD", "FP FastMCD", "FP PP MCD")

write.csv(table, file = paste("TP.FP.table.", key, ".csv", sep = ""))

table = cbind(1:n.in, zapsmall(purity.detmcd, 4), zapsmall(purity.fastmcd, 4), zapsmall(purity.ppmcd, 4))
table = as.data.frame(table)
colnames(table) = c("Bad pts. cnt.", "Bulk purity DetMCD", "Bulk purity FastMCD", "Bulk purity  PP MCD")

write.csv(table, file = paste("bulk.purity.table.", key, ".csv", sep = ""))