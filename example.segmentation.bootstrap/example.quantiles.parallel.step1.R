library(mclust)
library(PPcovMcd)

logdet <- function(A) {lg = determinant(A, logarithm = TRUE); as.numeric(lg$modulus*lg$sign)}

#setwd("E:/PPMCD.1.0/release.techno/example.segmentation.ht")
setwd("/media/wsadmin/My Passport/PPMCD.1.0/release.techno/example.segmentation.ht")

# source("example.quantiles.parallel.step1.R")
# source("example.quantiles.parallel.step2.R")

####

x.train = read.csv("segmentation.data", header = FALSE, skip = 5)
x.test  = read.csv("segmentation.test", header = FALSE, skip = 5)

x = rbind(x.train, x.test)

I.var = c(7:17, 19, 20)

for (j in I.var) {
  x[, j] = (x[, j] - mean(x[, j]))/sd(x[, j])
}

class  = c("SKY")

x = x[x[, 1] == class, I.var]
x = as.matrix(x)

dataset = x

p = ncol(dataset)
n = nrow(dataset)

####

parallel.flag = TRUE

if (parallel.flag) {
  library(doParallel)
  #HPC <- makeCluster(detectCores() - 4)
  HPC <- makeCluster(detectCores())
  registerDoParallel(HPC)
}

#nrep = 20000
nrep = 5000

p = 13L
n.in  = 330L
n.max = 2*n.in

ppmcd.mah.max.raw = matrix(0.0, nrow = nrep, ncol = n.in + 1)
ppmcd.mah.rnd.raw = matrix(0.0, nrow = nrep, ncol = n.in + 1)
ppmcd.logdet.raw  = matrix(0.0, nrow = nrep, ncol = n.in + 1)

fastmcd.mah.max.raw = matrix(0.0, nrow = nrep, ncol = n.in + 1)
fastmcd.mah.rnd.raw = matrix(0.0, nrow = nrep, ncol = n.in + 1)
fastmcd.logdet.raw  = matrix(0.0, nrow = nrep, ncol = n.in + 1)

detmcd.mah.max.raw = matrix(0.0, nrow = nrep, ncol = n.in + 1)
detmcd.mah.rnd.raw = matrix(0.0, nrow = nrep, ncol = n.in + 1)
detmcd.logdet.raw  = matrix(0.0, nrow = nrep, ncol = n.in + 1)

ptm <- proc.time()

if (parallel.flag) {
  for (n in n.in:n.max) {
    cat("n =", n, "\n")
    
    set.seed(n - n.in + 1)
    
    sim.ans = foreach(ind = 1:nrep, .inorder = FALSE, .packages = c("robustbase", "PPcovMcd"), .combine = "rbind") %dopar% {
      ans = rep(0.0, 9)
      
      x = dataset[sample.int(n.in, n, replace = TRUE), ]
      
      wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
      
      mah.std = mah.standard(x)
      
      ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            wgtFUN = wgtFUN)
      fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             wgtFUN = wgtFUN)
      detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", wgtFUN = wgtFUN)
      
      ans[1] = max(ppmcd$raw.mah)
      ans[2] = ppmcd$raw.mah[sample.int(n, 1)]
      ans[3] = logdet(ppmcd$raw.cov)
      
      ans[4] = max(fastmcd$raw.mah)
      ans[5] = fastmcd$raw.mah[sample.int(n, 1)]
      ans[6] = logdet(fastmcd$raw.cov)
      
      ans[7] = max(detmcd$raw.mah)
      ans[8] = detmcd$raw.mah[sample.int(n, 1)]
      ans[9] = logdet(detmcd$raw.cov)
      
      ans
    }
    
    ppmcd.mah.max.raw[, n - n.in + 1] = sim.ans[, 1]
    ppmcd.mah.rnd.raw[, n - n.in + 1] = sim.ans[, 2]
    ppmcd.logdet.raw[, n - n.in + 1]  = sim.ans[, 3]
    
    fastmcd.mah.max.raw[, n - n.in + 1] = sim.ans[, 4]
    fastmcd.mah.rnd.raw[, n - n.in + 1] = sim.ans[, 5]
    fastmcd.logdet.raw[, n - n.in + 1]  = sim.ans[, 6]
    
    detmcd.mah.max.raw[, n - n.in + 1] = sim.ans[, 7]
    detmcd.mah.rnd.raw[, n - n.in + 1] = sim.ans[, 8]
    detmcd.logdet.raw[, n - n.in + 1]  = sim.ans[, 9]
  } 
} else {
  for (i in 1:nrep) {
    set.seed(i)
    cat("rep =", i, "\n")
    
    for (n in n.in:n.max) {
      x = dataset[sample.int(n.in, n, replace = TRUE), ]
      
      wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
      
      mah.std = mah.standard(x)
      
      ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            raw.only = TRUE, wgtFUN = wgtFUN)
      fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             raw.only = TRUE, wgtFUN = wgtFUN)
      detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", raw.only = TRUE, wgtFUN = wgtFUN)
      
      ppmcd.mah.max.raw[i, n - n.in + 1] = max(ppmcd$raw.mah)
      ppmcd.mah.rnd.raw[i, n - n.in + 1] = ppmcd$raw.mah[sample.int(n, 1)]
      ppmcd.logdet.raw[i, n - n.in + 1]  = logdet(ppmcd$raw.cov)
      
      fastmcd.mah.max.raw[i, n - n.in + 1] = max(fastmcd$raw.mah)
      fastmcd.mah.rnd.raw[i, n - n.in + 1] = fastmcd$raw.mah[sample.int(n, 1)]
      fastmcd.logdet.raw[i, n - n.in + 1]  = logdet(fastmcd$raw.cov)
      
      detmcd.mah.max.raw[i, n - n.in + 1] = max(detmcd$raw.mah)
      detmcd.mah.rnd.raw[i, n - n.in + 1] = detmcd$raw.mah[sample.int(n, 1)]
      detmcd.logdet.raw[i, n - n.in + 1]  = logdet(detmcd$raw.cov)
    }  
  }  
}

if (parallel.flag)
    stopCluster(HPC)

print(proc.time() - ptm)

ppmcd.quantiles.max.raw.0.5 = rep(0.0, n.in + 1)
ppmcd.quantiles.rnd.raw.0.5 = rep(0.0, n.in + 1)
fastmcd.quantiles.max.raw.0.5 = rep(0.0, n.in + 1)
fastmcd.quantiles.rnd.raw.0.5 = rep(0.0, n.in + 1)
detmcd.quantiles.max.raw.0.5 = rep(0.0, n.in + 1)
detmcd.quantiles.rnd.raw.0.5 = rep(0.0, n.in + 1)

ppmcd.quantiles.max.raw.0.9 = rep(0.0, n.in + 1)
ppmcd.quantiles.rnd.raw.0.9 = rep(0.0, n.in + 1)
fastmcd.quantiles.max.raw.0.9 = rep(0.0, n.in + 1)
fastmcd.quantiles.rnd.raw.0.9 = rep(0.0, n.in + 1)
detmcd.quantiles.max.raw.0.9 = rep(0.0, n.in + 1)
detmcd.quantiles.rnd.raw.0.9 = rep(0.0, n.in + 1)

for (n in n.in:n.max) {
  ppmcd.quantiles.max.raw.0.5[n - n.in + 1] = quantile(ppmcd.mah.max.raw[, n - n.in + 1], 0.5)
  ppmcd.quantiles.rnd.raw.0.5[n - n.in + 1] = quantile(ppmcd.mah.rnd.raw[, n - n.in + 1], 0.5)
  fastmcd.quantiles.max.raw.0.5[n - n.in + 1] = quantile(fastmcd.mah.max.raw[, n - n.in + 1], 0.5)
  fastmcd.quantiles.rnd.raw.0.5[n - n.in + 1] = quantile(fastmcd.mah.rnd.raw[, n - n.in + 1], 0.5)
  detmcd.quantiles.max.raw.0.5[n - n.in + 1] = quantile(detmcd.mah.max.raw[, n - n.in + 1], 0.5)
  detmcd.quantiles.rnd.raw.0.5[n - n.in + 1] = quantile(detmcd.mah.rnd.raw[, n - n.in + 1], 0.5)
  
  ppmcd.quantiles.max.raw.0.9[n - n.in + 1] = quantile(ppmcd.mah.max.raw[, n - n.in + 1], 0.9)
  ppmcd.quantiles.rnd.raw.0.9[n - n.in + 1] = quantile(ppmcd.mah.rnd.raw[, n - n.in + 1], 0.9)
  fastmcd.quantiles.max.raw.0.9[n - n.in + 1] = quantile(fastmcd.mah.max.raw[, n - n.in + 1], 0.9)
  fastmcd.quantiles.rnd.raw.0.9[n - n.in + 1] = quantile(fastmcd.mah.rnd.raw[, n - n.in + 1], 0.9)
  detmcd.quantiles.max.raw.0.9[n - n.in + 1] = quantile(detmcd.mah.max.raw[, n - n.in + 1], 0.9)
  detmcd.quantiles.rnd.raw.0.9[n - n.in + 1] = quantile(detmcd.mah.rnd.raw[, n - n.in + 1], 0.9)
}

smooth.monotonic <- function(x, y) {
  df.tmp = data.frame(cbind(x, y))
  colnames(df.tmp) = c("x", "y")
  rownames(df.tmp) = NULL
  
  return(cobs::cobs(df.tmp$x, df.tmp$y, constraint = "decrease", nknots = 200)$fitted)
}

ppmcd.quantiles.max.raw.0.5 = smooth.monotonic(n.in:n.max, ppmcd.quantiles.max.raw.0.5)
ppmcd.quantiles.rnd.raw.0.5 = smooth.monotonic(n.in:n.max, ppmcd.quantiles.rnd.raw.0.5)
fastmcd.quantiles.max.raw.0.5 = smooth.monotonic(n.in:n.max, fastmcd.quantiles.max.raw.0.5)
fastmcd.quantiles.rnd.raw.0.5 = smooth.monotonic(n.in:n.max, fastmcd.quantiles.rnd.raw.0.5)
detmcd.quantiles.max.raw.0.5 = smooth.monotonic(n.in:n.max, detmcd.quantiles.max.raw.0.5)
detmcd.quantiles.rnd.raw.0.5 = smooth.monotonic(n.in:n.max, detmcd.quantiles.rnd.raw.0.5)

ppmcd.quantiles.max.raw.0.9 = smooth.monotonic(n.in:n.max, ppmcd.quantiles.max.raw.0.9)
ppmcd.quantiles.rnd.raw.0.9 = smooth.monotonic(n.in:n.max, ppmcd.quantiles.rnd.raw.0.9)
fastmcd.quantiles.max.raw.0.9 = smooth.monotonic(n.in:n.max, fastmcd.quantiles.max.raw.0.9)
fastmcd.quantiles.rnd.raw.0.9 = smooth.monotonic(n.in:n.max, fastmcd.quantiles.rnd.raw.0.9)
detmcd.quantiles.max.raw.0.9 = smooth.monotonic(n.in:n.max, detmcd.quantiles.max.raw.0.9)
detmcd.quantiles.rnd.raw.0.9 = smooth.monotonic(n.in:n.max, detmcd.quantiles.rnd.raw.0.9)

ppmcd.logdet.raw   = -smooth.monotonic(n.in:n.max, -(logdet(cov(mah.standard(dataset)$z.scores)) - colMeans(ppmcd.logdet.raw)))
fastmcd.logdet.raw = -smooth.monotonic(n.in:n.max, -(logdet(cov(mah.standard(dataset)$z.scores)) - colMeans(fastmcd.logdet.raw)))
detmcd.logdet.raw  = -smooth.monotonic(n.in:n.max, -(logdet(cov(mah.standard(dataset)$z.scores)) - colMeans(detmcd.logdet.raw)))

save(ppmcd.mah.max.raw, ppmcd.mah.rnd.raw, fastmcd.mah.max.raw, fastmcd.mah.rnd.raw, detmcd.mah.max.raw, detmcd.mah.rnd.raw,
     ppmcd.quantiles.max.raw.0.5, ppmcd.quantiles.rnd.raw.0.5, fastmcd.quantiles.max.raw.0.5, 
     fastmcd.quantiles.rnd.raw.0.5, detmcd.quantiles.max.raw.0.5, detmcd.quantiles.rnd.raw.0.5, 
     ppmcd.quantiles.max.raw.0.9, ppmcd.quantiles.rnd.raw.0.9, fastmcd.quantiles.max.raw.0.9, 
     fastmcd.quantiles.rnd.raw.0.9, detmcd.quantiles.max.raw.0.9, detmcd.quantiles.rnd.raw.0.9,
     ppmcd.logdet.raw, fastmcd.logdet.raw, detmcd.logdet.raw, file = "quantiles.raw.RData")

## Correction factor
file.name = paste("cor.factor.raw.graphs.eps", sep = "")
grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 12, height = 4, horizontal = FALSE)

par(mfrow = c(1, 3))

plot(n.in:n.max, exp(ppmcd.logdet.raw/p), type = "l", lty = 1, col = "black", lwd = 2, xlab = "Sample size n", ylab = "Raw covariance bias correction factor",
     xlim = c(n.in, n.max), ylim = c(0, 3), main = "PP MCD")
#lines(n.in:n.max, ppmcd.quantiles.rnd.raw.0.5, lty = 1, col = "black")
abline(h = 1.0, col = "black", lty = 3)
legend(x = "topright", legend = c("Bias correction factor", "Baseline = 1.0"), col = c("black", "black"), lty = c(1, 1), lwd = c(2, 1))

plot(n.in:n.max, exp(fastmcd.logdet.raw/p), type = "l", lty = 2, col = "blue", lwd = 2, xlab = "Sample size n", ylab = "Raw covariance bias correction factor",
     xlim = c(n.in, n.max), ylim = c(0, 3), main = "FastMCD")
#lines(n.in:n.max, fastmcd.quantiles.rnd.raw.0.5, lty = 2, col = "blue")
abline(h = 1.0, col = "blue", lty = 3)
legend(x = "topright", legend = c("Bias correction factor", "Baseline = 1.0"), col = c("blue", "blue"), lty = c(2, 2), lwd = c(2, 1))

plot(n.in:n.max, exp(detmcd.logdet.raw/p), type = "l", lty = 4, col = "red", lwd = 2, xlab = "Sample size n", ylab = "Raw covariance bias correction factor",
     xlim = c(n.in, n.max), ylim = c(0, 3), main = "DetMCD")
#lines(n.in:n.max, detmcd.quantiles.rnd.raw.0.5, lty = 4, col = "red")
abline(h = 1.0, col = "red", lty = 3)
legend(x = "topright", legend = c("Bias correction factor", "Baseline = 1.0"), col = c("red", "red"), lty = c(4, 4), lwd = c(2, 1))

grDevices::dev.off()

## 50-th quantile
file.name = paste("50.quantile.raw.graphs.eps", sep = "")
grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 12, height = 4, horizontal = FALSE)

par(mfrow = c(1, 3))

plot(n.in:n.max, ppmcd.quantiles.rnd.raw.0.5/exp(ppmcd.logdet.raw/p), 
     type = "l", lty = 1, col = "black", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 50-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 20), main = "PP MCD")
#lines(n.in:n.max, ppmcd.quantiles.rnd.raw.0.5, lty = 1, col = "black")
abline(h = qchisq(0.5, p), col = "black", lty = 3)
legend(x = "topright", legend = c("random MD 0.5", paste("chi-sq. 0.5 w/ df =", p)), col = c("black", "black"), lty = c(1, 3), lwd = c(2, 1))

plot(n.in:n.max, fastmcd.quantiles.rnd.raw.0.5/exp(fastmcd.logdet.raw/p), 
     type = "l", lty = 2, col = "blue", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 50-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 20), main = "FastMCD")
#lines(n.in:n.max, fastmcd.quantiles.rnd.raw.0.5, lty = 2, col = "blue")
abline(h = qchisq(0.5, p), col = "blue", lty = 3)
legend(x = "topright", legend = c("random MD 0.5", paste("chi-sq. 0.5 w/ df =", p)), col = c("blue", "blue"), lty = c(2, 3), lwd = c(2, 1))

plot(n.in:n.max, detmcd.quantiles.rnd.raw.0.5/exp(detmcd.logdet.raw/p), 
     type = "l", lty = 4, col = "red", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 50-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 20), main = "DetMCD")
#lines(n.in:n.max, detmcd.quantiles.rnd.raw.0.5, lty = 4, col = "red")
abline(h = qchisq(0.5, p), col = "red", lty = 3)
legend(x = "topright", legend = c("random MD 0.5", paste("chi-sq. 0.5 w/ df =", p)), col = c("red", "red"), lty = c(4, 3), lwd = c(2, 1))

grDevices::dev.off()

## 90-th quantile
file.name = paste("90.quantile.raw.graphs.eps", sep = "")
grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 12, height = 4, horizontal = FALSE)

par(mfrow = c(1, 3))

plot(n.in:n.max, ppmcd.quantiles.rnd.raw.0.9/exp(ppmcd.logdet.raw/p), 
     type = "l", lty = 1, col = "black", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 90-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 250), main = "PP MCD")
#lines(n.in:n.max, ppmcd.quantiles.rnd.raw.0.5, lty = 1, col = "black")
abline(h = qchisq(0.9, p), col = "black", lty = 3)
legend(x = "topright", legend = c("random MD 0.9", paste("chi-sq. 0.9 w/ df =", p)), col = c("black", "black"), lty = c(1, 3), lwd = c(2, 1))

plot(n.in:n.max, fastmcd.quantiles.rnd.raw.0.9/exp(fastmcd.logdet.raw/p), 
     type = "l", lty = 2, col = "blue", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 90-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 250), main = "FastMCD")
#lines(n.in:n.max, fastmcd.quantiles.rnd.raw.0.5, lty = 2, col = "blue")
abline(h = qchisq(0.9, p), col = "blue", lty = 3)
legend(x = "topright", legend = c("random MD 0.9", paste("chi-sq. 0.9 w/ df =", p)), col = c("blue", "blue"), lty = c(2, 3), lwd = c(2, 1))

plot(n.in:n.max, detmcd.quantiles.rnd.raw.0.9/exp(detmcd.logdet.raw/p), 
     type = "l", lty = 4, col = "red", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 90-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 250), main = "DetMCD")
#lines(n.in:n.max, detmcd.quantiles.rnd.raw.0.5, lty = 4, col = "red")
abline(h = qchisq(0.9, p), col = "red", lty = 3)
legend(x = "topright", legend = c("random MD 0.9", paste("chi-sq. 0.9 w/ df =", p)), col = c("red", "red"), lty = c(4, 3), lwd = c(2, 1))

grDevices::dev.off()