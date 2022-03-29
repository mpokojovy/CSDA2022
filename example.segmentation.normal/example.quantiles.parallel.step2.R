library(mclust)
library(PPcovMcd)
  
logdet <- function(A) {lg = determinant(A, logarithm = TRUE); as.numeric(lg$modulus*lg$sign)}

setwd("???")

####

parallel.flag = TRUE

if (parallel.flag) {
  library(doParallel)
  HPC <- makeCluster(detectCores())
  registerDoParallel(HPC)
}

#nrep = 20000
nrep = 5000

p = 13L
n.in  = 330L
n.max = 2*n.in

ppmcd.mah.max.rew = matrix(0.0, nrow = nrep, ncol = n.in + 1)
ppmcd.mah.rnd.rew = matrix(0.0, nrow = nrep, ncol = n.in + 1)
ppmcd.logdet.rew  = matrix(0.0, nrow = nrep, ncol = n.in + 1)

fastmcd.mah.max.rew = matrix(0.0, nrow = nrep, ncol = n.in + 1)
fastmcd.mah.rnd.rew = matrix(0.0, nrow = nrep, ncol = n.in + 1)
fastmcd.logdet.rew  = matrix(0.0, nrow = nrep, ncol = n.in + 1)

detmcd.mah.max.rew = matrix(0.0, nrow = nrep, ncol = n.in + 1)
detmcd.mah.rnd.rew = matrix(0.0, nrow = nrep, ncol = n.in + 1)
detmcd.logdet.rew  = matrix(0.0, nrow = nrep, ncol = n.in + 1)

ptm <- proc.time()

if (parallel.flag) {
  for (n in n.in:n.max) {
    cat("n =", n, "\n")
    
    set.seed(n - n.in + 1)
    
    wgtFUN.ppmcd   <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
    wgtFUN.fastmcd <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
    wgtFUN.detmcd  <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
    
    sim.ans = foreach(ind = 1:nrep, .inorder = FALSE, .packages = c("robustbase", "PPcovMcd"), .combine = "rbind") %dopar% {
      ans = rep(0.0, 9)
      
      x = matrix(rnorm(n*p), ncol = p)
      
      mah.std = mah.standard(x)

      ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            wgtFUN = wgtFUN.ppmcd)
      fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             wgtFUN = wgtFUN.fastmcd)
      detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", wgtFUN = wgtFUN.detmcd)
      
      ans[1] = max(ppmcd$mah)
      ans[2] = ppmcd$mah[sample.int(n, 1)]
      ans[3] = logdet(ppmcd$cov)
      
      ans[4] = max(fastmcd$mah)
      ans[5] = fastmcd$mah[sample.int(n, 1)]
      ans[6] = logdet(fastmcd$cov)
      
      ans[7] = max(detmcd$mah)
      ans[8] = detmcd$mah[sample.int(n, 1)]
      ans[9] = logdet(detmcd$cov)
      
      ans
    }
    
    ppmcd.mah.max.rew[, n - n.in + 1] = sim.ans[, 1]
    ppmcd.mah.rnd.rew[, n - n.in + 1] = sim.ans[, 2]
    ppmcd.logdet.rew[, n - n.in + 1]  = sim.ans[, 3]
    
    fastmcd.mah.max.rew[, n - n.in + 1] = sim.ans[, 4]
    fastmcd.mah.rnd.rew[, n - n.in + 1] = sim.ans[, 5]
    fastmcd.logdet.rew[, n - n.in + 1]  = sim.ans[, 6]
    
    detmcd.mah.max.rew[, n - n.in + 1] = sim.ans[, 7]
    detmcd.mah.rnd.rew[, n - n.in + 1] = sim.ans[, 8]
    detmcd.logdet.rew[, n - n.in + 1]  = sim.ans[, 9]
  } 
} else {
  for (i in 1:nrep) {
    set.seed(i)
    cat("rep =", i, "\n")
    
    for (n in n.in:n.max) {
      x = matrix(rnorm(n*p), ncol = p)
      
      wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
      
      mah.std = mah.standard(x)
      
      ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            wgtFUN = wgtFUN)
      fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             wgtFUN = wgtFUN)
      detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", wgtFUN = wgtFUN)
      
      ppmcd.mah.max.rew[i, n - n.in + 1] = max(ppmcd$mah)
      ppmcd.mah.rnd.rew[i, n - n.in + 1] = ppmcd$mah[sample.int(n, 1)]
      ppmcd.logdet.rew[i, n - n.in + 1]  = logdet(ppmcd$cov)
      
      fastmcd.mah.max.rew[i, n - n.in + 1] = max(fastmcd$mah)
      fastmcd.mah.rnd.rew[i, n - n.in + 1] = fastmcd$mah[sample.int(n, 1)]
      fastmcd.logdet.rew[i, n - n.in + 1]  = logdet(fastmcd$cov)
      
      detmcd.mah.max.rew[i, n - n.in + 1] = max(detmcd$mah)
      detmcd.mah.rnd.rew[i, n - n.in + 1] = detmcd$mah[sample.int(n, 1)]
      detmcd.logdet.rew[i, n - n.in + 1]  = logdet(detmcd$cov)
    }  
  }  
}

if (parallel.flag)
  stopCluster(HPC)

print(proc.time() - ptm)

ppmcd.quantiles.max.rew.0.95 = rep(0.0, n.in + 1)
ppmcd.quantiles.rnd.rew.0.95 = rep(0.0, n.in + 1)
fastmcd.quantiles.max.rew.0.95 = rep(0.0, n.in + 1)
fastmcd.quantiles.rnd.rew.0.95 = rep(0.0, n.in + 1)
detmcd.quantiles.max.rew.0.95 = rep(0.0, n.in + 1)
detmcd.quantiles.rnd.rew.0.95 = rep(0.0, n.in + 1)

for (n in n.in:n.max) {
  ppmcd.quantiles.max.rew.0.95[n - n.in + 1] = quantile(ppmcd.mah.max.rew[, n - n.in + 1], 0.95)
  ppmcd.quantiles.rnd.rew.0.95[n - n.in + 1] = quantile(ppmcd.mah.rnd.rew[, n - n.in + 1], 0.95)
  fastmcd.quantiles.max.rew.0.95[n - n.in + 1] = quantile(fastmcd.mah.max.rew[, n - n.in + 1], 0.95)
  fastmcd.quantiles.rnd.rew.0.95[n - n.in + 1] = quantile(fastmcd.mah.rnd.rew[, n - n.in + 1], 0.95)
  detmcd.quantiles.max.rew.0.95[n - n.in + 1] = quantile(detmcd.mah.max.rew[, n - n.in + 1], 0.95)
  detmcd.quantiles.rnd.rew.0.95[n - n.in + 1] = quantile(detmcd.mah.rnd.rew[, n - n.in + 1], 0.95)
}

smooth.monotonic <- function(x, y) {
  df.tmp = data.frame(cbind(x, y))
  colnames(df.tmp) = c("x", "y")
  rownames(df.tmp) = NULL
  
  return(cobs::cobs(df.tmp$x, df.tmp$y, constraint = "decrease", nknots = 200)$fitted)
}

ppmcd.quantiles.max.rew.0.95 = smooth.monotonic(n.in:n.max, ppmcd.quantiles.max.rew.0.95)
ppmcd.quantiles.rnd.rew.0.95 = smooth.monotonic(n.in:n.max, ppmcd.quantiles.rnd.rew.0.95)
fastmcd.quantiles.max.rew.0.95 = smooth.monotonic(n.in:n.max, fastmcd.quantiles.max.rew.0.95)
fastmcd.quantiles.rnd.rew.0.95 = smooth.monotonic(n.in:n.max, fastmcd.quantiles.rnd.rew.0.95)
detmcd.quantiles.max.rew.0.95 = smooth.monotonic(n.in:n.max, detmcd.quantiles.max.rew.0.95)
detmcd.quantiles.rnd.rew.0.95 = smooth.monotonic(n.in:n.max, detmcd.quantiles.rnd.rew.0.95)

ppmcd.logdet.rew   = -smooth.monotonic(n.in:n.max, colMeans(ppmcd.logdet.rew))
fastmcd.logdet.rew = -smooth.monotonic(n.in:n.max, colMeans(fastmcd.logdet.rew))
detmcd.logdet.rew  = -smooth.monotonic(n.in:n.max, colMeans(detmcd.logdet.rew))

save(ppmcd.mah.max.rew, ppmcd.mah.rnd.rew, fastmcd.mah.max.rew, fastmcd.mah.rnd.rew, detmcd.mah.max.rew, detmcd.mah.rnd.rew,
     ppmcd.quantiles.max.rew.0.95, ppmcd.quantiles.rnd.rew.0.95, fastmcd.quantiles.max.rew.0.95, 
     fastmcd.quantiles.rnd.rew.0.95, detmcd.quantiles.max.rew.0.95, detmcd.quantiles.rnd.rew.0.95, 
     ppmcd.logdet.rew, fastmcd.logdet.rew, detmcd.logdet.rew, file = "quantiles.rew.RData")

## Correction factor
file.name = paste("cor.factor.rew.graphs.eps", sep = "")
grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 12, height = 4, horizontal = FALSE)

par(mfrow = c(1, 3))

plot(n.in:n.max, exp(ppmcd.logdet.rew/p), type = "l", lty = 1, col = "black", lwd = 2, xlab = "Sample size n", ylab = "Reweighted covariance bias correction factor",
     xlim = c(n.in, n.max), ylim = c(0, 2.0), main = "PP MCD")
#lines(n.in:n.max, ppmcd.quantiles.rnd.raw.0.5, lty = 1, col = "black")
abline(h = 1.0, col = "black", lty = 3)
legend(x = "topright", legend = c("Bias correction factor", "Baseline = 1.0"), col = c("black", "black"), lty = c(1, 1), lwd = c(2, 1))

plot(n.in:n.max, exp(fastmcd.logdet.rew/p), type = "l", lty = 2, col = "blue", lwd = 2, xlab = "Sample size n", ylab = "Reweighted covariance bias correction factor",
     xlim = c(n.in, n.max), ylim = c(0, 2.0), main = "FastMCD")
#lines(n.in:n.max, fastmcd.quantiles.rnd.raw.0.5, lty = 2, col = "blue")
abline(h = 1.0, col = "blue", lty = 3)
legend(x = "topright", legend = c("Bias correction factor", "Baseline = 1.0"), col = c("blue", "blue"), lty = c(2, 2), lwd = c(2, 1))

plot(n.in:n.max, exp(detmcd.logdet.rew/p), type = "l", lty = 4, col = "red", lwd = 2, xlab = "Sample size n", ylab = "Reweighted covariance bias correction factor",
     xlim = c(n.in, n.max), ylim = c(0, 2.0), main = "DetMCD")
#lines(n.in:n.max, detmcd.quantiles.rnd.raw.0.5, lty = 4, col = "red")
abline(h = 1.0, col = "red", lty = 3)
legend(x = "topright", legend = c("Bias correction factor", "Baseline = 1.0"), col = c("red", "red"), lty = c(4, 4), lwd = c(2, 1))

grDevices::dev.off()

## 95-th quantile
file.name = paste("95.quantile.rew.graphs.eps", sep = "")
grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 12, height = 4, horizontal = FALSE)

par(mfrow = c(1, 3))

plot(n.in:n.max, ppmcd.quantiles.max.rew.0.95/exp(ppmcd.logdet.rew/p), type = "l", lty = 1, col = "black", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 95-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 40), main = "PP MCD")
lines(n.in:n.max, ppmcd.quantiles.rnd.rew.0.95, lty = 1, col = "black")
abline(h = qchisq(0.95, p), col = "black", lty = 3)
legend(x = "topright", legend = c("max MD 0.95", "random MD 0.95", paste("chi-sq. 0.95 w/ df =", p)), col = c("black", "black", "black"), lty = c(1, 1, 3), lwd = c(2, 1, 1))
  
plot(n.in:n.max, fastmcd.quantiles.max.rew.0.95/exp(fastmcd.logdet.rew/p), type = "l", lty = 2, col = "blue", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 95-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 40), main = "FastMCD")
lines(n.in:n.max, fastmcd.quantiles.rnd.rew.0.95, lty = 2, col = "blue")
abline(h = qchisq(0.95, p), col = "blue", lty = 3)
legend(x = "topright", legend = c("max MD 0.95", "random MD 0.95", paste("chi-sq. 0.95 w/ df =", p)), col = c("blue", "blue", "blue"), lty = c(2, 2, 3), lwd = c(2, 1, 1))

plot(n.in:n.max, detmcd.quantiles.max.rew.0.95/exp(detmcd.logdet.rew/p), type = "l", lty = 4, col = "red", lwd = 2, xlab = "Sample size n", ylab = "Sq. Mahalanobis distance 95-th quantile",
     xlim = c(n.in, n.max), ylim = c(0, 40), main = "DetMCD")
lines(n.in:n.max, detmcd.quantiles.rnd.rew.0.95, lty = 4, col = "red")
abline(h = qchisq(0.95, p), col = "red", lty = 3)
legend(x = "topright", legend = c("max MD 0.95", "random MD 0.95", paste("chi-sq. 0.95 w/ df =", p)), col = c("red", "red", "red"), lty = c(4, 4, 3), lwd = c(2, 1, 1))

grDevices::dev.off()