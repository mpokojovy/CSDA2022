library(mclust)
library(PPcovMcd)

setwd("???")

asymptotic.quantiles.flag = FALSE

if (!asymptotic.quantiles.flag) {
  load("quantiles.raw.RData")
  load("quantiles.rew.RData")
}

####

x.train = read.csv("segmentation.data", header = FALSE, skip = 5)
x.test  = read.csv("segmentation.test", header = FALSE, skip = 5)

x = rbind(x.train, x.test)

#I.var = c(2:3, 7:17, 19, 20)
I.var = c(7:17, 19, 20)
#I.var = c(7:20)

for (j in I.var) {
  x[, j] = (x[, j] - mean(x[, j]))/sd(x[, j])
}

####

class.in  = c("SKY")
class.out = c("GRASS")
key = "SKY.vs.GRASS"

# class.in  = c("SKY")
# class.out = setdiff(unique(x[, 1]), c(class.in, "GRASS"))
# key = "SKY.vs.rest"

x = x[is.element(x[, 1], union(class.in, class.out)), ]
y = ifelse(is.element(x[, 1], class.in), 0, 1)

x = x[, I.var]

I.in  = which(y == 0)
I.out = which(y == 1)

set.seed(123)

n.in  = length(I.in)

I.in  = sample(I.in,  n.in)
I.out = sample(I.out, n.in)

dataset = x[c(I.in, I.out), ]

p     = ncol(dataset)
n.max = nrow(dataset)

###

n.array = 330 + c(100) # c(50, 100, 150, 200, 250, 300, 330)

if (key == "SKY.vs.GRASS") {
  n.pairs.array = matrix(c(330 + 253, 330 + 303, 330 + 308, 330 + 321,
                           330 + 254, 330 + 304, 330 + 309, 330 + 322), ncol = 2)
} else {
  n.pairs.array = matrix(c(330 + 129, 330 + 175, 330 + 192, 330 + 305, 330 + 311,
                           330 + 130, 330 + 176, 330 + 193, 330 + 306, 330 + 312), ncol = 2)
}

n.min = 330

for (n in n.array) {
  file.name = paste(n - n.min, ".bad.", key, ".eps", sep = "")
  grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 10, height = 4, horizontal = FALSE)
  
  x = as.matrix(dataset)
  rownames(x) = NULL
  colnames(x) = NULL
  
  p = ncol(x)
  
  set.seed(n - n.min)
  
  x = x[1:n, ]
  
  mah.std = mah.standard(x)
  
  if (asymptotic.quantiles.flag)
    wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
  else
    wgtFUN <- function(dist) (dist < ppmcd.quantiles.rnd.raw.0.9[n - n.in + 1]*median(dist)/ppmcd.quantiles.rnd.raw.0.5[n - n.in + 1])
  ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            wgtFUN = wgtFUN)
  
  if (asymptotic.quantiles.flag)
    wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
  else
    wgtFUN <- function(dist) (dist < fastmcd.quantiles.rnd.raw.0.9[n - n.in + 1]*median(dist)/fastmcd.quantiles.rnd.raw.0.5[n - n.in + 1])
  fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             wgtFUN = wgtFUN)
  
  if (asymptotic.quantiles.flag)
    wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
  else
    wgtFUN <- function(dist) (dist < detmcd.quantiles.rnd.raw.0.9[n - n.in + 1]*median(dist)/detmcd.quantiles.rnd.raw.0.5[n - n.in + 1])
  detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", wgtFUN = wgtFUN)
  
  par(mfrow = c(1, 3))
  
  md.max = if (key == "SKY.vs.GRASS") max(ppmcd$mah, fastmcd$mah, detmcd$mah)*1.1 else quantile(c(ppmcd$mah, fastmcd$mah, detmcd$mah), 0.95)*1.1
  
  # PP MCD
  if (max(ppmcd$mah) > ppmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
    TP = sum(ppmcd$mah[(n.min + 1):n] > ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    FP = sum(ppmcd$mah[1:n.min] > ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
  } else {
    TP = 0.0
    FP = 0.0
  }
  
  str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
  plot(ppmcd$mah/exp(ppmcd.logdet.raw[n - n.min + 1]/p), col = "black", pch = 1, ylab = "Sq. Mahalanobis dist.", main = paste("PP MCD \n", str), ylim = c(0, md.max))
  abline(a = ppmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(ppmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "black", lty = 1, lwd = 2)
  abline(a = ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1]/exp(ppmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "black", lty = 1, lwd = 1)
  abline(v = n.min, col = "black", lty = 3)
  legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
         col = c("black", "black", "black"), lty = c(NA, 1, 1), lwd = c(NA, 2, 1), pch = c(1, NA, NA))
  
  # FastMCD
  if (max(fastmcd$mah) > fastmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
    TP = sum(fastmcd$mah[(n.min + 1):n] > fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    FP = sum(fastmcd$mah[1:n.min] > fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
  } else {
    TP = 0.0
    FP = 0.0
  }
  
  str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
  plot(fastmcd$mah/exp(fastmcd.logdet.raw[n - n.min + 1]/p), col = "blue", pch = 0, ylab = "Sq. Mahalanobis dist.", main = paste("FastMCD \n", str), ylim = c(0, md.max))
  abline(a = fastmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(fastmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "blue", lty = 2, lwd = 2)
  abline(a = fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1]/exp(fastmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "blue", lty = 2, lwd = 1)
  abline(v = n.min, col = "black", lty = 3)
  legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
         col = c("blue", "blue", "blue"), lty = c(NA, 2, 2), lwd = c(NA, 2, 1), pch = c(0, NA, NA))
  
  # DetMCD
  if (max(detmcd$mah) > detmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
    TP = sum(detmcd$mah[(n.min + 1):n] > detmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    FP = sum(detmcd$mah[1:n.min] > detmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
  } else {
    TP = 0.0
    FP = 0.0
  }
  
  str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
  plot(detmcd$mah/exp(detmcd.logdet.raw[n - n.min + 1]/p), col = "red", pch = 2, ylab = "Sq. Mahalanobis dist.", main = paste("DetMCD \n", str), ylim = c(0, md.max))
  abline(a = detmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(detmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "red", lty = 4, lwd = 2)
  abline(a = detmcd.quantiles.rnd.rew.0.95[n - n.min + 1]/exp(detmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "red", lty = 4, lwd = 1)
  abline(v = n.min, col = "black", lty = 3)
  legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
         col = c("red", "red", "red"), lty = c(NA, 4, 4), lwd = c(NA, 2, 1), pch = c(2, NA, NA))
  
  grDevices::dev.off()
}

for (i in 1:nrow(n.pairs.array)) {
  file.name = paste(n.pairs.array[i, 1] - n.min, ".", n.pairs.array[i, 2] - n.min, ".bad.", key, ".eps", sep = "")
  grDevices::postscript(file = paste(getwd(), "/fig/", file.name, sep = ""), width = 10, height = 8, horizontal = TRUE)
  
  par(mfrow = c(2, 3))
  
  for (n in n.pairs.array[i, ]) {
    x = as.matrix(dataset)
    
    p = ncol(x)
    
    set.seed(n - n.min)
    
    x = x[1:n, ]
    
    mah.std = mah.standard(x)
    
    if (asymptotic.quantiles.flag)
      wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
    else
      wgtFUN <- function(dist) (dist < ppmcd.quantiles.rnd.raw.0.9[n - n.in + 1]*median(dist)/ppmcd.quantiles.rnd.raw.0.5[n - n.in + 1])
    ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            wgtFUN = wgtFUN)
    
    if (asymptotic.quantiles.flag)
      wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
    else
      wgtFUN <- function(dist) (dist < fastmcd.quantiles.rnd.raw.0.9[n - n.in + 1]*median(dist)/fastmcd.quantiles.rnd.raw.0.5[n - n.in + 1])
    fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             wgtFUN = wgtFUN)
    
    if (asymptotic.quantiles.flag)
      wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
    else
      wgtFUN <- function(dist) (dist < detmcd.quantiles.rnd.raw.0.9[n - n.in + 1]*median(dist)/detmcd.quantiles.rnd.raw.0.5[n - n.in + 1])
    detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", wgtFUN = wgtFUN)
    
    #md.max = max(ppmcd$mah, fastmcd$mah, detmcd$mah)*1.1
    
    md.max = if (key == "SKY.vs.GRASS") max(ppmcd$mah)/exp(ppmcd.logdet.raw[n - n.min + 1]/p)*1.1 else ppmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(ppmcd.logdet.raw[n - n.min + 1]/p)*1.25
    
    # PP MCD
    if (max(ppmcd$mah) > ppmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
      TP = sum(ppmcd$mah[(n.min + 1):n] > ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
      FP = sum(ppmcd$mah[1:n.min] > ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    } else {
      TP = 0.0
      FP = 0.0
    }
    
    str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
    plot(ppmcd$mah/exp(ppmcd.logdet.raw[n - n.min + 1]/p), col = "black", pch = 1, xlab = "Index", ylab = "Sq. Mahalanobis dist.", main = paste("PP MCD \n", str), ylim = c(0, md.max))
    abline(a = ppmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(ppmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "black", lty = 1, lwd = 2)
    abline(a = ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1]/exp(ppmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "black", lty = 1, lwd = 1)
    abline(v = n.min, col = "black", lty = 3)
    legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
           col = c("black", "black", "black"), lty = c(NA, 1, 1), lwd = c(NA, 2, 1), pch = c(1, NA, NA))
    
    # md.max = max(fastmcd$mah)*1.1
    md.max = if (key == "SKY.vs.GRASS") max(fastmcd$mah)/exp(fastmcd.logdet.raw[n - n.min + 1]/p)*1.1 else fastmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(ppmcd.logdet.raw[n - n.min + 1]/p)*1.25
    
    # FastMCD
    if (max(fastmcd$mah) > fastmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
      TP = sum(fastmcd$mah[(n.min + 1):n] > fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
      FP = sum(fastmcd$mah[1:n.min] > fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    } else {
      TP = 0.0
      FP = 0.0
    }
    
    str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
    plot(fastmcd$mah/exp(fastmcd.logdet.raw[n - n.min + 1]/p), col = "blue", pch = 0, xlab = "Index", ylab = "Sq. Mahalanobis dist.", main = paste("FastMCD \n", str), ylim = c(0, md.max))
    abline(a = fastmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(fastmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "blue", lty = 2, lwd = 2)
    abline(a = fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1]/exp(fastmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "blue", lty = 2, lwd = 1)
    abline(v = n.min, col = "black", lty = 3)
    legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
           col = c("blue", "blue", "blue"), lty = c(NA, 2, 2), lwd = c(NA, 2, 1), pch = c(0, NA, NA))
    
    # DetMCD
    
    # md.max = max(detmcd$mah)*1.1
    md.max = if (key == "SKY.vs.GRASS") max(detmcd$mah)/exp(detmcd.logdet.raw[n - n.min + 1]/p)*1.1 else detmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(ppmcd.logdet.raw[n - n.min + 1]/p)*1.25
    
    if (max(detmcd$mah) > detmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
      TP = sum(detmcd$mah[(n.min + 1):n] > detmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
      FP = sum(detmcd$mah[1:n.min] > detmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    } else {
      TP = 0.0
      FP = 0.0
    }
    
    str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
    plot(detmcd$mah/exp(detmcd.logdet.raw[n - n.min + 1]/p), col = "red", pch = 2, xlab = "Index", ylab = "Sq. Mahalanobis dist.", main = paste("DetMCD \n", str), ylim = c(0, md.max))
    abline(a = detmcd.quantiles.max.rew.0.95[n - n.min + 1]/exp(detmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "red", lty = 4, lwd = 2)
    abline(a = detmcd.quantiles.rnd.rew.0.95[n - n.min + 1]/exp(detmcd.logdet.raw[n - n.min + 1]/p), b = 0, col = "red", lty = 4, lwd = 1)
    abline(v = n.min, col = "black", lty = 3)
    legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
           col = c("red", "red", "red"), lty = c(NA, 4, 4), lwd = c(NA, 2, 1), pch = c(2, NA, NA))
  }
  
  grDevices::dev.off()
}