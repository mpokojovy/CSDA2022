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

n.array = 100 + c(50) # c(50, 100, 150, 200, 250, 300, 330)

n.pairs.array = matrix(c(150, 177, 188,
                         151, 178, 189), ncol = 2)
  

####

n.min = 100

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
  
  wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))

  ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            wgtFUN = wgtFUN)
  fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             wgtFUN = wgtFUN)
  detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", wgtFUN = wgtFUN)
  
  par(mfrow = c(1, 3))
  
  md.max = max(ppmcd$mah, fastmcd$mah, detmcd$mah)*1.1
  
  # PP MCD
  if (max(ppmcd$mah) > ppmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
    TP = sum(ppmcd$mah[(n.min + 1):n] > ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    FP = sum(ppmcd$mah[1:n.min] > ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
  } else {
    TP = 0.0
    FP = 0.0
  }
  
  str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
  plot(ppmcd$mah, col = "black", pch = 1, ylab = "Sq. Mahalanobis dist.", main = paste("PP MCD \n", str), ylim = c(0, md.max))
  abline(a = ppmcd.quantiles.max.rew.0.95[n - n.min + 1], b = 0, col = "black", lty = 1, lwd = 2)
  abline(a = ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1], b = 0, col = "black", lty = 1, lwd = 1)
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
  plot(fastmcd$mah, col = "blue", pch = 0, ylab = "Sq. Mahalanobis dist.", main = paste("FastMCD \n", str), ylim = c(0, md.max))
  abline(a = fastmcd.quantiles.max.rew.0.95[n - n.min + 1], b = 0, col = "blue", lty = 2, lwd = 2)
  abline(a = fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1], b = 0, col = "blue", lty = 2, lwd = 1)
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
  plot(detmcd$mah, col = "red", pch = 2, ylab = "Sq. Mahalanobis dist.", main = paste("DetMCD \n", str), ylim = c(0, md.max))
  abline(a = detmcd.quantiles.max.rew.0.95[n - n.min + 1], b = 0, col = "red", lty = 4, lwd = 2)
  abline(a = detmcd.quantiles.rnd.rew.0.95[n - n.min + 1], b = 0, col = "red", lty = 4, lwd = 1)
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
    
    wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))

    ppmcd   = PPcovMcd(mah.std$z.scores, nsamp = "PP",            wgtFUN = wgtFUN)
    fastmcd = PPcovMcd(mah.std$z.scores, nsamp = 500,             wgtFUN = wgtFUN)
    detmcd  = PPcovMcd(mah.std$z.scores, nsamp = "deterministic", wgtFUN = wgtFUN)
    
    md.max = max(ppmcd$mah, fastmcd$mah, detmcd$mah)*1.1
    
    # PP MCD
    if (max(ppmcd$mah) > ppmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
      TP = sum(ppmcd$mah[(n.min + 1):n] > ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
      FP = sum(ppmcd$mah[1:n.min] > ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    } else {
      TP = 0.0
      FP = 0.0
    }
    
    str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
    plot(ppmcd$mah, col = "black", pch = 1, xlab = "Index", ylab = "Sq. Mahalanobis dist.", main = paste("PP MCD \n", str), ylim = c(0, md.max))
    abline(a = ppmcd.quantiles.max.rew.0.95[n - n.min + 1], b = 0, col = "black", lty = 1, lwd = 2)
    abline(a = ppmcd.quantiles.rnd.rew.0.95[n - n.min + 1], b = 0, col = "black", lty = 1, lwd = 1)
    abline(v = n.min, col = "black", lty = 3)
    legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
           col = c("black", "black", "black"), lty = c(NA, 1, 1), lwd = c(NA, 2, 1), pch = c(1, NA, NA))
    
    # md.max = max(fastmcd$mah)*1.1
    md.max = max(ppmcd$mah, fastmcd$mah, detmcd$mah)*1.1
    
    # FastMCD
    if (max(fastmcd$mah) > fastmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
      TP = sum(fastmcd$mah[(n.min + 1):n] > fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
      FP = sum(fastmcd$mah[1:n.min] > fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    } else {
      TP = 0.0
      FP = 0.0
    }
    
    str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
    plot(fastmcd$mah, col = "blue", pch = 0, xlab = "Index", ylab = "Sq. Mahalanobis dist.", main = paste("FastMCD \n", str), ylim = c(0, md.max))
    abline(a = fastmcd.quantiles.max.rew.0.95[n - n.min + 1], b = 0, col = "blue", lty = 2, lwd = 2)
    abline(a = fastmcd.quantiles.rnd.rew.0.95[n - n.min + 1], b = 0, col = "blue", lty = 2, lwd = 1)
    abline(v = n.min, col = "black", lty = 3)
    legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
           col = c("blue", "blue", "blue"), lty = c(NA, 2, 2), lwd = c(NA, 2, 1), pch = c(0, NA, NA))
    
    # DetMCD
    
    # md.max = max(detmcd$mah)*1.1
    md.max = max(ppmcd$mah, fastmcd$mah, detmcd$mah)*1.1
    
    if (max(detmcd$mah) > detmcd.quantiles.max.rew.0.95[n - n.min + 1]) {
      TP = sum(detmcd$mah[(n.min + 1):n] > detmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
      FP = sum(detmcd$mah[1:n.min] > detmcd.quantiles.rnd.rew.0.95[n - n.min + 1])
    } else {
      TP = 0.0
      FP = 0.0
    }
    
    str = paste("TPR = ", TP, "/", n - n.min, " = ", zapsmall(TP/(n - n.min), digits = 3), ", FPR = ", FP, "/", n.min, " = ", zapsmall(FP/n.min, digits = 3), sep = "")
    plot(detmcd$mah, col = "red", pch = 2, xlab = "Index", ylab = "Sq. Mahalanobis dist.", main = paste("DetMCD \n", str), ylim = c(0, md.max))
    abline(a = detmcd.quantiles.max.rew.0.95[n - n.min + 1], b = 0, col = "red", lty = 4, lwd = 2)
    abline(a = detmcd.quantiles.rnd.rew.0.95[n - n.min + 1], b = 0, col = "red", lty = 4, lwd = 1)
    abline(v = n.min, col = "black", lty = 3)
    legend(x = "topleft", legend = c("sq. MD", "max MD 0.95", "rnd MD 0.95"), 
           col = c("red", "red", "red"), lty = c(NA, 4, 4), lwd = c(NA, 2, 1), pch = c(2, NA, NA))
  }
  
  grDevices::dev.off()
}