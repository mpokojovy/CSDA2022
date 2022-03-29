#setwd("???")
source("sim_aux.R")

nsim = 5L

# Just tell the program to run parallel
#parallel.flag = FALSE
parallel.flag = TRUE

library(doParallel)
library(PPcovMcd)

if (parallel.flag) {
  #HPC <- makeCluster(detectCores())
  HPC <- makeCluster(68L)
  registerDoParallel(HPC)
}

################################

## Perform simulation
mcd.nsamp = 500L

#rho.mult.array = c(0.75, sqrt(0.9))
rho.mult.array = c(0.75)

eps.array = c(0.0, 0.1, 0.2, 0.3, 0.4)
n.p.array = t(matrix(c(100, 100, 200, 400, 600,
                         2,   5,  10,  40,  60), nrow = 2, byrow = TRUE))
delta.array = c(1.0, 0.001)
psi.array   = c(5.0)

r.range <- function(p, delta = 1.0, alpha = 0.01) {
  r.min = ceiling(1.2*sqrt(qchisq(1 - alpha, df = p)/p) + delta*sqrt(qchisq(1 - alpha, df = p)/p))
  ans = c(seq(r.min, 20, length.out = 7), seq(30, 250, by = 10))
}
n.row = length(rho.mult.array)*length(eps.array)*nrow(n.p.array)*(length(r.range(p = 1))*length(delta.array) + length(psi.array))

WS.avg.names = sapply(1:7, function(i) paste("WS", i, ".avg", sep = ""))

nam = c("nsim", "n", "p", "eps", "rho.mult", "rho.mult.sq", "r", "delta", "psi",
        "obj.avg", "obj.std", "err.mu.avg", "err.mu.std",
        "err.Sigma1.avg", "err.Sigma1.std", "err.Sigma2.avg", "err.Sigma2.std", #"err.Sigma3.avg", "err.Sigma3.std", "err.Sigma4.avg", "err.Sigma4.std",
        "breakdown.chance", "impurity.avg", "impurity.std", "n.outlier.avg", "n.outlier.std", "bulk.size.avg", "bulk.size.std",
        "n.C.steps.avg", "n.C.steps.std", "time.avg", "time.std", WS.avg.names)

n.col = length(nam)

fastmcd.db       = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(fastmcd.db)       = nam
detmcd.db        = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(detmcd.db)        = nam
detmcd.std.db    = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(detmcd.std.db)    = nam
ppmcd.db         = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(ppmcd.db)         = nam
ppmcd.std.db     = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(ppmcd.std.db)     = nam
ppplusmcd.db     = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(ppplusmcd.db)     = nam
ppplusmcd.std.db = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(ppplusmcd.std.db) = nam

row.ind = 0L
ptm <- proc.time()
for (rho.mult in rho.mult.array)
for (i.n.p in 1:nrow(n.p.array))
for (i.eps in 1:length(eps.array))
{
  print(paste("rho.mult = ", rho.mult, ", scenario = ", LETTERS[i.n.p], ", eps = ", eps.array[i.eps], sep = ""))

  for (i.delta in 1:length(delta.array)) {
    r.array = r.range(p = n.p.array[i.n.p, 2], delta = delta.array[i.delta])

    for (i.r in 1:length(r.array))
    {
      cat("delta = ", delta.array[i.delta], ", r = ", r.array[i.r], "\n", sep = "")

      row.ind = row.ind + 1L

      fastmcd.db[row.ind, ]        = doSimulation(mcd.nsamp,              nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                 eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
      detmcd.db[row.ind, ]         = doSimulation("deterministic",        nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                 eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
      detmcd.std.db[row.ind, ]     = doSimulation("deterministic",        nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                 eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = TRUE)
      ppmcd.db[row.ind, ]          = doSimulation("PP",                   nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                 eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
      ppmcd.std.db[row.ind, ]      = doSimulation("PP",                   nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                 eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = TRUE)
      ppplusmcd.db[row.ind, ]      = doSimulation("PP.plus",               nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                  eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
      ppplusmcd.std.db[row.ind, ]  = doSimulation("PP.plus",               nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                  eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = TRUE)
    }
  }

  for (psi in psi.array)
  {
    row.ind = row.ind + 1L
    
    fastmcd.db[row.ind, ]        = doSimulation(mcd.nsamp,              nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = TRUE)
    detmcd.db[row.ind, ]         = doSimulation("deterministic",        nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = FALSE)
    detmcd.std.db[row.ind, ]     = doSimulation("deterministic",        nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = TRUE)
    ppmcd.db[row.ind, ]          = doSimulation("PP",                   nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = FALSE)
    ppmcd.std.db[row.ind, ]      = doSimulation("PP",                   nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = TRUE)
    ppplusmcd.db[row.ind, ]      = doSimulation("PP.plus",              nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = FALSE)
    ppplusmcd.std.db[row.ind, ]  = doSimulation("PP.plus",              nsim = nsim, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                                eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = TRUE)
  }
}

print(proc.time() - ptm)

save(data.chunk, fastmcd.db, detmcd.db, detmcd.std.db, ppmcd.db, ppmcd.std.db, ppplusmcd.db, ppplusmcd.std.db,
     file = paste(getwd(), "/small.n.simulation.data.chunk=", data.chunk, ".RData", sep = ""))

if (parallel.flag)
  stopCluster(HPC)

# write.csv(detmcd.db,        file = "6pack.csv")
# write.csv(detmcd.std.db,    file = "6pack.standardized.csv")
# write.csv(ppplusmcd.db,     file = "6+1pack.csv")
# write.csv(ppplusmcd.std.db, file = "6+1pack.standardized.csv")
