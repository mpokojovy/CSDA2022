setwd("???")
source("sim_aux.R")

nsim = 5L

# Just tell the program to run parallel
parallel.flag = TRUE
library(doParallel)
HPC <- makeCluster(detectCores())
registerDoParallel(HPC)

################################

## Perform simulation
mcd.nsamp = 500L

rho.mult.array = c(0.75, sqrt(0.9))

eps.array   = c(0.1, 0.4)
n.array     = c(1, 3, 5, 10, 15, 20, 30, 40)*1000
p.array     = c(10, 60)
delta.array = c(1.0)
r.array     = c(100.0, 200.0)

n.row = length(rho.mult.array)*length(eps.array)*length(n.array)*length(p.array)*length(delta.array)*length(r.array)

nam = c("nsim", "n", "p", "eps", "rho.mult", "rho.mult.sq", "r", "delta", "psi",
        "obj.avg", "obj.std", "err.mu.avg", "err.mu.std",
        "err.Sigma1.avg", "err.Sigma1.std", "err.Sigma2.avg", "err.Sigma2.std", #"err.Sigma3.avg", "err.Sigma3.std", "err.Sigma4.avg", "err.Sigma4.std",
        "breakdown.chance", "impurity.avg", "impurity.std", "n.outlier.avg", "n.outlier.std", "bulk.size.avg", "bulk.size.std",
        "n.C.steps.avg", "n.C.steps.std", "time.avg", "time.std")

n.col = length(nam)

fastmcd.db    = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(fastmcd.db)    = nam
detmcd.db     = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(detmcd.db)     = nam
detmcd.std.db = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(detmcd.std.db) = nam
ppmcd.db      = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(ppmcd.db)      = nam
ppmcd.std.db  = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(ppmcd.std.db)  = nam

row.ind = 0L
ptm <- proc.time()
for (rho.mult in rho.mult.array)
for (n in n.array)
for (p in p.array)
for (i.eps in 1:length(eps.array))
{
  for (i.delta in 1:length(delta.array))
  for (i.r in 1:length(r.array))
  {
    cat("rho.mult = ", rho.mult, ", n = ", n, ", p = ", p, ", eps = ", eps.array[i.eps], ", delta = ", delta.array[i.delta], ", r = ", r.array[i.r], "\n", sep = "")

    row.ind = row.ind + 1L

    fastmcd.db[row.ind, ]    = doSimulation(mcd.nsamp,       nsim = nsim, p = p, n = n,
                                            eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
    detmcd.db[row.ind, ]     = doSimulation("deterministic", nsim = nsim, p = p, n = n,
                                            eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
    detmcd.std.db[row.ind, ] = doSimulation("deterministic", nsim = nsim, p = p, n = n,
                                            eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = TRUE)
    ppmcd.db[row.ind, ]      = doSimulation("PP",            nsim = nsim, p = p, n = n,
                                            eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
    ppmcd.std.db[row.ind, ]  = doSimulation("PP",            nsim = nsim, p = p, n = n,
                                            eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = TRUE)
  }
}

print(proc.time() - ptm)

save(data.chunk, fastmcd.db, detmcd.db, detmcd.std.db, ppmcd.db, ppmcd.std.db, file = paste(getwd(), "/large.n.simulation.RData", sep = ""))

if (parallel.flag)
  stopCluster(HPC)
