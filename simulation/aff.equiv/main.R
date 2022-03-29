setwd("???")
source("sim_aux.R")

set.seed(1)

nsim = 20L # number of random samples
nrep = 50L # number of random transformations per sample

# Just tell the program to run parallel
parallel.flag = TRUE
#parallel.flag = FALSE

if (parallel.flag) {
  library(doParallel)
  #HPC <- makeCluster(detectCores())
  HPC <- makeCluster(detectCores() - 4)
  registerDoParallel(HPC)
} else {
  ##!!!
  #library(PPcovMcd)
}

################################

## Perform simulation
mcd.nsamp = 500L

rho.mult.array = c(0.75, sqrt(0.9))

eps.array = c(0.0, 0.1, 0.2, 0.3, 0.4)
n.p.array = t(matrix(c(100, 100, 200, 400, 600,
                         2,   5,  10,  40,  60), nrow = 2, byrow = TRUE))
delta.array = c(1.0, 0.001)
psi.array   = c(5.0)

r.array = c(100.0, 200.0)
n.row = length(rho.mult.array)*length(eps.array)*nrow(n.p.array)*(length(r.array)*length(delta.array) + length(psi.array))

nam = c("nsim", "n", "p", "eps", "rho.mult", "rho.mult.sq", "r", "delta", "psi",
        "d_mu.avg", "d_mu.std", "d_Sigma.avg", "d_Sigma.std", "nsing")

n.col = length(nam)

fastmcd.db    = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(fastmcd.db)    = nam
detmcd.db     = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(detmcd.db)     = nam
detmcd.std.db = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(detmcd.std.db) = nam
ppmcd.db      = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(ppmcd.db)      = nam
ppmcd.std.db  = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(ppmcd.std.db)  = nam

row.ind = 0L
ptm <- proc.time()
for (rho.mult in rho.mult.array)
for (i.n.p in 1:nrow(n.p.array))
for (i.eps in 1:length(eps.array))
{
  print(paste("rho.mult = ", rho.mult, ", scenario = ", LETTERS[i.n.p], ", eps = ", eps.array[i.eps], sep = ""))

  for (i.delta in 1:length(delta.array)) {
    for (i.r in 1:length(r.array))
    {
      cat("delta = ", delta.array[i.delta], ", r = ", r.array[i.r], "\n", sep = "")

      row.ind = row.ind + 1L

      # fastmcd.db[row.ind, ]    = doSimulation(mcd.nsamp,       nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
      #                                         eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
      # detmcd.db[row.ind, ]     = doSimulation("deterministic", nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
      #                                         eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
      # detmcd.std.db[row.ind, ] = doSimulation("deterministic", nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
      #                                         eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = TRUE)
      ppmcd.db[row.ind, ]      = doSimulation("PP",            nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                              eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = FALSE)
      ppmcd.std.db[row.ind, ]  = doSimulation("PP",            nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                              eps = eps.array[i.eps], rho.mult = rho.mult, r = r.array[i.r], delta = delta.array[i.delta], psi = NA, standardize = TRUE)
    }
  }

  for (psi in psi.array)
  {
    row.ind = row.ind + 1L
    
      # fastmcd.db[row.ind, ]    = doSimulation(mcd.nsamp,       nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
      #                                      eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = TRUE)
      # detmcd.db[row.ind, ]     = doSimulation("deterministic", nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
      #                                      eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = FALSE)
      # detmcd.std.db[row.ind, ] = doSimulation("deterministic", nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
      #                                      eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = TRUE)
    ppmcd.db[row.ind, ]      = doSimulation("PP",            nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                         eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = FALSE)
    ppmcd.std.db[row.ind, ]  = doSimulation("PP",            nsim = nsim, nrep = nrep, p = n.p.array[i.n.p, 2], n = n.p.array[i.n.p, 1],
                                         eps = eps.array[i.eps], rho.mult = rho.mult, r = NA, delta = NA, psi = psi, standardize = TRUE)
  }
}

print(proc.time() - ptm)

save(fastmcd.db, ppmcd.db, ppmcd.std.db, file = paste(getwd(), "/aff.equiv.RData", sep = ""))

#save(data.chunk, fastmcd.db, detmcd.db, detmcd.std.db, ppmcd.db, ppmcd.std.db, file = paste(getwd(), "/large.n.simulation.RData", sep = ""))
#save(fastmcd.db, detmcd.db, detmcd.std.db, ppmcd.db, ppmcd.std.db, file = paste(getwd(), "/aff.equiv.RData", sep = ""))

if (parallel.flag)
    stopCluster(HPC)
