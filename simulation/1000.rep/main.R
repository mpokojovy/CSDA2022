#setwd("???")

source("sim_aux.R")

# make sure worst.r.RData is in the directory
load("worst.r.RData")

nsim = 5L

# Just tell the program to run parallel
parallel.flag = FALSE
library(doParallel)
HPC <- makeCluster(detectCores())
registerDoParallel(HPC)

################################

## Perform simulation
mcd.nsamp = 500L

rho.mult.array = c(0.75, sqrt(0.9))

eps.array = c(0.0, 0.1, 0.2, 0.3, 0.4)
n.p.array = t(matrix(c(100, 100, 200, 400, 600,
                         2,   5,  10,  40,  60), nrow = 2, byrow = TRUE))
delta.array = c(1.0, 0.001)
psi.array   = c(5.0)

nam = c("nsim", "n", "p", "eps", "rho.mult", "rho.mult.sq", "r", "delta", "psi",
        "obj.avg", "obj.std", "err.mu.avg", "err.mu.std",
        "err.Sigma1.avg", "err.Sigma1.std", "err.Sigma2.avg", "err.Sigma2.std",
        "breakdown.chance", "impurity.avg", "impurity.std", "n.outlier.avg", "n.outlier.std", "bulk.size.avg", "bulk.size.std",
        "n.C.steps.avg", "n.C.steps.std", "time.avg", "time.std")

n.col = length(nam)

fastmcd.db    = NULL
detmcd.db     = NULL
detmcd.std.db = NULL
ppmcd.db      = NULL
ppmcd.std.db  = NULL

ptm <- proc.time()
for (rho.mult in rho.mult.array)
for (i.n.p in 1:nrow(n.p.array))
for (i.eps in 1:length(eps.array))
{
  status.string = paste("rho.mult = ", rho.mult, ", scenario = ", LETTERS[i.n.p], ", eps = ", eps.array[i.eps], sep = "")

  for (i.delta in 1:length(delta.array)) {
    p = n.p.array[i.n.p, 2]
    n = n.p.array[i.n.p, 1]
    eps = eps.array[i.eps] 
    rho.mult = rho.mult
    delta = delta.array[i.delta]
    
    row.ind = which((worst.r.obj.db[, "p"] == p) & (worst.r.obj.db[, "n"] == n) & (worst.r.obj.db[, "eps"] == eps) & 
                    (worst.r.obj.db[, "rho.mult"] == rho.mult) & (worst.r.obj.db[, "delta"] == delta))
    
    if (eps < 1E-6) {
      fastmcd.r.range    = 0.0
      detmcd.r.range     = 0.0
      detmcd.std.r.range = 0.0
      ppmcd.r.range      = 0.0
      ppmcd.std.r.range  = 0.0
      
    } else {
      fastmcd.r.range    = as.vector(unique(c(worst.r.obj.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.err.mu.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.err.Sigma1.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.time.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              250.0)))
      
      detmcd.r.range     = as.vector(unique(c(worst.r.obj.db[row.ind, c("r.detmcd",  "r.ppmcd")],
                                              worst.r.err.mu.db[row.ind, c("r.detmcd",  "r.ppmcd")],
                                              worst.r.err.Sigma1.db[row.ind, c("r.detmcd",  "r.ppmcd")],
                                              worst.r.time.db[row.ind, c("r.detmcd",  "r.ppmcd")],
                                              250.0)))
      
      detmcd.std.r.range = as.vector(unique(c(worst.r.obj.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.err.mu.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.err.Sigma1.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.time.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              250.0)))
      
      ppmcd.r.range      = as.vector(unique(c(worst.r.obj.db[row.ind, c("r.detmcd",  "r.ppmcd")],
                                              worst.r.err.mu.db[row.ind, c("r.detmcd",  "r.ppmcd")],
                                              worst.r.err.Sigma1.db[row.ind, c("r.detmcd",  "r.ppmcd")],
                                              worst.r.time.db[row.ind, c("r.detmcd",  "r.ppmcd")],
                                              250.0)))
      
      ppmcd.std.r.range  = as.vector(unique(c(worst.r.obj.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.err.mu.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.err.Sigma1.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              worst.r.time.db[row.ind, c("r.fastmcd", "r.detmcd.std", "r.ppmcd.std")],
                                              250.0)))
    }
    
    for (r in fastmcd.r.range) {
      cat(status.string, ", delta = ", delta, ": computing FastMCD for r = ", r, "\n", sep = "")
      fastmcd.db    = rbind(fastmcd.db,    unlist(doSimulation(mcd.nsamp, nsim = nsim, p = p, n = n,
                                                               eps = eps, rho.mult = rho.mult, r = r, delta = delta, psi = NA, standardize = FALSE)))
    }
    for (r in detmcd.r.range) {
      cat(status.string, ", delta = ", delta, ": computing DetMCD for r = ", r, "\n", sep = "")
      detmcd.db     = rbind(detmcd.db,     unlist(doSimulation("deterministic", nsim = nsim, p = p, n = n,
                                                               eps = eps, rho.mult = rho.mult, r = r, delta = delta, psi = NA, standardize = FALSE)))
    }
    for (r in detmcd.std.r.range) {
      cat(status.string, ", delta = ", delta, ": computing DetMCD std for r = ", r, "\n", sep = "")
      detmcd.std.db = rbind(detmcd.std.db, unlist(doSimulation("deterministic", nsim = nsim, p = p, n = n,
                                                               eps = eps, rho.mult = rho.mult, r = r, delta = delta, psi = NA, standardize = TRUE)))
    }
    for (r in ppmcd.r.range) {
      cat(status.string, ", delta = ", delta, ": computing PP MCD for r = ", r, "\n", sep = "")
      ppmcd.db      = rbind(ppmcd.db,      unlist(doSimulation("PP", nsim = nsim, p = p, n = n,
                                                               eps = eps, rho.mult = rho.mult, r = r, delta = delta, psi = NA, standardize = FALSE)))
    }
    for (r in ppmcd.std.r.range) {
      cat(status.string, ", delta = ", delta, ": computing PP MCD std for r = ", r, "\n", sep = "")
      ppmcd.std.db  = rbind(ppmcd.std.db,  unlist(doSimulation("PP", nsim = nsim, p = p, n = n,
                                                               eps = eps, rho.mult = rho.mult, r = r, delta = delta, psi = NA, standardize = TRUE)))
    }
    
    if (eps > 1E-6) {
      for (psi in psi.array) {
        cat(status.string, ", delta = ", delta, ": computing FastMCD for psi = ", psi, "\n", sep = "")
        fastmcd.db    = rbind(fastmcd.db,    unlist(doSimulation(mcd.nsamp, nsim = nsim, p = p, n = n,
                                                                 eps = eps, rho.mult = rho.mult, r = NA, delta = delta, psi = psi, standardize = FALSE)))
        
        cat(status.string, ", delta = ", delta, ": computing DetMCD for psi = ", psi, "\n", sep = "")
        detmcd.db     = rbind(detmcd.db,     unlist(doSimulation("deterministic", nsim = nsim, p = p, n = n,
                                                                 eps = eps, rho.mult = rho.mult, r = NA, delta = delta, psi = psi, standardize = FALSE)))
        
        cat(status.string, ", delta = ", delta, ": computing DetMCD std for psi = ", psi, "\n", sep = "")
        detmcd.std.db = rbind(detmcd.std.db, unlist(doSimulation("deterministic", nsim = nsim, p = p, n = n,
                                                                 eps = eps, rho.mult = rho.mult, r = NA, delta = delta, psi = psi, standardize = TRUE)))
        
        cat(status.string, ", delta = ", delta, ": computing PP MCD for psi = ", psi, "\n", sep = "")
        ppmcd.db      = rbind(ppmcd.db,      unlist(doSimulation("PP", nsim = nsim, p = p, n = n,
                                                                 eps = eps, rho.mult = rho.mult, r = NA, delta = delta, psi = psi, standardize = FALSE)))
        
        cat(status.string, ", delta = ", delta, ": computing PP MCD std for psi = ", psi, "\n", sep = "")
        ppmcd.std.db  = rbind(ppmcd.std.db,  unlist(doSimulation("PP", nsim = nsim, p = p, n = n,
                                                                 eps = eps, rho.mult = rho.mult, r = NA, delta = delta, psi = psi, standardize = TRUE)))
      }
    }
  }
}

fastmcd.db    = data.frame(fastmcd.db)
detmcd.db     = data.frame(detmcd.db)
detmcd.std.db = data.frame(detmcd.std.db)
ppmcd.db      = data.frame(ppmcd.db)
ppmcd.std.db  = data.frame(ppmcd.std.db)

names(fastmcd.db)    = nam
names(detmcd.db)     = nam
names(detmcd.std.db) = nam
names(ppmcd.db)      = nam
names(ppmcd.std.db)  = nam

print(proc.time() - ptm)

save(data.chunk, fastmcd.db, detmcd.db, detmcd.std.db, ppmcd.db, ppmcd.std.db, file = paste(getwd(), "/1000.rep.data.chunk=", data.chunk, ".RData", sep = ""))

if (parallel.flag)
  stopCluster(HPC)
