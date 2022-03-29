setwd("???")

load("output/small.n.simulation.data.RData")

##
rho.mult.array = c(0.75, sqrt(0.9))

eps.array = c(0.0, 0.1, 0.2, 0.3, 0.4)
n.p.array = t(matrix(c(100, 100, 200, 400, 600,
                       2,   5,  10,  40,  60), nrow = 2, byrow = TRUE))
delta.array = c(1.0, 0.001)

nam = c("n", "p", "eps", "rho.mult", "rho.mult.sq", "delta", "r.fastmcd", "r.detmcd", "r.detmcd.std", "r.ppmcd", "r.ppmcd.std")

n.row = length(rho.mult.array)*length(eps.array)*nrow(n.p.array)*length(delta.array)
n.col = length(nam)

##
worst.r.obj.db        = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(worst.r.obj.db)        = nam
worst.r.err.mu.db     = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(worst.r.err.mu.db)     = nam
worst.r.err.Sigma1.db = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(worst.r.err.Sigma1.db) = nam
worst.r.time.db       = data.frame(matrix(0.0, nrow = n.row, ncol = n.col)); names(worst.r.time.db)       = nam

##
row.ind = 0L

for (rho.mult in rho.mult.array)
for (i.n.p in 1:nrow(n.p.array))
for (i.eps in 1:length(eps.array))
for (i.delta in 1:length(delta.array))
{
  row.ind = row.ind + 1L
  
  n     = n.p.array[i.n.p, 1]
  p     = n.p.array[i.n.p, 2]
  eps   = eps.array[i.eps]
  delta = delta.array[i.delta]
  
  row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) & (fastmcd.db[, "n"] == n) & (fastmcd.db[, "p"] == p) &
                      (fastmcd.db[, "eps"] == eps) & (fastmcd.db[, "delta"] == delta))
  
  r.array = fastmcd.db[row.entries, "r"]
  
  ## no affine equivariance
  worst.r.obj.db[row.ind, c("n", "p", "eps", "rho.mult", "rho.mult.sq", "delta")]        = c(n, p, eps, rho.mult, rho.mult^2, delta)
  worst.r.err.mu.db[row.ind, c("n", "p", "eps", "rho.mult", "rho.mult.sq", "delta")]     = c(n, p, eps, rho.mult, rho.mult^2, delta)
  worst.r.err.Sigma1.db[row.ind, c("n", "p", "eps", "rho.mult", "rho.mult.sq", "delta")] = c(n, p, eps, rho.mult, rho.mult^2, delta)
  worst.r.time.db[row.ind, c("n", "p", "eps", "rho.mult", "rho.mult.sq", "delta")]       = c(n, p, eps, rho.mult, rho.mult^2, delta)
  
  worst.r.obj.db[row.ind, c("r.fastmcd", "r.detmcd", "r.detmcd.std", "r.ppmcd", "r.ppmcd.std")] = c(r.array[which.max(fastmcd.db[row.entries, "obj.avg"])],
                                                                                                    r.array[which.max(detmcd.db[row.entries, "obj.avg"])],
                                                                                                    r.array[which.max(detmcd.std.db[row.entries, "obj.avg"])],
                                                                                                    r.array[which.max(ppmcd.db[row.entries, "obj.avg"])],
                                                                                                    r.array[which.max(ppmcd.std.db[row.entries, "obj.avg"])])
  
  worst.r.err.mu.db[row.ind, c("r.fastmcd", "r.detmcd", "r.detmcd.std", "r.ppmcd", "r.ppmcd.std")] = c(r.array[which.max(fastmcd.db[row.entries, "err.mu.avg"])],
                                                                                                       r.array[which.max(detmcd.db[row.entries, "err.mu.avg"])],
                                                                                                       r.array[which.max(detmcd.std.db[row.entries, "err.mu.avg"])],
                                                                                                       r.array[which.max(ppmcd.db[row.entries, "err.mu.avg"])],
                                                                                                       r.array[which.max(ppmcd.std.db[row.entries, "err.mu.avg"])])
  
  worst.r.err.Sigma1.db[row.ind, c("r.fastmcd", "r.detmcd", "r.detmcd.std", "r.ppmcd", "r.ppmcd.std")] = c(r.array[which.max(fastmcd.db[row.entries, "err.Sigma1.avg"])],
                                                                                                           r.array[which.max(detmcd.db[row.entries, "err.Sigma1.avg"])],
                                                                                                           r.array[which.max(detmcd.std.db[row.entries, "err.Sigma1.avg"])],
                                                                                                           r.array[which.max(ppmcd.db[row.entries, "err.Sigma1.avg"])],
                                                                                                           r.array[which.max(ppmcd.std.db[row.entries, "err.Sigma1.avg"])])
  
  worst.r.time.db[row.ind, c("r.fastmcd", "r.detmcd", "r.detmcd.std", "r.ppmcd", "r.ppmcd.std")] = c(r.array[which.max(fastmcd.db[row.entries, "time.avg"])],
                                                                                                     r.array[which.max(detmcd.db[row.entries, "time.avg"])],
                                                                                                     r.array[which.max(detmcd.std.db[row.entries, "time.avg"])],
                                                                                                     r.array[which.max(ppmcd.db[row.entries, "time.avg"])],
                                                                                                     r.array[which.max(ppmcd.std.db[row.entries, "time.avg"])])
}

save(worst.r.obj.db, worst.r.err.mu.db, worst.r.err.Sigma1.db,  worst.r.time.db, file = "worst.r.RData")
