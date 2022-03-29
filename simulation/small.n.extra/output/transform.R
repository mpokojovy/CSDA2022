setwd("???")

no.chunk = 50

load(paste("small.n.simulation.data.chunk=1.RData", sep = ""))

detmcd.DB        = detmcd.db*0
detmcd.std.DB    = detmcd.std.db*0
fastmcd.DB       = fastmcd.db*0
ppmcd.DB         = ppmcd.db*0
ppmcd.std.DB     = ppmcd.std.db*0
ppplusmcd.DB     = ppplusmcd.db*0
ppplusmcd.std.DB = ppplusmcd.std.db*0

## Preset fields
for (varname in c("n", "p", "eps", "rho.mult", "rho.mult.sq", "r", "delta", "psi")) {
  detmcd.DB[, varname]        = detmcd.db[, varname]
  detmcd.std.DB[, varname]    = detmcd.std.db[, varname]
  fastmcd.DB[, varname]       = fastmcd.db[, varname]
  ppmcd.DB[, varname]         = ppmcd.db[, varname]
  ppmcd.std.DB[, varname]     = ppmcd.std.db[, varname]
  ppplusmcd.DB[, varname]     = ppplusmcd.db[, varname]
  ppplusmcd.std.DB[, varname] = ppplusmcd.std.db[, varname]
}

WS.avg.names = sapply(1:7, function(i) paste("WS", i, ".avg", sep = ""))

## Handle average entries
for (i in 1:no.chunk) {
  rm(detmcd.db, detmcd.std.db, fastmcd.db, ppmcd.db, ppmcd.std.db, ppplusmcd.db, ppplusmcd.std.db)
  
  load(paste("small.n.simulation.data.chunk=", i, ".RData", sep = ""))
  
  detmcd.DB$nsim        = detmcd.DB$nsim        + detmcd.db$nsim
  detmcd.std.DB$nsim    = detmcd.std.DB$nsim    + detmcd.std.db$nsim
  fastmcd.DB$nsim       = fastmcd.DB$nsim       + fastmcd.db$nsim
  ppmcd.DB$nsim         = ppmcd.DB$nsim         + ppmcd.db$nsim
  ppmcd.std.DB$nsim     = ppmcd.std.DB$nsim     + ppmcd.std.db$nsim
  ppplusmcd.DB$nsim     = ppplusmcd.DB$nsim     + ppplusmcd.db$nsim
  ppplusmcd.std.DB$nsim = ppplusmcd.std.DB$nsim + ppplusmcd.std.db$nsim
  
  for (varname in c("obj.avg", "err.mu.avg", "err.Sigma1.avg", "err.Sigma2.avg", 
                    "breakdown.chance", "impurity.avg", "n.outlier.avg", "bulk.size.avg", "n.C.steps.avg", "time.avg",
                    WS.avg.names)) {
    detmcd.DB[, varname]        = detmcd.DB[, varname]        + detmcd.db[, varname]*detmcd.db$nsim
    detmcd.std.DB[, varname]    = detmcd.std.DB[, varname]    + detmcd.std.db[, varname]*detmcd.std.db$nsim
    fastmcd.DB[, varname]       = fastmcd.DB[, varname]       + fastmcd.db[, varname]*fastmcd.db$nsim
    ppmcd.DB[, varname]         = ppmcd.DB[, varname]         + ppmcd.db[, varname]*ppmcd.db$nsim
    ppmcd.std.DB[, varname]     = ppmcd.std.DB[, varname]     + ppmcd.std.db[, varname]*ppmcd.std.db$nsim
    ppplusmcd.DB[, varname]     = ppplusmcd.DB[, varname]     + ppplusmcd.db[, varname]*ppplusmcd.db$nsim
    ppplusmcd.std.DB[, varname] = ppplusmcd.std.DB[, varname] + ppplusmcd.std.db[, varname]*ppplusmcd.std.db$nsim
  }
}

for (varname in c("obj.avg", "err.mu.avg", "err.Sigma1.avg", "err.Sigma2.avg", 
                  "breakdown.chance", "impurity.avg", "n.outlier.avg", "bulk.size.avg", "n.C.steps.avg", "time.avg",
                  WS.avg.names)) {
  detmcd.DB[, varname]        = detmcd.DB[, varname]/detmcd.DB$nsim
  detmcd.std.DB[, varname]    = detmcd.std.DB[, varname]/detmcd.std.DB$nsim
  fastmcd.DB[, varname]       = fastmcd.DB[, varname]/fastmcd.DB$nsim
  ppmcd.DB[, varname]         = ppmcd.DB[, varname]/ppmcd.DB$nsim
  ppmcd.std.DB[, varname]     = ppmcd.std.DB[, varname]/ppmcd.std.DB$nsim
  ppplusmcd.DB[, varname]     = ppplusmcd.DB[, varname]/ppplusmcd.DB$nsim
  ppplusmcd.std.DB[, varname] = ppplusmcd.std.DB[, varname]/ppplusmcd.std.DB$nsim
}

## Handle variance entries

I.avg = c("obj.avg", "err.mu.avg", "err.Sigma1.avg", "err.Sigma2.avg", "impurity.avg", "n.outlier.avg", "bulk.size.avg", "n.C.steps.avg", "time.avg")
I.std = c("obj.std", "err.mu.std", "err.Sigma1.std", "err.Sigma2.std", "impurity.std", "n.outlier.std", "bulk.size.std", "n.C.steps.std", "time.std")

for (i in 1:no.chunk) {
  rm(detmcd.db, detmcd.std.db, fastmcd.db, ppmcd.db, ppmcd.std.db, ppplusmcd.db, ppplusmcd.std.db)
  
  load(paste("small.n.simulation.data.chunk=", i, ".RData", sep = ""))
  
  for (j in 1:length(I.avg)) {
    varname.avg = I.avg[j]
    varname.std = I.std[j]
  
    detmcd.DB[, varname.std]        = detmcd.DB[, varname.std]        + detmcd.db[, varname.std]^2*(detmcd.db$nsim - 1)               + detmcd.db[, varname.avg]^2*detmcd.db$nsim
    detmcd.std.DB[, varname.std]    = detmcd.std.DB[, varname.std]    + detmcd.std.db[, varname.std]^2*(detmcd.std.db$nsim - 1)       + detmcd.std.db[, varname.avg]^2*detmcd.std.db$nsim
    fastmcd.DB[, varname.std]       = fastmcd.DB[, varname.std]       + fastmcd.db[, varname.std]^2*(fastmcd.db$nsim - 1)             + fastmcd.db[, varname.avg]^2*fastmcd.db$nsim
    ppmcd.DB[, varname.std]         = ppmcd.DB[, varname.std]         + ppmcd.db[, varname.std]^2*(ppmcd.db$nsim - 1)                 + ppmcd.db[, varname.avg]^2*ppmcd.db$nsim
    ppmcd.std.DB[, varname.std]     = ppmcd.std.DB[, varname.std]     + ppmcd.std.db[, varname.std]^2*(ppmcd.std.db$nsim - 1)         + ppmcd.std.db[, varname.avg]^2*ppmcd.std.db$nsim
    ppplusmcd.DB[, varname.std]     = ppplusmcd.DB[, varname.std]     + ppplusmcd.db[, varname.std]^2*(ppplusmcd.db$nsim - 1)         + ppplusmcd.db[, varname.avg]^2*ppplusmcd.db$nsim
    ppplusmcd.std.DB[, varname.std] = ppplusmcd.std.DB[, varname.std] + ppplusmcd.std.db[, varname.std]^2*(ppplusmcd.std.db$nsim - 1) + ppplusmcd.std.db[, varname.avg]^2*ppplusmcd.std.db$nsim
  }
}

for (j in 1:length(I.avg)) {
  varname.avg = I.avg[j]
  varname.std = I.std[j]
  
  detmcd.DB[, varname.std]        = sqrt(pmax(0, detmcd.DB[, varname.std]/(detmcd.DB$nsim - 1)               - detmcd.DB[, varname.avg]^2*detmcd.DB$nsim/(detmcd.DB$nsim - 1)))
  detmcd.std.DB[, varname.std]    = sqrt(pmax(0, detmcd.std.DB[, varname.std]/(detmcd.std.DB$nsim - 1)       - detmcd.std.DB[, varname.avg]^2*detmcd.std.DB$nsim/(detmcd.std.DB$nsim - 1)))
  fastmcd.DB[, varname.std]       = sqrt(pmax(0, fastmcd.DB[, varname.std]/(fastmcd.DB$nsim - 1)             - fastmcd.DB[, varname.avg]^2*fastmcd.DB$nsim/(fastmcd.DB$nsim - 1)))
  ppmcd.DB[, varname.std]         = sqrt(pmax(0, ppmcd.DB[, varname.std]/(ppmcd.DB$nsim - 1)                 - ppmcd.DB[, varname.avg]^2*ppmcd.DB$nsim/(ppmcd.DB$nsim - 1)))
  ppmcd.std.DB[, varname.std]     = sqrt(pmax(0, ppmcd.std.DB[, varname.std]/(ppmcd.std.DB$nsim - 1)         - ppmcd.std.DB[, varname.avg]^2*ppmcd.std.DB$nsim/(ppmcd.std.DB$nsim - 1)))
  ppplusmcd.DB[, varname.std]     = sqrt(pmax(0, ppplusmcd.DB[, varname.std]/(ppplusmcd.DB$nsim - 1)         - ppplusmcd.DB[, varname.avg]^2*ppplusmcd.DB$nsim/(ppplusmcd.DB$nsim - 1)))
  ppplusmcd.std.DB[, varname.std] = sqrt(pmax(0, ppplusmcd.std.DB[, varname.std]/(ppplusmcd.std.DB$nsim - 1) - ppplusmcd.std.DB[, varname.avg]^2*ppplusmcd.std.DB$nsim/(ppplusmcd.std.DB$nsim - 1)))
}

detmcd.db        = detmcd.DB
detmcd.std.db    = detmcd.std.DB
fastmcd.db       = fastmcd.DB
ppmcd.db         = ppmcd.DB
ppmcd.std.db     = ppmcd.std.DB
ppplusmcd.db     = ppplusmcd.DB
ppplusmcd.std.db = ppplusmcd.std.DB

rm(detmcd.DB, detmcd.std.DB, fastmcd.DB, ppmcd.DB, ppmcd.std.DB, ppplusmcd.DB, ppplusmcd.std.DB)
rm(i, j, varname.avg, varname.std, I.avg, I.std, varname, data.chunk)

# Mind the order!
save(fastmcd.db, detmcd.db, detmcd.std.db, ppmcd.db, ppmcd.std.db, ppplusmcd.db, ppplusmcd.std.db,
     file = "small.n.simulation.data.RData")

write.csv(detmcd.db,        file = "6pack.csv")
write.csv(detmcd.std.db,    file = "6pack.standardized.csv")
write.csv(ppplusmcd.db,     file = "6+1pack.csv")
write.csv(ppplusmcd.std.db, file = "6+1pack.standardized.csv")
