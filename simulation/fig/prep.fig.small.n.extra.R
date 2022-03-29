setwd("???")

load("small.n.simulation.data.RData")

rm(ppmcd.db, detmcd.db)

ppmcd.bd.db   = ppmcd.std.db
fastmcd.bd.db = fastmcd.db
detmcd.bd.db  = detmcd.std.db

rm(ppmcd.std.db, fastmcd.db, detmcd.std.db)

load("small.n.simulation.data.extra.RData")

mcd.nsamp = 500L

eps.array = c(0.0, 0.1, 0.2, 0.3, 0.4)
n.p.array = t(matrix(c(100, 100, 200, 400, 600,
                       2,   5,  10,  40,  60), nrow = 2, byrow = TRUE))
delta.array = c(1.0, 0.001)
psi.array   = c(5.0)

r.range <- function(p, delta = 1.0, alpha = 0.01) {
  r.min = ceiling(1.2*sqrt(qchisq(1 - alpha, df = p)/p) + delta*sqrt(qchisq(1 - alpha, df = p)/p))
  ans = c(seq(r.min, 20, length.out = 7), seq(30, 250, by = 10))
}

## Selection frequency plots

rho.mult.array = c(0.75)

case.array = c("A", "B", "C", "D", "E")

legend.txt = list()
for (i.n.p in 1:nrow(n.p.array)) {
  n = n.p.array[i.n.p, 1]
  p = n.p.array[i.n.p, 2]
  legend.txt = c(legend.txt, bquote(n==.(n)~"and"~p==.(p)))
}

plot.ind = 0L

for (std.flag in c(FALSE, TRUE))
for (n.WS in c(6L, 7L)) {
  for (rho.mult in rho.mult.array) {
    dir.name = "small.sample.extra/"
    if (!dir.exists(dir.name))
      dir.create(dir.name)
    
    dir.name = paste("small.sample.extra/rho.mult=", zapsmall(rho.mult, digits = 3), sep = "")
    
    if (!dir.exists(dir.name))
      dir.create(dir.name)
    
    dir.name = file.path(getwd(), dir.name)
    
    for (i.eps in 1:length(eps.array))
    for (i.delta in 1:length(delta.array)) {
      plot.ind = plot.ind + 1L
        
      n = n.p.array[i.n.p, 1]
      p = n.p.array[i.n.p, 2]
        
      file.name = paste(if (n.WS == 6L) "6pack." else "6+1pack.",
                        if (std.flag) "sphered." else "unsphered.",
                        "WS.plots.", "delta=", delta.array[i.delta], ".eps=", eps.array[i.eps], sep = "")
        
      grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 8, height = 5)
        
      WS.tab = NULL
      WS.avg.names = sapply(1:n.WS, function(i) paste("WS", i, ".avg", sep = ""))
      
      ppmcd.breakdown  = rep(FALSE, nrow(n.p.array))
      detmcd.breakdown = rep(FALSE, nrow(n.p.array))
        
      for (i.n.p in 1:nrow(n.p.array)) {
        if ((std.flag == FALSE) && (n.WS == 6L))
          mcd.db = detmcd.db
        else if ((std.flag == TRUE) && (n.WS == 6L))
          mcd.db = detmcd.std.db
        else if ((std.flag == FALSE) && (n.WS == 7L))
          mcd.db = ppplusmcd.db
        else if ((std.flag == TRUE) && (n.WS == 7L))
          mcd.db = ppplusmcd.std.db
        
        row.entries = which((mcd.db[, "rho.mult"] == rho.mult) &
                            (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                            (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta]))
        
        prob = colMeans(mcd.db[row.entries, WS.avg.names])
          
        WS.tab = rbind(WS.tab, prob)
        
        if (std.flag) {
          mcd.db = ppmcd.std.db
          ppmcd.r.max = max(mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                         (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                         (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta])), "r"])
          mcd.db = ppmcd.bd.db
          ppmcd.breakdown[i.n.p] = mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                                (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                                (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta]) &
                                                (mcd.db[, "r"] ==  ppmcd.r.max)), "breakdown.chance"] > 0.01
          
          mcd.db = detmcd.std.db
          detmcd.r.max = max(mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                          (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                          (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta])), "r"])
          mcd.db = detmcd.bd.db
          detmcd.breakdown[i.n.p] = mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                                 (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                                 (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta]) &
                                                 (mcd.db[, "r"] ==  detmcd.r.max)), "breakdown.chance"] > 0.01
        }
      }
      
      colnames(WS.tab) = 1:n.WS
      rownames(WS.tab) = letters[1:nrow(n.p.array)]
        
      WS.tab = as.table(WS.tab)
        
      cols = NULL
      
      for (i.WS in 1:n.WS) {
        bd = if (i.WS < 7) detmcd.breakdown else ppmcd.breakdown
          
        alphas = ifelse(bd, 0.10, 1.00)
        
        col.vec = c(rgb(0.0000, 0.4470, 0.7410, alphas[1]),
                    rgb(0.8500, 0.3250, 0.0980, alphas[2]),
                    rgb(0.9290, 0.6940, 0.1250, alphas[3]),
                    rgb(0.4940, 0.1840, 0.5560, alphas[4]),
                    rgb(0.4660, 0.6740, 0.1880, alphas[5]))
        
        cols = cbind(cols, col.vec)
      }
      
      col.vec = c(rgb(0.0000, 0.4470, 0.7410),
                  rgb(0.8500, 0.3250, 0.0980),
                  rgb(0.9290, 0.6940, 0.1250),
                  rgb(0.4940, 0.1840, 0.5560),
                  rgb(0.4660, 0.6740, 0.1880))
      
      par(mar = c(4, 4, 3, 3))
      barplot(WS.tab, col = cols, beside = TRUE, ylim = c(0, max(WS.tab, 0.6)),
              xlab = "Warmstart index", ylab = "Relative frequency optimal")
      box(which = "plot", lty = "solid")
      legend("topleft", legend = as.expression(legend.txt), fill = col.vec)
      
      grDevices::dev.off()
    }
  }
}

# Breakdown chance plots

rho.mult.array = c(0.75, sqrt(0.9))

delta.array = rev(delta.array)

legend.txt = c("PP MCD", "Fast MCD", "DetMCD (sphered)")

plot.ind = 0L

horiz.flag = FALSE

for (rho.mult in rho.mult.array) {
  dir.name = "small.sample.extra/"
  if (!dir.exists(dir.name))
    dir.create(dir.name)
      
  dir.name = paste("small.sample.extra/rho.mult=", zapsmall(rho.mult, digits = 3), sep = "")
      
  if (!dir.exists(dir.name))
    dir.create(dir.name)
      
  dir.name = file.path(getwd(), dir.name)
  
  file.name = paste("breakdown.chance.rho.mult=", zapsmall(rho.mult, digits = 3), sep = "")
  
  if (horiz.flag)
    grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 6)
  else
    grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 10, height = 12)
  
  par(mfcol = if (horiz.flag) c(2, 4) else c(4, 2))
  par(mar = c(5, 5, 2, 2))
  
  for (i.delta in 1:length(delta.array))
  for (i.eps in 2:length(eps.array)) {
    plot.ind = plot.ind + 1L
    
    eps   = eps.array[i.eps]
    delta = delta.array[i.delta]
    
    if (horiz.flag)
      par(mfg = c(i.delta, i.eps - 1L))
    else
      par(mfg = c(i.eps - 1, i.delta))
      
    ppmcd.bdp   = rep(0.0, nrow(n.p.array))
    fastmcd.bdp = rep(0.0, nrow(n.p.array))
    detmcd.bdp  = rep(0.0, nrow(n.p.array))
      
    for (i.n.p in 1:nrow(n.p.array)) {
      n = n.p.array[i.n.p, 1]
      p = n.p.array[i.n.p, 2]
        
      mcd.db = ppmcd.bd.db
      ppmcd.r.max = max(mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                     (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                     (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta])), "r"])
      ppmcd.bdp[i.n.p] = mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                      (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                      (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta]) &
                                      (mcd.db[, "r"] ==  ppmcd.r.max)), "breakdown.chance"]
        
      mcd.db = fastmcd.bd.db
      fastmcd.r.max = max(mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                       (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                       (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta])), "r"])
      fastmcd.bdp[i.n.p] = mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                        (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                        (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta]) &
                                        (mcd.db[, "r"] ==  ppmcd.r.max)), "breakdown.chance"]
      
      mcd.db = detmcd.bd.db
      detmcd.r.max = max(mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                      (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                      (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta])), "r"])
      detmcd.bdp[i.n.p] = mcd.db[which((mcd.db[, "rho.mult"] == rho.mult) &
                                       (mcd.db[, "n"] == n.p.array[i.n.p, 1]) & (mcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                       (mcd.db[, "eps"] == eps.array[i.eps]) & (mcd.db[, "delta"] == delta.array[i.delta]) &
                                       (mcd.db[, "r"] ==  detmcd.r.max)), "breakdown.chance"]
    }
      
    bdp.table = rbind(ppmcd.bdp, fastmcd.bdp, detmcd.bdp)
    colnames(bdp.table) = LETTERS[1:nrow(n.p.array)]
    rownames(bdp.table) = c("PP", "Fast", "Det")
      
    col.vec = c(rgb(0.0, 0.0, 0.0, 0.5), rgb(0.0, 0.0, 1.0, 0.5), rgb(1.0, 0.0, 0.0, 0.5))
    
    cont.txt = paste("of", ifelse(delta == 1.0, "cluster", "point"), "contamination")
    
    ind.symb = paste("(", letters[plot.ind], ")", sep = "")
    
    barplot(bdp.table, beside = TRUE, ylim = c(0, 1.4), col = col.vec,
            xlab = "Scenario", ylab = "Breakdown chance",
            main = bquote(.(ind.symb)~epsilon==.(eps)~.(cont.txt)))
    box(which = "plot", lty = "solid")
    legend("topleft", legend = legend.txt, fill = col.vec)
  }
    
  grDevices::dev.off()
}

delta.array = rev(delta.array)

try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
try(dev.off(), silent = TRUE)
