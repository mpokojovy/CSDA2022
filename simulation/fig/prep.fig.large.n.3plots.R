setwd("???")

load("large.n.simulation.data.RData")

mcd.nsamp = 500L

rho.mult.array = c(0.75, sqrt(0.9))

eps.array   = c(0.1, 0.4)
n.array     = c(1, 3, 5, 10, 15, 20, 30, 40)*1000
p.array     = c(10, 60)
delta.array = c(1.0)
r.array     = c(100.0, 200.0)

## Do plots
case.array = c("A", "B", "C", "D", "E")
stat.types = c("err.mu.avg", "err.Sigma1.avg", "err.Sigma2.avg", "breakdown.chance", "impurity.avg", "n.outlier.avg", "bulk.size.avg", "time.avg", "obj.avg") #"err.Sigma3.avg", "err.Sigma4.avg")
y.labels   = c(bquote(e[mu]), bquote(e[Sigma]), bquote(e[Sigma]^2),
               bquote("breakdown chance"), bquote("avg. rew. bulk impurity"), bquote("avg. # of outliers in rew. bulk"), bquote("avg. rew. bulk size"), bquote("avg. run time (sec)"), bquote("avg. objective (covariance determinant)")) #bquote(e[Sigma]^3), bquote(e[Sigma]^4))

plot.ind = 0L
for (rho.mult in rho.mult.array)
{
  dir.name = "large.sample/"
  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = paste("large.sample/rho.mult=", zapsmall(rho.mult, digits = 3), sep = "")

  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = file.path(getwd(), dir.name)

  for (p in p.array)
  for (r in r.array)
  for (delta in delta.array)
  for (eps in eps.array)
  for (i.stat.type in 1:length(stat.types))
  {
    plot.ind = plot.ind + 1L

    file.name = paste(stat.types[i.stat.type], ".p=", p, ".delta=", delta, ".eps=", eps, ".r=", r, sep = "")

    grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 8, height = 6)
    
    #grDevices::setEPS()
    #grDevices::postscript(file = paste(dir.name, "/", file.name, ".eps", sep = ""), horiz = FALSE, onefile = TRUE, width = 8, height = 8, paper = "special")
    #grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), paper = "letter")

    row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) & (fastmcd.db[, "p"] == p) & (fastmcd.db[, "r"] == r) &
                        (fastmcd.db[, "eps"] == eps) & (fastmcd.db[, "delta"] == delta))

    # Error measures
    fastmcd.err    = fastmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.err     = detmcd.std.db[row.entries, stat.types[i.stat.type]]
    ppmcd.err      = ppmcd.std.db[row.entries, stat.types[i.stat.type]]
    # if (sum(is.na(fastmcd.err)) > 0)
    #   print(paste(stat.types[i.stat.type], ".delta=", delta.array[i.delta], ".eps=", eps.array[i.eps], ".case=", case.array[i.n.p],
    #               ".n=", n.p.array[i.n.p, 1], ".p=", n.p.array[i.n.p, 2], sep = ""))

    cont.type = if (delta == 1.0) "cluster contamination" else "point contamination"
    
    min.y = min(fastmcd.err, detmcd.err, ppmcd.err)
    max.y = max(fastmcd.err, detmcd.err, ppmcd.err)
    
    max.y = min.y + 1.35*(max.y - min.y)
    
    par(mar = c(5, 5, 5, 5))
    y.lab = as.expression(y.labels[i.stat.type])
    plot(x = NULL, y = NULL, xlim = c(min(n.array), max(n.array)), ylim = c(min.y, max.y),
         xlab = bquote(n), ylab = bquote(.(y.labels[i.stat.type])),
         main = bquote(p~"="~.(p)*":"~.(cont.type)~"with"~epsilon~"="~.(eps)*","~r~"="~.(r)))
    
    lines(n.array, ppmcd.err,      lty = 1, lwd = 2, col = "black")
    lines(n.array, fastmcd.err,    lty = 2, lwd = 2, col = "blue")
    lines(n.array, detmcd.err,     lty = 4, lwd = 2, col = "red")
    
    legend("topleft", legend = c("PP MCD", paste("FastMCD(", mcd.nsamp, ")", sep = ""), "DetMCD"),
           lty = c(1, 2, 4), lwd = c(2, 2, 2), col = c("black", "blue", "red"))
    
    # cont.type = if (delta.array[i.delta] == 1) "cluster contamination" else "point contamination"
    # 
    # min.y = min(fastmcd.err, detmcd.err, detmcd.std.err, ppmcd.err, ppmcd.std.err)
    # max.y = max(fastmcd.err, detmcd.err, detmcd.std.err, ppmcd.err, ppmcd.std.err)
    # 
    # max.y = min.y + 1.35*(max.y - min.y)
    # 
    # par(mar = c(5, 5, 5, 5))
    # y.lab = as.expression(y.labels[i.stat.type])
    # plot(x = NULL, y = NULL, xlim = c(min(r.array), max(r.array)), ylim = c(min.y, max.y),
    #      xlab = bquote("r"), ylab = bquote(.(y.labels[i.stat.type])),
    #      main = bquote("Setting" ~ .(case.array[i.n.p])*" (n ="~.(n.p.array[i.n.p, 1])*", p ="~.(n.p.array[i.n.p, 2])*"): "~.(cont.type)~"with"~epsilon~"="~.(eps.array[i.eps])))
    # 
    # lines(r.array, ppmcd.err,      lty = 1, lwd = 1, col = "black")
    # lines(r.array, ppmcd.std.err,  lty = 1, lwd = 2, col = "black")
    # lines(r.array, detmcd.err,     lty = 2, lwd = 1, col = "red")
    # lines(r.array, detmcd.std.err, lty = 2, lwd = 2, col = "red")
    # lines(r.array, fastmcd.err,    lty = 4, lwd = 2, col = "blue")
    # 
    # legend("topleft", legend = c("PP MCD", "PP MCD (std'd)", "DetMCD", "DetMCD (std'd)", paste("FastMCD(", mcd.nsamp, ") (std'd)", sep = "")),
    #        lty = c(1, 1, 2, 2, 4), lwd = c(1, 2, 1, 2, 2), col = c("black", "black", "red", "red", "blue"))

    grDevices::dev.off()
  }
}