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
y.labels   = c(bquote("avg." ~ e[mu]), bquote("avg." ~ e[Sigma]), bquote("avg." ~ e[Sigma]^2),
               bquote("breakdown chance"), bquote("avg. bulk impurity"), bquote("avg. # of outliers in bulk"), bquote("avg. bulk size"), bquote("avg. run time (sec)"), bquote("avg. objective (covariance determinant)")) #bquote(e[Sigma]^3), bquote(e[Sigma]^4))

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

    grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 6)

    row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) & (fastmcd.db[, "p"] == p) & (fastmcd.db[, "r"] == r) &
                        (fastmcd.db[, "eps"] == eps) & (fastmcd.db[, "delta"] == delta))

    # Error measures
    fastmcd.err    = fastmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.err     = detmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.std.err = detmcd.std.db[row.entries, stat.types[i.stat.type]]
    ppmcd.err      = ppmcd.db[row.entries, stat.types[i.stat.type]]
    ppmcd.std.err  = ppmcd.std.db[row.entries, stat.types[i.stat.type]]

    cont.type = if (delta == 1.0) "cluster contamination" else "point contamination"

    min.y = min(fastmcd.err, detmcd.std.err, ppmcd.std.err)
    max.y = max(fastmcd.err, detmcd.std.err, ppmcd.std.err)

    max.y = min.y + 1.40*(max.y - min.y)

    par(mar = c(4, 4, 3, 3))
    y.lab = as.expression(y.labels[i.stat.type])
    plot(x = NULL, y = NULL, xlim = c(min(n.array), max(n.array)), ylim = c(min.y, max.y),
         xlab = bquote(n), ylab = bquote(.(y.labels[i.stat.type])),
         main = bquote(p~"="~.(p)*":"~.(cont.type)~"with"~epsilon~"="~.(eps)*","~r~"="~.(r)))

    lines(n.array, ppmcd.std.err,  lty = 1, lwd = 1, col = "black")
    lines(n.array, detmcd.std.err, lty = 2, lwd = , col = "red")
    lines(n.array, fastmcd.err,    lty = 4, lwd = 1, col = "blue")

    points(n.array, ppmcd.std.err,  lty = 1, lwd = 1, pch = 1, col = "black")
    points(n.array, detmcd.std.err, lty = 2, lwd = 1, pch = 2, col = "red")
    points(n.array, fastmcd.err,    lty = 4, lwd = 1, pch = 0, col = "blue")

    legend("topleft", legend = c("PP MCD", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
           lty = c(1, 2, 4), lwd = c(1, 1, 1), pch = c(1, 2, 0), col = c("black", "red", "blue"))

    grDevices::dev.off()
  }
}

plot.ind = 0L
for (rho.mult in rho.mult.array)
{
  dir.name = "large.sample/"
  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = paste("large.sample/rho.mult=", zapsmall(rho.mult, digits = 3), ".extra", sep = "")

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

    grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 6)

    row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) & (fastmcd.db[, "p"] == p) & (fastmcd.db[, "r"] == r) &
                        (fastmcd.db[, "eps"] == eps) & (fastmcd.db[, "delta"] == delta))

    # Error measures
    fastmcd.err    = fastmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.err     = detmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.std.err = detmcd.std.db[row.entries, stat.types[i.stat.type]]
    ppmcd.err      = ppmcd.db[row.entries, stat.types[i.stat.type]]
    ppmcd.std.err  = ppmcd.std.db[row.entries, stat.types[i.stat.type]]

    cont.type = if (delta == 1.0) "cluster contamination" else "point contamination"

    min.y = min(detmcd.err, detmcd.std.err)
    max.y = max(detmcd.err, detmcd.std.err)

    max.y = min.y + 1.40*(max.y - min.y)

    par(mar = c(4, 4, 3, 3))
    y.lab = as.expression(y.labels[i.stat.type])
    plot(x = NULL, y = NULL, xlim = c(min(n.array), max(n.array)), ylim = c(min.y, max.y),
         xlab = bquote(n), ylab = bquote(.(y.labels[i.stat.type])),
         main = bquote(p~"="~.(p)*":"~.(cont.type)~"with"~epsilon~"="~.(eps)*","~r~"="~.(r)))

    lines(n.array, detmcd.err,     lty = 2, lwd = 1, col = "red")
    lines(n.array, detmcd.std.err, lty = 2, lwd = , col = "red")

    lines(n.array,  detmcd.err,     lty = 2, lwd = 1, col = "red")
    points(n.array, detmcd.std.err, lty = 2, lwd = 1, pch = 2, col = "red")

    legend("topleft", legend = c("DetMCD (sphered)", "DetMCD (not sphered)"),
           lty = c(2, 2), lwd = c(1, 2), pch = c(2, NA), col = c("red", "red"))

    grDevices::dev.off()
  }
}

##

conf = rbind(c(0.75, 10, 100.0, 0.4, 1.0),
             c(0.75, 60, 100.0, 0.4, 1.0))
colnames(conf) = c("rho.mult", "p", "r", "eps", "delta")

dir.name = "large.sample/"
if (!dir.exists(dir.name))
  dir.create(dir.name)

dir.name = "large.sample"
dir.name = file.path(getwd(), dir.name)

file.name = "large.sample.figure"
grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 10, height = 7)

par(mfrow = c(2, 2))

#for (i.stat.type in c(1, 2, 8))

ind = 0

for (i.stat.type in c(2, 8))
for (i in 1:nrow(conf)) {
  ind = ind + 1L

  row.entries = which((fastmcd.db[, "rho.mult"] == conf[i, "rho.mult"]) &
                        (fastmcd.db[, "r"] == conf[i, "r"]) & (fastmcd.db[, "p"] == conf[i, "p"]) &
                        (fastmcd.db[, "eps"] == conf[i, "eps"]) & (fastmcd.db[, "delta"] == conf[i, "delta"]))

  r.array = fastmcd.db[row.entries, "r"]

  # Error measures
  fastmcd.err    = fastmcd.db[row.entries, stat.types[i.stat.type]]
  detmcd.err     = detmcd.db[row.entries, stat.types[i.stat.type]]
  detmcd.std.err = detmcd.std.db[row.entries, stat.types[i.stat.type]]
  ppmcd.err      = ppmcd.db[row.entries, stat.types[i.stat.type]]
  ppmcd.std.err  = ppmcd.std.db[row.entries, stat.types[i.stat.type]]

  if (sum(is.na(fastmcd.err)) > 0)
    print(paste(stat.types[i.stat.type], ".delta=", delta.array[i.delta], ".eps=", eps.array[i.eps], ".case=", case.array[i.n.p],
                ".n=", n.p.array[i.n.p, 1], ".p=", n.p.array[i.n.p, 2], sep = ""))

  cont.type = if (conf[i, "delta"] == 1) "cluster contamination" else "point contamination"

  min.y = min(fastmcd.err, detmcd.std.err, ppmcd.std.err)
  max.y = max(fastmcd.err, detmcd.std.err, ppmcd.std.err)

  max.y = min.y + 1.6*(max.y - min.y)

  par(mar = c(4, 4, 3, 3))
  y.lab = as.expression(y.labels[i.stat.type])
  plot(x = NULL, y = NULL, xlim = c(min(n.array), max(n.array)), ylim = c(min.y, max.y),
       xlab = bquote(n), ylab = bquote(.(y.labels[i.stat.type])),
       main = bquote("("*.(letters[ind])*")"~p~"="~.(conf[i, "p"])*":"~.(cont.type)~"with"~epsilon~"="~.(conf[i, "eps"])*","~r~"="~.(conf[i, "r"])))

  lines(n.array, ppmcd.std.err,  lty = 1, lwd = 1, col = "black")
  lines(n.array, detmcd.std.err, lty = 2, lwd = , col = "red")
  lines(n.array, fastmcd.err,    lty = 4, lwd = 1, col = "blue")

  points(n.array, ppmcd.std.err,  lty = 1, lwd = 1, pch = 1, col = "black")
  points(n.array, detmcd.std.err, lty = 2, lwd = 1, pch = 2, col = "red")
  points(n.array, fastmcd.err,    lty = 4, lwd = 1, pch = 0, col = "blue")

  legend("topleft", legend = c("PP MCD", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
         lty = c(1, 2, 4), lwd = c(1, 1, 1), pch = c(1, 2, 0), col = c("black", "red", "blue"))
}

grDevices::dev.off()
