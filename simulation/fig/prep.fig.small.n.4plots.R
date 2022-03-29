setwd("???")

load("small.n.simulation.data.RData")

mcd.nsamp = 500L

rho.mult.array = c(0.75, sqrt(0.9))

eps.array = c(0.0, 0.1, 0.2, 0.3, 0.4)
n.p.array = t(matrix(c(100, 100, 200, 400, 600,
                       2,   5,  10,  40,  60), nrow = 2, byrow = TRUE))
delta.array = c(1.0, 0.001)
psi.array   = c(5.0)

r.range <- function(p, delta = 1.0, alpha = 0.01) {
  r.min = ceiling(1.2*sqrt(qchisq(1 - alpha, df = p)/p) + delta*sqrt(qchisq(1 - alpha, df = p)/p))
  ans = c(seq(r.min, 20, length.out = 7), seq(30, 250, by = 10))
}

## Do plots
case.array = c("A", "B", "C", "D", "E")
stat.types = c("err.mu.avg", "err.Sigma1.avg", "err.Sigma2.avg", "breakdown.chance", "impurity.avg", "n.outlier.avg", "bulk.size.avg", "time.avg", "obj.avg") #"err.Sigma3.avg", "err.Sigma4.avg")
y.labels   = c(bquote("avg." ~ e[mu]), bquote("avg." ~ e[Sigma]), bquote("avg." ~ e[Sigma]^2),
               bquote("breakdown chance"), bquote("avg. bulk impurity"), bquote("avg. # of outliers in bulk"), bquote("avg. bulk size"), bquote("avg. run time (sec)"), bquote("avg. objective (covariance determinant)")) #bquote(e[Sigma]^3), bquote(e[Sigma]^4))

plot.ind = 0L
for (rho.mult in rho.mult.array)
{
  dir.name = "small.sample/"
  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = paste("small.sample/rho.mult=", zapsmall(rho.mult, digits = 3), sep = "")

  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = file.path(getwd(), dir.name)

  for (i.n.p in 1:nrow(n.p.array))
  for (i.eps in 1:length(eps.array))
  for (i.delta in 1:length(delta.array))
  for (i.stat.type in 1:length(stat.types))
  {
    plot.ind = plot.ind + 1L

    file.name = paste(stat.types[i.stat.type], ".delta=", delta.array[i.delta], ".eps=", eps.array[i.eps], ".case=", case.array[i.n.p],
                      ".n=", n.p.array[i.n.p, 1], ".p=", n.p.array[i.n.p, 2], sep = "")

    grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 6)

    row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) &
                         (fastmcd.db[, "n"] == n.p.array[i.n.p, 1]) & (fastmcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                         (fastmcd.db[, "eps"] == eps.array[i.eps]) & (fastmcd.db[, "delta"] == delta.array[i.delta]))

    row.entries = which((detmcd.db[, "rho.mult"] == rho.mult) &
                          (detmcd.db[, "n"] == n.p.array[i.n.p, 1]) & (detmcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                          (detmcd.db[, "eps"] == eps.array[i.eps]) & (detmcd.db[, "delta"] == delta.array[i.delta]))

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

    cont.type = if (delta.array[i.delta] == 1) "cluster contamination" else "point contamination"

    min.y = min(fastmcd.err, detmcd.err, detmcd.std.err, ppmcd.std.err)
    max.y = max(fastmcd.err, detmcd.err, detmcd.std.err, ppmcd.std.err)

    max.y = min.y + 1.40*(max.y - min.y)

    par(mar = c(4, 4, 3, 3))
    y.lab = as.expression(y.labels[i.stat.type])
    plot(x = NULL, y = NULL, xlim = c(min(r.array), max(r.array)), ylim = c(min.y, max.y),
         xlab = bquote("r"), ylab = bquote(.(y.labels[i.stat.type])),
         main = bquote("Setting" ~ .(case.array[i.n.p])*" (n ="~.(n.p.array[i.n.p, 1])*", p ="~.(n.p.array[i.n.p, 2])*"): "~.(cont.type)~"with"~epsilon~"="~.(eps.array[i.eps])))

    lines(r.array, ppmcd.std.err,  lty = 1, lwd = 1, col = "black")
    lines(r.array, detmcd.std.err, lty = 2, lwd = , col = "red")
    lines(r.array, fastmcd.err,    lty = 4, lwd = 1, col = "blue")

    points(r.array, ppmcd.std.err,  lty = 1, lwd = 1, pch = 1, col = "black")
    points(r.array, detmcd.std.err, lty = 2, lwd = 1, pch = 2, col = "red")
    points(r.array, fastmcd.err,    lty = 4, lwd = 1, pch = 0, col = "blue")

    legend("topleft", legend = c("PP MCD", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
           lty = c(1, 2, 4), lwd = c(1, 1, 1), pch = c(1, 2, 0), col = c("black", "red", "blue"))

    grDevices::dev.off()
  }
}

plot.ind = 0L
for (rho.mult in rho.mult.array)
{
  dir.name = "small.sample/"
  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = paste("small.sample/rho.mult=", zapsmall(rho.mult, digits = 3), sep = "")

  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = file.path(getwd(), dir.name)

  for (i.eps in 1:length(eps.array))
  for (i.psi in 1:length(psi.array))
  for (i.stat.type in 1:length(stat.types))
  {
    plot.ind = plot.ind + 1L

    file.name = paste(stat.types[i.stat.type], ".psi=", psi.array[i.psi], ".eps=", eps.array[i.eps], sep = "")

    grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 6)

    row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) &
                          (fastmcd.db[, "eps"] == eps.array[i.eps]) & (fastmcd.db[, "psi"] == psi.array[i.psi]))

    # Error measures
    fastmcd.err    = fastmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.err     = detmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.std.err = detmcd.std.db[row.entries, stat.types[i.stat.type]]
    ppmcd.err      = ppmcd.db[row.entries, stat.types[i.stat.type]]
    ppmcd.std.err  = ppmcd.std.db[row.entries, stat.types[i.stat.type]]

    cont.type = "Radial contamination"

    AthruE = factor(LETTERS[1:5])

    min.y = min(fastmcd.err, detmcd.std.err, ppmcd.std.err)
    max.y = max(fastmcd.err, detmcd.std.err, ppmcd.std.err)

    max.y = min.y + 1.4*(max.y - min.y)

    par(mar = c(4, 4, 3, 3))
    y.lab = as.expression(y.labels[i.stat.type])

    plot(AthruE, y = NULL, ylim = c(min.y, max.y),
         xlab = bquote("Scenario"), xaxt = "n", ylab = bquote(.(y.labels[i.stat.type])),
         main = bquote(.(cont.type)~"with"~epsilon~"="~.(eps.array[i.eps])), type = "n")

    axis(side = 1, at = AthruE, labels = AthruE)

    lines(AthruE, ppmcd.std.err,  lty = 1, lwd = 1, col = "black")
    lines(AthruE, detmcd.std.err, lty = 2, lwd = , col = "red")
    lines(AthruE, fastmcd.err,    lty = 4, lwd = 1, col = "blue")

    points(AthruE, ppmcd.std.err,  lwd = 1, pch = 1, col = "black")
    points(AthruE, detmcd.std.err, lwd = 1, pch = 2, col = "red")
    points(AthruE, fastmcd.err,    lwd = 1, pch = 0, col = "blue")

    legend("topleft", legend = c("PP MCD", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
           lty = c(1, 2, 4), lwd = c(1, 1, 1), pch = c(1, 2, 0), col = c("black", "red", "blue"))

    grDevices::dev.off()
  }
}

# DetMCD plots

plot.ind = 0L
for (rho.mult in rho.mult.array)
{
  dir.name = "small.sample/"
  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = paste("small.sample/rho.mult=", zapsmall(rho.mult, digits = 3), ".extra", sep = "")

  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = file.path(getwd(), dir.name)

  for (i.n.p in 1:nrow(n.p.array))
    for (i.eps in 1:length(eps.array))
      for (i.delta in 1:length(delta.array))
        for (i.stat.type in 1:length(stat.types))
        {
          plot.ind = plot.ind + 1L

          file.name = paste(stat.types[i.stat.type], ".delta=", delta.array[i.delta], ".eps=", eps.array[i.eps], ".case=", case.array[i.n.p],
                            ".n=", n.p.array[i.n.p, 1], ".p=", n.p.array[i.n.p, 2], sep = "")

          grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 6)

          row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) &
                                (fastmcd.db[, "n"] == n.p.array[i.n.p, 1]) & (fastmcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                (fastmcd.db[, "eps"] == eps.array[i.eps]) & (fastmcd.db[, "delta"] == delta.array[i.delta]))

          row.entries = which((detmcd.db[, "rho.mult"] == rho.mult) &
                                (detmcd.db[, "n"] == n.p.array[i.n.p, 1]) & (detmcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                                (detmcd.db[, "eps"] == eps.array[i.eps]) & (detmcd.db[, "delta"] == delta.array[i.delta]))

          r.array = fastmcd.db[row.entries, "r"]

          # Error measures
          detmcd.err     = detmcd.db[row.entries, stat.types[i.stat.type]]
          detmcd.std.err = detmcd.std.db[row.entries, stat.types[i.stat.type]]

          if (sum(is.na(fastmcd.err)) > 0)
            print(paste(stat.types[i.stat.type], ".delta=", delta.array[i.delta], ".eps=", eps.array[i.eps], ".case=", case.array[i.n.p],
                        ".n=", n.p.array[i.n.p, 1], ".p=", n.p.array[i.n.p, 2], sep = ""))

          cont.type = if (delta.array[i.delta] == 1) "cluster contamination" else "point contamination"

          min.y = min(detmcd.err, detmcd.std.err)
          max.y = max(detmcd.err, detmcd.std.err)

          max.y = min.y + 1.40*(max.y - min.y)

          par(mar = c(4, 4, 3, 3))
          y.lab = as.expression(y.labels[i.stat.type])
          plot(x = NULL, y = NULL, xlim = c(min(r.array), max(r.array)), ylim = c(min.y, max.y),
               xlab = bquote("r"), ylab = bquote(.(y.labels[i.stat.type])),
               main = bquote("Setting" ~ .(case.array[i.n.p])*" (n ="~.(n.p.array[i.n.p, 1])*", p ="~.(n.p.array[i.n.p, 2])*"): "~.(cont.type)~"with"~epsilon~"="~.(eps.array[i.eps])))

          lines(r.array, detmcd.err,     lty = 2, lwd = 2, col = "red")
          lines(r.array, detmcd.std.err, lty = 2, lwd = 1, col = "red")

          #lines(r.array,  detmcd.err,     lty = 2, lwd = 1, col = "red")
          points(r.array, detmcd.std.err, lty = 2, lwd = 1, pch = 2, col = "red")

          legend("topleft", legend = c("DetMCD (sphered)", "DetMCD (not sphered)"),
                 lty = c(2, 2), lwd = c(1, 2), pch = c(2, NA), col = c("red", "red"))

          grDevices::dev.off()
        }
}

plot.ind = 0L
for (rho.mult in rho.mult.array)
{
  dir.name = "small.sample/"
  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = paste("small.sample/rho.mult=", zapsmall(rho.mult, digits = 3), ".extra", sep = "")

  if (!dir.exists(dir.name))
    dir.create(dir.name)

  dir.name = file.path(getwd(), dir.name)

  for (i.eps in 1:length(eps.array))
  for (i.psi in 1:length(psi.array))
  for (i.stat.type in 1:length(stat.types))
  {
    plot.ind = plot.ind + 1L

    file.name = paste(stat.types[i.stat.type], ".psi=", psi.array[i.psi], ".eps=", eps.array[i.eps], sep = "")

    grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 6)

    row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) &
                          (fastmcd.db[, "eps"] == eps.array[i.eps]) & (fastmcd.db[, "psi"] == psi.array[i.psi]))

    # Error measures
    fastmcd.err    = fastmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.err     = detmcd.db[row.entries, stat.types[i.stat.type]]
    detmcd.std.err = detmcd.std.db[row.entries, stat.types[i.stat.type]]
    ppmcd.err      = ppmcd.db[row.entries, stat.types[i.stat.type]]
    ppmcd.std.err  = ppmcd.std.db[row.entries, stat.types[i.stat.type]]

    cont.type = "Radial contamination"

    AthruE = factor(LETTERS[1:5])

    min.y = min(detmcd.err, detmcd.std.err)
    max.y = max(detmcd.err, detmcd.std.err)

    max.y = min.y + 1.4*(max.y - min.y)

    par(mar = c(4, 4, 3, 3))
    y.lab = as.expression(y.labels[i.stat.type])

    plot(AthruE, y = NULL, ylim = c(min.y, max.y),
         xlab = bquote("Scenario"), xaxt = "n", ylab = bquote(.(y.labels[i.stat.type])),
         main = bquote(.(cont.type)~"with"~epsilon~"="~.(eps.array[i.eps])), type = "n")

    axis(side = 1, at = AthruE, labels = AthruE)

    lines(AthruE, detmcd.err,     lty = 2, lwd = 2, col = "red")
    lines(AthruE, detmcd.std.err, lty = 2, lwd = 1, col = "red")

    points(AthruE, detmcd.err,     lty = 2, lwd = 1, pch = NA, col = "red")
    points(AthruE, detmcd.std.err, lty = 2, lwd = 1, pch = 2,  col = "red")

    legend("topleft", legend = c("DetMCD (sphered)", "DetMCD (not sphered)"),
          lty = c(2, 2), lwd = c(1, 2), pch = c(2, NA), col = c("red", "red"))

    grDevices::dev.off()
  }
}

##

#i.stat.type = 2
i.stat.type = 1

conf = rbind(c(1, 0.75, 100,  2, 0.1, 0.001),
             c(1, 0.75, 100,  2, 0.4, 0.001),
             c(2, 0.75, 100,  5, 0.1, 0.001),
             c(2, 0.75, 100,  5, 0.4, 0.001),
             c(3, 0.75, 200, 10, 0.1, 0.001),
             c(3, 0.75, 200, 10, 0.4, 0.001),
             c(4, 0.75, 400, 40, 0.1, 0.001),
             c(4, 0.75, 400, 40, 0.4, 0.001),
             c(5, 0.75, 600, 60, 0.1, 0.001),
             c(5, 0.75, 600, 60, 0.4, 0.001))

# conf = rbind(c(1, 0.75, 100,  2, 0.1, 1),
#              c(1, 0.75, 100,  2, 0.4, 1),
#              c(2, 0.75, 100,  5, 0.1, 1),
#              c(2, 0.75, 100,  5, 0.4, 1),
#              c(3, 0.75, 200, 10, 0.1, 1),
#              c(3, 0.75, 200, 10, 0.4, 1),
#              c(4, 0.75, 400, 40, 0.1, 1),
#              c(4, 0.75, 400, 40, 0.4, 1),
#              c(5, 0.75, 600, 60, 0.1, 1),
#              c(5, 0.75, 600, 60, 0.4, 1))

# conf = rbind(#c(1, 0.75, 100,  2, 0.1, 0.001),
#              #c(1, 0.75, 100,  2, 0.4, 0.001),
#              c(2, 0.75, 100,  5, 0.1, 0.001),
#              c(2, 0.75, 100,  5, 0.4, 0.001),
#              #c(3, 0.75, 200, 10, 0.1, 0.001),
#              #c(3, 0.75, 200, 10, 0.4, 0.001),
#              #c(4, 0.75, 400, 40, 0.1, 0.001),
#              #c(4, 0.75, 400, 40, 0.4, 0.001),
#              c(5, 0.75, 600, 60, 0.1, 0.001),
#              c(5, 0.75, 600, 60, 0.4, 0.001))
# conf = rbind(#c(1, 0.75, 100,  2, 0.1, 1),
#              #c(1, 0.75, 100,  2, 0.4, 1),
#              c(2, 0.75, 100,  5, 0.1, 1),
#              c(2, 0.75, 100,  5, 0.4, 1),
#              #c(3, 0.75, 200, 10, 0.1, 1),
#              #c(3, 0.75, 200, 10, 0.4, 1),
#              #c(4, 0.75, 400, 40, 0.1, 1),
#              #c(4, 0.75, 400, 40, 0.4, 1),
#              c(5, 0.75, 600, 60, 0.1, 1),
#              c(5, 0.75, 600, 60, 0.4, 1))
# conf = rbind(c(5, 0.75, 600, 60, 0.1, 0.001),
#              c(5, 0.75, 600, 60, 0.4, 0.001))

colnames(conf) = c("setting", "rho.mult", "n", "p", "eps", "delta")

dir.name = "small.sample/"
if (!dir.exists(dir.name))
  dir.create(dir.name)

dir.name = "small.sample"
dir.name = file.path(getwd(), dir.name)

file.name = "small.sample.figure"
grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 10, height = 20) #width = 12, height = 20)
#grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 4) #width = 12, height = 20)

par(mfrow = c(5, 2))

for (i in 1:nrow(conf)) {
  row.entries = which((fastmcd.db[, "rho.mult"] == conf[i, "rho.mult"]) &
                        (fastmcd.db[, "n"] == conf[i, "n"]) & (fastmcd.db[, "p"] == conf[i, "p"]) &
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

  min.y = min(fastmcd.err, detmcd.err, detmcd.std.err, ppmcd.err, ppmcd.std.err)
  max.y = max(fastmcd.err, detmcd.err, detmcd.std.err, ppmcd.err, ppmcd.std.err)

  max.y = min.y + 1.6*(max.y - min.y)

  par(mar = c(4, 4, 3, 3))
  y.lab = as.expression(y.labels[i.stat.type])
  plot(x = NULL, y = NULL, xlim = c(min(r.array), max(r.array)), ylim = c(min.y, max.y),
       xlab = bquote("r"), ylab = bquote(.(y.labels[i.stat.type])),
       main = bquote("(" * .(letters[i]) * ")" ~ "Setting" ~ .(LETTERS[conf[i, "setting"]])*" (n ="~.(conf[i, "n"])*", p ="~.(conf[i, "p"])*"): "~.(cont.type)~"with"~epsilon~"="~.(conf[i, "eps"])))
       #sub = paste("(", letters[i], ")", sep = ""))

  lines(r.array, ppmcd.std.err,  lty = 1, lwd = 1, col = "black")
  #lines(r.array, detmcd.err,     lty = 2, lwd = 1, col = "red")
  lines(r.array, detmcd.std.err, lty = 2, lwd = , col = "red")
  lines(r.array, fastmcd.err,    lty = 4, lwd = 1, col = "blue")

  points(r.array, ppmcd.std.err,  lty = 1, lwd = 1, pch = 1, col = "black")
  #lines(r.array, detmcd.err,     lty = 2, lwd = 1, col = "red")
  points(r.array, detmcd.std.err, lty = 2, lwd = 1, pch = 2, col = "red")
  points(r.array, fastmcd.err,    lty = 4, lwd = 1, pch = 0, col = "blue")

  #legend("topleft", legend = c("PP MCD", "DetMCD (not sphered)", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
  #       lty = c(1, 2, 2, 4), lwd = c(2, 1, 2, 2), col = c("black", "red", "red", "blue"))

  legend("topleft", legend = c("PP MCD", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
         lty = c(1, 2, 4), lwd = c(1, 1, 1), pch = c(1, 2, 0), col = c("black", "red", "blue"))
}

grDevices::dev.off()

# for (i in 1:nrow(conf)) {
#   row.entries = which((fastmcd.db[, "rho.mult"] == conf[i, "rho.mult"]) &
#                         (fastmcd.db[, "n"] == conf[i, "n"]) & (fastmcd.db[, "p"] == conf[i, "p"]) &
#                         (fastmcd.db[, "eps"] == conf[i, "eps"]) & (fastmcd.db[, "delta"] == conf[i, "delta"]))
#
#   r.array = fastmcd.db[row.entries, "r"]
#
#   # Error measures
#   fastmcd.err    = fastmcd.db[row.entries, stat.types[i.stat.type]]
#   detmcd.err     = detmcd.db[row.entries, stat.types[i.stat.type]]
#   detmcd.std.err = detmcd.std.db[row.entries, stat.types[i.stat.type]]
#   ppmcd.err      = ppmcd.db[row.entries, stat.types[i.stat.type]]
#   ppmcd.std.err  = ppmcd.std.db[row.entries, stat.types[i.stat.type]]
#
#   if (sum(is.na(fastmcd.err)) > 0)
#     print(paste(stat.types[i.stat.type], ".delta=", delta.array[i.delta], ".eps=", eps.array[i.eps], ".case=", case.array[i.n.p],
#                 ".n=", n.p.array[i.n.p, 1], ".p=", n.p.array[i.n.p, 2], sep = ""))
#
#   cont.type = if (conf[i, "delta"] == 1) "cluster contamination" else "point contamination"
#
#   min.y = min(detmcd.err, detmcd.std.err)
#   max.y = max(detmcd.err, detmcd.std.err)
#
#   max.y = min.y + 1.6*(max.y - min.y)
#
#   par(mar = c(4, 4, 3, 3))
#   y.lab = as.expression(y.labels[i.stat.type])
#   plot(x = NULL, y = NULL, xlim = c(min(r.array), max(r.array)), ylim = c(min.y, max.y),
#        xlab = bquote("r"), ylab = bquote(.(y.labels[i.stat.type])),
#        main = bquote("(" * .(letters[i]) * ")" ~ "Setting" ~ .(LETTERS[conf[i, "setting"]])*" (n ="~.(conf[i, "n"])*", p ="~.(conf[i, "p"])*"): "~.(cont.type)~"with"~epsilon~"="~.(conf[i, "eps"])))
#        #sub = paste("(", letters[i], ")", sep = ""))
#
#   lines(r.array, detmcd.err,     lty = 2, lwd = 2, col = "red")
#   lines(r.array, detmcd.std.err, lty = 2, lwd = 1, col = "red")
#
#   points(r.array, detmcd.std.err, lty = 2, lwd = 1, pch = 2, col = "red")
#
#   #legend("topleft", legend = c("PP MCD", "DetMCD (not sphered)", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
#   #       lty = c(1, 2, 2, 4), lwd = c(2, 1, 2, 2), col = c("black", "red", "red", "blue"))
#
#   legend("topleft", legend = c("DetMCD (sphered)", "DetMCD (not sphered)"),
#          lty = c(2, 2), lwd = c(1, 2), pch = c(2, NA), col = c("red", "red"))
# }
#
# grDevices::dev.off()

## Radial

i.stat.type = 2

dir.name = "small.sample/"
if (!dir.exists(dir.name))
  dir.create(dir.name)

dir.name = "small.sample"
dir.name = file.path(getwd(), dir.name)

file.name = "small.sample.radial.figure"
#grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 10, height = 3.5) #width = 12, height = 20)
grDevices::pdf(file = paste(dir.name, "/", file.name, ".pdf", sep = ""), width = 12, height = 4) #width = 12, height = 20)

par(mfrow = c(1, 2))

#for (i.n.p in 1:nrow(n.p.array))
i = 0L
for (eps in c(0.1, 0.4))
for (rho.mult in c(0.75))
for (psi in c(5.0)) {
  i = i + 1L
  row.entries = which((fastmcd.db[, "rho.mult"] == rho.mult) &
                      #(fastmcd.db[, "n"] == n.p.array[i.n.p, 1]) & (fastmcd.db[, "p"] == n.p.array[i.n.p, 2]) &
                      (fastmcd.db[, "eps"] == eps) & (fastmcd.db[, "psi"] == psi))

  # Error measures
  fastmcd.err    = fastmcd.db[row.entries, stat.types[i.stat.type]]
  detmcd.err     = detmcd.db[row.entries, stat.types[i.stat.type]]
  detmcd.std.err = detmcd.std.db[row.entries, stat.types[i.stat.type]]
  ppmcd.err      = ppmcd.db[row.entries, stat.types[i.stat.type]]
  ppmcd.std.err  = ppmcd.std.db[row.entries, stat.types[i.stat.type]]

  cont.type = "Radial contamination"

  AthruE = factor(LETTERS[1:5])

  min.y = min(fastmcd.err, detmcd.err, detmcd.std.err, ppmcd.std.err)
  max.y = max(fastmcd.err, detmcd.err, detmcd.std.err, ppmcd.std.err)

  max.y = min.y + 1.6*(max.y - min.y)

  par(mar = c(4, 4, 3, 3))
  y.lab = as.expression(y.labels[i.stat.type])
  # plot(x = NULL, y = NULL, xlim = c(min(AthruE), max(AthruE)), ylim = c(min.y, max.y),
  #      xlab = bquote("Scenario"), ylab = bquote(.(y.labels[i.stat.type])),
  #      main = bquote("(" * .(letters[i]) * ")" ~ "Setting" ~ .(cont.type)~"with"~epsilon~"="~.(conf[i, "eps"])))

  plot(AthruE, y = NULL, ylim = c(min.y, max.y),
       xlab = bquote("Scenario"), xaxt = "n", ylab = bquote(.(y.labels[i.stat.type])),
       main = bquote("(" * .(letters[i]) * ")" ~ .(cont.type)~"with"~epsilon~"="~.(eps)), type = "n")

  axis(side = 1, at = AthruE, labels = AthruE)

  lines(AthruE, ppmcd.std.err,  lty = 1, lwd = 1, col = "black")
  #lines(r.array, detmcd.err,     lty = 2, lwd = 1, col = "red")
  lines(AthruE, detmcd.std.err, lty = 2, lwd = , col = "red")
  lines(AthruE, fastmcd.err,    lty = 4, lwd = 1, col = "blue")

  points(AthruE, ppmcd.std.err,  lwd = 1, pch = 1, col = "black")
  #lines(r.array, detmcd.err,     lty = 2, lwd = 1, col = "red")
  points(AthruE, detmcd.std.err, lwd = 1, pch = 2, col = "red")
  points(AthruE, fastmcd.err,    lwd = 1, pch = 0, col = "blue")

  #legend("topleft", legend = c("PP MCD", "DetMCD (not sphered)", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
  #       lty = c(1, 2, 2, 4), lwd = c(2, 1, 2, 2), col = c("black", "red", "red", "blue"))

  legend("topleft", legend = c("PP MCD", "DetMCD (sphered)", paste("FastMCD(", mcd.nsamp, ")", sep = "")),
         lty = c(1, 2, 4), lwd = c(1, 1, 1), pch = c(1, 2, 0), col = c("black", "red", "blue"))
}

grDevices::dev.off()
