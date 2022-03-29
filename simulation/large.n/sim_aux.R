## Data generation function
drawSample <- function(n, p, eps = 0.5, rho.mult = 0.75, r = NA, delta = NA, psi = NA) {
  # Generate correlated "clean" data
  m = min(floor(n*eps), n - floor((n + p + 1)/2)) # number of "bad" pts

  y = matrix(rnorm(n*p), ncol = p)
  
  G <- function(rho) {
    G = matrix(rho, nrow = p, ncol = p)
    diag(G) = 1
    return(G)
  }

  R2 <- function(x, ind = 1) {
    Ind = setdiff(1:p, ind)
    Cor = cor(x)
    return(sum(solve(Cor[Ind, Ind], matrix(Cor[Ind, ind], ncol = 1))*Cor[Ind, ind]))
  }

  rho.ast = if (R2(y[1:(n - m), ]) >= rho.mult^2) 0.0 else uniroot(function(rho) R2(y[1:(n - m), ] %*% G(rho)) - rho.mult^2, interval = c(0.0, rho.mult))$root

  if (m > 0) {
    if (is.na(psi))
    {
      b = rep(1.0, p)
      z = rnorm(p)
      a0 = z - sum(z*b)*b/p
      a0 = a0/sqrt(sum(a0^2))

      y[(n - m + 1):n, ] = sweep(y[(n - m + 1):n, ]*delta, MARGIN = 2, STATS = (r*sqrt(p))*a0, FUN = "+")
    } else {
      ## Add radial contamination
      cutoff = qchisq(0.8, df = p)
      alpha0 = pchisq(cutoff/psi, df = p)
      Q <- function(alpha) psi*qchisq(alpha*(1 - alpha0) + alpha0, df = p) # Q = sol to F(r^2/psi) - F(r0^2/psi))/(1 - F(r0^2/psi)) = alpha

      runif.const = runif(1, min = 0.0, max = 1.0)
      
      for (i in (n - m + 1):n) {
        xi = y[i, ]
        y[i, ] = xi*sqrt(Q(runif.const)/sum(xi^2))
      }
    }}

  ## Prepare output
  ans = NULL
  cor = G(rho.ast)
  ans$dat = y %*% cor
  ans$rho = rho.ast
  ans$cor = cor

  return(ans)
}

## Get relevant statistics from MCD
getSummary <- function(mcd, G, m) {
  p = ncol(G)
  n = mcd$n.obs

  nam = c("obj", "err.mu", "err.Sigma1", "err.Sigma2", "breakdown", "impurity", "n.outlier", "bulk.size", "n.C.steps")

  G.inv = solve(G)

  obj     = exp(mcd$crit + 2.0*as.numeric(determinant(G.inv)$modulus))

  center  = G.inv %*% matrix(mcd$center, ncol = 1)
  scatter = G.inv %*% mcd$cov %*% G.inv

  eig  = eigen(scatter, symmetric = TRUE, only.values = TRUE)

  err.Sigma1 = max(eig$values)/min(eig$values)
  err.Sigma2 = max(abs(eig$values - 1.0), abs(1/eig$values - 1.0))

  breakdown = if (m == 0L) 0 else (sum(mcd$mcd.wt[(n - m + 1):n]) > 0)
  impurity  = if (m == 0L) 0 else sum(mcd$mcd.wt[(n - m + 1):n])/max(1, sum(mcd$mcd.wt)) # outliers are assumed to come at the end of the sample
  n.outlier = if (m == 0L) 0 else sum(mcd$mcd.wt[(n - m + 1):n])
  bulk.size = sum(mcd$mcd.wt)

  ans = c(obj,
          sum(center^2),
          log10(err.Sigma1), log10(err.Sigma2),
          breakdown, impurity, n.outlier, bulk.size,
          sum(mcd$n.csteps))

  names(ans) = nam

  return(ans)
}

## Apply MCD
applyMCD <- function(mcd.type = "PP", p, n, eps = 0.5, rho.mult = 0.75, r = NA, delta = NA, psi = NA, standardize = TRUE) {
  m = min(floor(n*eps), n - floor((n + p + 1)/2))

  samp = drawSample(n, p, eps, rho.mult, r, delta, psi)

  #wgtFUN <- function(dist) (dist < qchisq(0.975, p))
  wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))

  ptm <- proc.time()
  if (standardize) {
    mah.std = mah.standard(samp$dat)

    if (mcd.type == "PP")
      # no need to standardize again
      mcd = PPcovMcd(mah.std$z.score, nsamp = mcd.type, wgtFUN = wgtFUN, PP.standardization = "none")
    else
      # not standardized in the first place
      mcd = PPcovMcd(mah.std$z.score, nsamp = mcd.type, wgtFUN = wgtFUN)

    mcd$crit = mcd$crit + 2.0*as.numeric(determinant(mah.std$root.cov)$modulus)
    mcd$center = as.vector(mah.std$loc + mah.std$root.cov %*% matrix(mcd$center, ncol = 1))
    mcd$raw.center = as.vector(mah.std$loc + mah.std$root.cov %*% matrix(mcd$raw.center, ncol = 1))
    mcd$cov = mah.std$root.cov %*% mcd$cov %*% mah.std$root.cov
    mcd$raw.cov = mah.std$root.cov %*% mcd$raw.cov %*% mah.std$root.cov
  } else {
    mcd = PPcovMcd(samp$dat, nsamp = mcd.type, wgtFUN = wgtFUN)
  }
  time = unname((proc.time() - ptm)[3])

  mcd.sum = getSummary(mcd, samp$cor, m);

  nam = c(names(mcd.sum), "time")
  ans = c(mcd.sum, time)
  names(ans) = nam

  return(ans)
}

## Do simulation
doSimulation <- function(mcd.type = "PP", nsim = 100, p, n, eps = 0.5, rho.mult = 0.75, r = NA, delta = NA, psi = NA, standardize = TRUE) {
  if (parallel.flag) {
    sim.resp = foreach(ind = 1:nsim, .inorder = FALSE, .packages = c("robustbase", "PPcovMcd"),
                       .export = c("drawSample", "getSummary", "applyMCD"), .combine = "rbind") %dopar% 
                       applyMCD(mcd.type, p, n, eps, rho.mult, r, delta, psi, standardize)
    
    row.names(sim.resp) = NULL
  } else {
    sim.resp = replicate(nsim, applyMCD(mcd.type, p, n, eps, rho.mult, r, delta, psi, standardize), simplify = "array")
  }

  nam = c("nsim", "n", "p", "eps", "rho.mult", "rho.mult.sq", "r", "delta", "psi",
          "obj.avg", "obj.std", "err.mu.avg", "err.mu.std",
          "err.Sigma1.avg", "err.Sigma1.std", "err.Sigma2.avg", "err.Sigma2.std", #"err.Sigma3.avg", "err.Sigma3.std", "err.Sigma4.avg", "err.Sigma4.std",
          "breakdown.chance", "impurity.avg", "impurity.std", "n.outlier.avg", "n.outlier.std", "bulk.size.avg", "bulk.size.std",
          "n.C.steps.avg", "n.C.steps.std", "time.avg", "time.std")

  if (parallel.flag) {
    # for parallel
    ans = c(nsim, n, p, eps, rho.mult, rho.mult^2, r, delta, psi,
            mean(sim.resp[, "obj"]), sd(sim.resp[, "obj"]), mean(sim.resp[, "err.mu"]), sd(sim.resp[, "err.mu"]),
            mean(sim.resp[, "err.Sigma1"]), sd(sim.resp[, "err.Sigma1"]), mean(sim.resp[, "err.Sigma2"]), sd(sim.resp[, "err.Sigma2"]),
            mean(sim.resp[, "breakdown"]), mean(sim.resp[, "impurity"]), sd(sim.resp[, "impurity"]),
            mean(sim.resp[, "n.outlier"]), sd(sim.resp[, "n.outlier"]), mean(sim.resp[, "bulk.size"]), sd(sim.resp[, "bulk.size"]),
            mean(sim.resp[, "n.C.steps"]), sd(sim.resp[, "n.C.steps"]), mean(sim.resp[, "time"]), sd(sim.resp[, "time"]))
  } else {
    # for non-parallel
    ans = c(nsim, n, p, eps, rho.mult, rho.mult^2, r, delta, psi,
            mean(sim.resp["obj", ]), sd(sim.resp["obj", ]), mean(sim.resp["err.mu", ]), sd(sim.resp["err.mu", ]),
            mean(sim.resp["err.Sigma1", ]), sd(sim.resp["err.Sigma1", ]), mean(sim.resp["err.Sigma2", ]), sd(sim.resp["err.Sigma2", ]),
            #mean(sim.resp["err.Sigma3", ]), sd(sim.resp["err.Sigma3", ]), mean(sim.resp["err.Sigma4", ]), sd(sim.resp["err.Sigma4", ]),
            mean(sim.resp["breakdown", ]), mean(sim.resp["impurity", ]), sd(sim.resp["impurity", ]),
            mean(sim.resp["n.outlier", ]), sd(sim.resp["n.outlier", ]), mean(sim.resp["bulk.size", ]), sd(sim.resp["bulk.size", ]),
            mean(sim.resp["n.C.steps", ]), sd(sim.resp["n.C.steps", ]), mean(sim.resp["time", ]), sd(sim.resp["time", ]))
  }

  # ans = c(nsim, n, p, eps, rho.mult, rho.mult^2, r, delta, psi, runif(length(nam) - 9))
  names(ans) = nam

  return(ans)
}