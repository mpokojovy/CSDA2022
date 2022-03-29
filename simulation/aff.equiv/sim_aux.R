#require("PPcovMcd")

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
sampleMeasures <- function(mcd.type = "PP", nrep = 50, p, n, eps = 0.5, rho.mult = 0.75, r = NA, delta = NA, psi = NA, standardize = TRUE) {
  m = min(floor(n*eps), n - floor((n + p + 1)/2))

  samp = drawSample(n, p, eps, rho.mult, r, delta, psi)
  x    = samp$dat
  
  #wgtFUN <- function(dist) (dist < qchisq(0.975, p))
  wgtFUN <- function(dist) (dist < qchisq(0.9, p)*median(dist)/qchisq(0.5, p))
  
  # Compute MCD for the original sample
  if (standardize) {
    mah.std = mah.standard(x)

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
    mcd = PPcovMcd(x, nsamp = mcd.type, wgtFUN = wgtFUN, PP.standardization = "diagonal")
  }
  
  mu0    = mcd$center
  Sigma0 = mcd$cov
  
  svd.ans = svd(Sigma0)
  inv.root.Sigma0 = svd.ans$u %*% diag(1/sqrt(svd.ans$d)) %*% t(svd.ans$v)
  rm(svd.ans)
  
  ##
  d_mu    = rep(0.0, nrep)
  d_Sigma = rep(0.0, nrep)

  imax = 1000L
  
  i     = 0L
  nsing = 0L
  
  while (i < imax) {
    if (i >= nrep) break
    
    d     = runif(p, min = 0.0, max = 1.0)
    D     = matrix(0.0, nrow = p, ncol = p); diag(D)     = d
    D.inv = matrix(0.0, nrow = p, ncol = p); diag(D.inv) = 1.0/d
    
    U = pracma::randortho(p)
    
    A     = U %*% D
    A.inv = D.inv %*% t(U)
    
    y = x %*% t(A)
    
    if (rcond(cov(y)) < 1E-12) {
      nsing = nsing + 1L
      next
    } else 
      i = i + 1L
    
    if (standardize) {
      mah.std = mah.standard(y)
      
      mcd = PPcovMcd(mah.std$z.score, nsamp = mcd.type, wgtFUN = wgtFUN, PP.standardization = "none")
      
      mcd$crit = mcd$crit + 2.0*as.numeric(determinant(mah.std$root.cov)$modulus)
      mcd$center = as.vector(mah.std$loc + mah.std$root.cov %*% matrix(mcd$center, ncol = 1))
      mcd$raw.center = as.vector(mah.std$loc + mah.std$root.cov %*% matrix(mcd$raw.center, ncol = 1))
      mcd$cov = mah.std$root.cov %*% mcd$cov %*% mah.std$root.cov
      mcd$raw.cov = mah.std$root.cov %*% mcd$raw.cov %*% mah.std$root.cov
    } else {
      mcd = PPcovMcd(y, nsamp = mcd.type, wgtFUN = wgtFUN, PP.standardization = "diagonal")
    }
    
    mu_A    = A.inv %*% mcd$center
    Sigma_A = A.inv %*% mcd$cov %*% t(A.inv)
    
    d_mu[i]    = sqrt(sum((mu_A - mu0)^2))
    
    d_Sigma[i] = base::kappa(inv.root.Sigma0 %*% Sigma_A %*% inv.root.Sigma0, exact = TRUE)
  }
  
  ans        = c(mean(d_mu), mean(d_Sigma), nsing)
  names(ans) = c("d_mu", "d_Sigma", "nsing")
  return(ans)
}

## Do simulation
doSimulation <- function(mcd.type = "PP", nsim = 100, nrep = 50, p, n, eps = 0.5, rho.mult = 0.75, r = NA, delta = NA, psi = NA, standardize = TRUE) {
  if (parallel.flag) {
    sim.resp = foreach(ind = 1:nsim, .inorder = FALSE, .packages = c("robustbase", "PPcovMcd"),
                       .export = c("drawSample", "getSummary", "sampleMeasures"), .combine = "rbind") %dopar%
                       sampleMeasures(mcd.type, nrep, p, n, eps, rho.mult, r, delta, psi, standardize)

    row.names(sim.resp) = NULL
  } else {
    sim.resp = replicate(nsim, sampleMeasures(mcd.type, nrep, p, n, eps, rho.mult, r, delta, psi, standardize), simplify = "array")
  }

  nam = c("nsim", "n", "p", "eps", "rho.mult", "rho.mult.sq", "r", "delta", "psi",
          "d_mu.avg", "d_mu.std", "d_Sigma.avg", "d_Sigma.std", "nsing")

  if (parallel.flag) {
    # for parallel
    ans = c(nsim, n, p, eps, rho.mult, rho.mult^2, r, delta, psi,
            mean(sim.resp[, "d_mu"]), sd(sim.resp[, "d_mu"]), mean(sim.resp[, "d_Sigma"]), sd(sim.resp[, "d_Sigma"]), sum(sim.resp[, "nsing"]))
  } else {
    # for non-parallel
    ans = c(nsim, n, p, eps, rho.mult, rho.mult^2, r, delta, psi,
            mean(sim.resp["d_mu", ]), sd(sim.resp["d_mu", ]), mean(sim.resp["d_Sigma", ]), sd(sim.resp["d_Sigma", ]), sum(sim.resp["nsing", ]))
  }

  names(ans) = nam

  return(ans)
}