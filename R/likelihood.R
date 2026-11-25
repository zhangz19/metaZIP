# The four likelihoods retain the Gamma baseline models in the original util.R.
# q is the probability of the count component, not the structural-zero component.

.log_add <- function(a, b) {
  m <- pmax(a, b)
  ans <- m + log1p(exp(-abs(a - b)))
  ans[!is.finite(m)] <- m[!is.finite(m)]
  ans
}

.row_log_sum <- function(a) {
  m <- apply(a, 1L, max)
  ans <- m + log(rowSums(exp(a - m)))
  ans[!is.finite(m)] <- m[!is.finite(m)]
  ans
}

# Golub-Welsch rule for a STANDARD NORMAL expectation: sum(w * f(z)).
.normal_rule <- function(n) {
  jacobi <- matrix(0, n, n)
  k <- seq_len(n - 1L)
  jacobi[cbind(k, k + 1L)] <- sqrt(k)
  jacobi[cbind(k + 1L, k)] <- sqrt(k)
  ee <- eigen(jacobi, symmetric = TRUE)
  w <- ee$vectors[1L, ]^2
  keep <- w > 0
  list(nodes = ee$values[keep], weights = w[keep] / sum(w[keep]))
}

.single_log_kernel <- function(dat, alpha, rate, tau, q) {
  nr <- length(dat$y)
  nk <- length(tau)
  mean_count <- exp(log(dat$ni) + log(rate) + outer(dat$x, tau))
  ans <- matrix(stats::dnbinom(dat$y, size = alpha,
                               mu = mean_count, log = TRUE), nr, nk)
  if (q == 1) return(ans)
  ans <- ans + log(q)
  iz <- dat$y == 0
  if (any(iz)) ans[iz, ] <- .log_add(ans[iz, , drop = FALSE], log1p(-q))
  ans
}

.paired_log_kernel <- function(dat, alpha, rate, tau, q) {
  nr <- length(dat$y1)
  nk <- length(tau)
  log_e1 <- outer(log(dat$n1), tau, "+")
  log_e0 <- matrix(log(dat$n0), nr, nk)
  log_et <- .log_add(log_e1, log_e0)
  total <- dat$y1 + dat$y0
  # Shared-Gamma joint pmf = NB(total) * Binomial(y1 | total).
  joint <- matrix(stats::dnbinom(total, size = alpha,
                                 mu = exp(log(rate) + log_et), log = TRUE), nr, nk)
  joint <- joint + lchoose(total, dat$y1) +
    dat$y1 * (log_e1 - log_et) + dat$y0 * (log_e0 - log_et)
  if (q == 1) return(joint)
  ans <- 2 * log(q) + joint
  z0 <- dat$y0 == 0
  z1 <- dat$y1 == 0
  mix_weight <- log(q) + log1p(-q)
  if (any(z0)) {
    g1 <- matrix(stats::dnbinom(dat$y1, size = alpha,
                                mu = exp(log(rate) + log_e1), log = TRUE), nr, nk)
    ans[z0, ] <- .log_add(ans[z0, , drop = FALSE],
                          mix_weight + g1[z0, , drop = FALSE])
  }
  if (any(z1)) {
    g0 <- matrix(stats::dnbinom(dat$y0, size = alpha,
                                mu = dat$n0 * rate, log = TRUE), nr, nk)
    ans[z1, ] <- .log_add(ans[z1, , drop = FALSE],
                          mix_weight + g0[z1, , drop = FALSE])
  }
  zz <- z0 & z1
  if (any(zz)) ans[zz, ] <- .log_add(ans[zz, , drop = FALSE], 2 * log1p(-q))
  ans
}

.metaZIP_logprob <- function(dat, par, rule) {
  alpha <- exp(par[["log_alpha"]])
  rate <- exp(par[["log_rate"]])
  sigma <- par[["sigma"]]
  tau <- if (sigma == 0) par[["mu"]] else par[["mu"]] + sigma * rule$nodes
  kernel <- if (dat$design == "single") .single_log_kernel else .paired_log_kernel
  ll <- kernel(dat, alpha, rate, tau, par[["q"]])
  if (sigma == 0) return(as.vector(ll))
  .row_log_sum(sweep(ll, 2L, log(rule$weights), "+"))
}
