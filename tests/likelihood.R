# Deterministic checks of the published likelihoods and numerical integration.
library(metaZIP)

get_internal <- function(name) getFromNamespace(name, "metaZIP")
normal_rule <- get_internal(".normal_rule")
logprob <- get_internal(".metaZIP_logprob")
single_kernel <- get_internal(".single_log_kernel")
paired_kernel <- get_internal(".paired_log_kernel")
log_add <- get_internal(".log_add")
row_log_sum <- get_internal(".row_log_sum")
close_to <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(length(actual) == length(expected),
            all(is.finite(actual)), all(is.finite(expected)),
            max(abs(actual - expected)) <= tolerance)
}

rule <- normal_rule(61L)
close_to(sum(rule$weights), 1)
close_to(sum(rule$weights * rule$nodes), 0)
close_to(sum(rule$weights * rule$nodes^2), 1)
close_to(sum(rule$weights * rule$nodes^4), 3)
stopifnot(all(rule$weights > 0))
close_to(log_add(-1000, -1001), -1000 + log1p(exp(-1)))
stopifnot(identical(log_add(-Inf, -Inf), -Inf))
close_to(row_log_sum(matrix(c(-1000, -1001), nrow = 1L)),
         -1000 + log1p(exp(-1)))

alpha <- 0.8
beta <- 140
rate <- alpha / beta
mu <- 0.3
q <- 0.65
par <- c(log_alpha = log(alpha), log_rate = log(rate),
         mu = mu, sigma = 0, q = q)

# Independent closed-form Gamma integral, distinct from the NB/binomial engine.
gamma_probability <- function(y, exposure) {
  exp(alpha * log(beta) + lgamma(sum(y) + alpha) - lgamma(alpha) -
        (sum(y) + alpha) * log(beta + sum(exposure)) +
        sum(y * log(exposure) - lgamma(y + 1)))
}

paired <- list(y1 = c(0, 0, 3, 3), y0 = c(0, 2, 0, 2),
               n1 = c(250, 160, 400, 300), n0 = c(190, 210, 170, 450),
               design = "paired")
# Varying exposures and three effect values catch row/column recycling errors.
taus <- c(-0.4, 0.3, 0.8)
expected <- matrix(NA_real_, 4, 3)
for (j in seq_along(taus)) for (i in seq_along(paired$y1)) {
  y1 <- paired$y1[i]
  y0 <- paired$y0[i]
  e1 <- paired$n1[i] * exp(taus[j])
  e0 <- paired$n0[i]
  expected[i, j] <- q^2 * gamma_probability(c(y1, y0), c(e1, e0)) +
    q * (1 - q) * ((y0 == 0) * gamma_probability(y1, e1) +
                    (y1 == 0) * gamma_probability(y0, e0)) +
    (1 - q)^2 * (y1 == 0 && y0 == 0)
}
actual <- exp(paired_kernel(paired, alpha, rate, taus, q))
stopifnot(identical(dim(actual), c(4L, 3L)))
close_to(actual, expected)

# q=1 recovers the paired Gamma-Poisson likelihood, retaining dependence.
paired_plain <- exp(paired_kernel(paired, alpha, rate, mu, 1))
for (i in seq_along(paired$y1)) {
  close_to(paired_plain[i], gamma_probability(c(paired$y1[i], paired$y0[i]),
                                             c(paired$n1[i] * exp(mu), paired$n0[i])))
}
independent_zero <- gamma_probability(0, paired$n1[1] * exp(mu)) *
  gamma_probability(0, paired$n0[1])
stopifnot(abs(paired_plain[1] - independent_zero) > 0.01)

# Reference values from independent SciPy adaptive integration over a Normal z.
# Rows cover both-zero, control-positive, treatment-positive, and both-positive.
fixture <- paired
fixture$n1 <- rep(250, 4)
fixture$n0 <- rep(190, 4)
reference <- c(0.44442150923597984, 0.035439509607855094,
               0.029418730848995288, 0.006681908405915444)
par["sigma"] <- 0.5
close_to(exp(logprob(fixture, par, rule)), reference)
par["sigma"] <- 1.2
close_to(exp(logprob(fixture, par, rule)),
         c(0.44736087521253004, 0.03851143405457065,
           0.025188411063901087, 0.004512959822260233))
par["sigma"] <- 1e-4
close_to(exp(logprob(fixture, par, rule)),
         c(0.44343942722027213, 0.03441134785550084,
           0.030746941805525056, 0.007773227018429342))
par["sigma"] <- 0
close_to(logprob(paired, par, rule),
         as.vector(paired_kernel(paired, alpha, rate, mu, q)))

# Single-arm formulas and adaptive integration, including one-row input.
single <- list(y = c(0, 0, 2, 3), x = c(0, 1, 0, 1),
               ni = c(40, 250, 190, 300), design = "single")
expected_single <- sapply(taus, function(tau) {
  p <- mapply(function(y, n, x) gamma_probability(y, n * exp(tau * x)),
              single$y, single$ni, single$x)
  q * p + (1 - q) * (single$y == 0)
})
close_to(exp(single_kernel(single, alpha, rate, taus, q)), expected_single)
close_to(exp(single_kernel(single, alpha, rate, mu, 1)),
         stats::dnbinom(single$y, size = alpha,
                        mu = single$ni * rate * exp(mu * single$x)))
par["sigma"] <- 0.7
integrated <- vapply(seq_along(single$y), function(i) {
  stats::integrate(function(z) {
    count_prob <- stats::dnbinom(single$y[i], size = alpha,
      mu = single$ni[i] * rate * exp((mu + 0.7 * z) * single$x[i]))
    (q * count_prob + (1 - q) * (single$y[i] == 0)) * stats::dnorm(z)
  }, lower = -12, upper = 12, rel.tol = 1e-11, abs.tol = 1e-12)$value
}, numeric(1))
close_to(exp(logprob(single, par, rule)), integrated)
for (dat in list(single, paired)) {
  par["sigma"] <- 0.5
  batch <- logprob(dat, par, rule)
  for (i in seq_along(batch)) {
    one <- lapply(dat, function(z) if (is.numeric(z)) z[i] else z)
    close_to(logprob(one, par, rule), batch[i])
  }
}

# Fast observed-data check of nesting, dimensions and effect uncertainty.
eus <- utils::read.csv(system.file("extdata", "eus_fna.csv", package = "metaZIP"))
fit <- metaZIPs(eus$y, eus$x, eus$ni, method = c("POIF", "ZIPF"),
                control = list(warn = FALSE))
stopifnot(fit$out["ZIPF", "Log Likelihood"] >=
            fit$out["POIF", "Log Likelihood"] - 1e-8,
          identical(as.numeric(fit$out$npara), c(3, 4)),
          is.finite(fit$out["POIF", "SE"]), fit$out["POIF", "SE"] > 0)
# The exact sigma=0/q=1 limit must use the effect's information, not sigma as SE.
boundary <- fit$fits$POIF
boundary$model <- "ZIPR"
boundary$npara <- 5L
boundary <- get_internal(".fit_information")(boundary, fit$data, rule, fit$control)
close_to(boundary$se, fit$out["POIF", "SE"], tolerance = 1e-7)
stopifnot(grepl("sigma=0", boundary$boundary, fixed = TRUE),
          grepl("q=1", boundary$boundary, fixed = TRUE),
          grepl("held fixed", boundary$se_type, fixed = TRUE),
          identical(colnames(boundary$vcov), c("log_alpha", "log_rate", "mu")))
