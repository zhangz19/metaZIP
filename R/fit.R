.model_indices <- list(POIF = 1:3, POIR = 1:4, ZIPF = c(1:3, 5L), ZIPR = 1:5)

.metaZIP_control <- function(control) {
  defaults <- list(nquad = 61L, maxit = 1500L, reltol = 1e-10,
                   quad_tol = 1e-5, max_refine = 2L, sigma_max = 5,
                   alpha_min = 1e-6, alpha_max = 1e6, mu_bound = 20,
                   boundary_tol = 1e-5, warn = TRUE)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!names(control) %in% names(defaults)))))
    stop("control must be a named list of documented control settings.", call. = FALSE)
  ans <- utils::modifyList(defaults, control)
  for (nm in setdiff(names(ans), "warn")) {
    if (!is.numeric(ans[[nm]]) || length(ans[[nm]]) != 1L ||
        !is.finite(ans[[nm]])) stop("Invalid control setting: ", nm, call. = FALSE)
  }
  if (ans$nquad < 15 || ans$nquad > 301 || ans$nquad %% 1 != 0 ||
      ans$maxit < 1 || ans$maxit %% 1 != 0 ||
      ans$max_refine < 0 || ans$max_refine > 3 || ans$max_refine %% 1 != 0 ||
      ans$reltol <= 0 || ans$quad_tol <= 0 || ans$sigma_max <= 0 ||
      ans$alpha_min <= 0 || ans$alpha_max <= ans$alpha_min ||
      ans$mu_bound <= 0 || ans$boundary_tol <= 0 ||
      !is.logical(ans$warn) || length(ans$warn) != 1L || is.na(ans$warn))
    stop("Invalid range for a control setting; see ?metaZIPs.", call. = FALSE)
  ans
}

.metaZIP_validate <- function(values, counts, sizes) {
  if (!length(values) || any(vapply(values, function(z)
      !is.numeric(z) || !is.null(dim(z)) || !length(z) || any(!is.finite(z)), logical(1))))
    stop("Inputs must be nonempty, finite numeric vectors; handle missing values explicitly.",
         call. = FALSE)
  if (length(unique(lengths(values))) != 1L)
    stop("All input vectors must have the same length.", call. = FALSE)
  for (nm in counts) if (any(values[[nm]] < 0 | values[[nm]] %% 1 != 0))
    stop(nm, " must contain nonnegative integer event counts.", call. = FALSE)
  for (nm in sizes) if (any(values[[nm]] <= 0 | values[[nm]] %% 1 != 0))
    stop(nm, " must contain positive integer sample sizes.", call. = FALSE)
}

.metaZIP_start <- function(dat, para0, ctl) {
  if (dat$design == "single") {
    r0 <- sum(dat$y[dat$x == 0]) / sum(dat$ni[dat$x == 0])
    r1 <- sum(dat$y[dat$x == 1]) / sum(dat$ni[dat$x == 1])
  } else {
    r0 <- sum(dat$y0) / sum(dat$n0)
    r1 <- sum(dat$y1) / sum(dat$n1)
  }
  if (r0 == 0 || r1 == 0)
    stop("A finite log-risk-ratio fit requires at least one event in each group overall. ",
         "Individual zero-event studies are supported.", call. = FALSE)
  p <- c(log_alpha = 0, log_rate = log(r0), mu = log(r1 / r0), sigma = 0, q = 1)
  if (!is.null(para0)) {
    if (is.list(para0)) para0 <- unlist(para0, use.names = TRUE)
    if (!is.numeric(para0) || is.null(names(para0)) || any(!is.finite(para0)) ||
        anyDuplicated(names(para0)))
      stop("para0 must be a finite named numeric vector (or named list).", call. = FALSE)
    permitted <- c("alpha", "alph", "beta", "bet", "tau", "mu", "tau(mu)",
                   "sigma", "sig", "q", "c", "d")
    if (any(!names(para0) %in% permitted))
      stop("Unrecognized name in para0.", call. = FALSE)
    get_start <- function(keys, fallback) {
      k <- keys[keys %in% names(para0)]
      if (length(k)) unname(para0[[k[[1L]]]]) else fallback
    }
    alpha <- get_start(c("alpha", "alph"), exp(p[["log_alpha"]]))
    beta <- get_start(c("beta", "bet"), alpha / exp(p[["log_rate"]]))
    if (alpha <= 0 || beta <= 0) stop("Starting alpha and beta must be positive.", call. = FALSE)
    p[1:2] <- c(log(alpha), log(alpha / beta))
    p[["mu"]] <- get_start(c("mu", "tau", "tau(mu)"), p[["mu"]])
    p[["sigma"]] <- get_start(c("sigma", "sig"), 0)
    if (any(c("c", "d") %in% names(para0))) {
      if (!all(c("c", "d") %in% names(para0)) || any(para0[c("c", "d")] <= 0))
        stop("Legacy starts require both c and d to be positive.", call. = FALSE)
      p[["q"]] <- para0[["c"]] / (para0[["c"]] + para0[["d"]])
    }
    p[["q"]] <- get_start("q", p[["q"]])
    if (p[["sigma"]] < 0 || p[["sigma"]] > ctl$sigma_max ||
        p[["q"]] <= 0 || p[["q"]] > 1)
      stop("Starting sigma or q is outside its allowed range.", call. = FALSE)
  }
  p
}

.fit_one_model <- function(model, starts, nested, dat, rule, ctl) {
  idx <- .model_indices[[model]]
  lower <- c(log(ctl$alpha_min), log(1e-12), -ctl$mu_bound, 0, 1e-8)
  upper <- c(log(ctl$alpha_max), log(1e3), ctl$mu_bound, ctl$sigma_max, 1)
  base <- c(log_alpha = 0, log_rate = 0, mu = 0, sigma = 0, q = 1)
  objective <- function(x) {
    p <- base
    p[idx] <- x
    val <- -sum(.metaZIP_logprob(dat, p, rule))
    if (is.finite(val)) val else 1e100
  }
  candidates <- list()
  for (p in starts) {
    x <- pmin(upper[idx], pmax(lower[idx], p[idx]))
    names(x) <- names(base)[idx]
    fit <- tryCatch(stats::nlminb(x, objective, lower = lower[idx], upper = upper[idx],
        control = list(iter.max = ctl$maxit, eval.max = 3L * ctl$maxit,
                       rel.tol = ctl$reltol)), error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$objective) && fit$objective < 1e99) {
      full <- base
      full[idx] <- fit$par
      candidates[[length(candidates) + 1L]] <- list(par = full,
          objective = fit$objective, convergence = fit$convergence,
          message = fit$message, origin = "optimized")
    }
  }
  for (fit in nested) {
    fit$objective <- -sum(.metaZIP_logprob(dat, fit$par, rule))
    fit$origin <- "nested boundary"
    candidates[[length(candidates) + 1L]] <- fit
  }
  if (!length(candidates)) stop("Optimization failed for ", model, ".", call. = FALSE)
  vals <- vapply(candidates, `[[`, numeric(1), "objective")
  best <- candidates[[which.min(vals)]]
  # Prefer the exact nested boundary when objective differences are negligible.
  for (fit in nested) {
    val <- -sum(.metaZIP_logprob(dat, fit$par, rule))
    if (is.finite(val) && val <= best$objective + 1e-8) {
      best <- fit
      best$objective <- val
      best$origin <- "nested boundary"
    }
  }
  best$model <- model
  best$npara <- length(idx)
  best$attempts <- length(starts)
  best$lower <- lower
  best$upper <- upper
  best
}

.fit_model_set <- function(dat, requested, initial, rule, ctl, warm = NULL) {
  fits <- list()
  change <- function(p, ...) { z <- c(...); p[names(z)] <- z; p }
  run <- function(model, starts, nested = list()) {
    if (!is.null(warm[[model]])) starts <- c(list(warm[[model]]$par), starts)
    .fit_one_model(model, starts, nested, dat, rule, ctl)
  }
  p <- change(initial, sigma = 0, q = 1)
  fits$POIF <- run("POIF", list(p, change(p, log_alpha = log(.2)),
                               change(p, log_alpha = log(5))))
  if (any(c("POIR", "ZIPR") %in% requested)) {
    p <- fits$POIF$par
    starts <- lapply(c(.15, .6, 1.2), function(s) change(p, sigma = min(s, ctl$sigma_max)))
    if (initial[["sigma"]] > 0) starts <- c(list(change(initial, q = 1)), starts)
    fits$POIR <- run("POIR", starts, list(fits$POIF))
  }
  if (any(c("ZIPF", "ZIPR") %in% requested)) {
    p <- fits$POIF$par
    starts <- lapply(c(.35, .65, .9), function(q) change(p, q = q))
    if (initial[["q"]] < 1) starts <- c(list(change(initial, sigma = 0)), starts)
    fits$ZIPF <- run("ZIPF", starts, list(fits$POIF))
  }
  if ("ZIPR" %in% requested) {
    starts <- c(lapply(c(.4, .75, .95), function(q) change(fits$POIR$par, q = q)),
                lapply(c(.2, .7, 1.2), function(s)
                  change(fits$ZIPF$par, sigma = min(s, ctl$sigma_max))))
    if (initial[["sigma"]] > 0 || initial[["q"]] < 1) starts <- c(list(initial), starts)
    fits$ZIPR <- run("ZIPR", starts, list(fits$POIR, fits$ZIPF, fits$POIF))
  }
  fits
}

.fit_information <- function(fit, dat, rule, ctl) {
  p <- fit$par
  idx <- .model_indices[[fit$model]]
  boundary <- character()
  if (4L %in% idx && p[["sigma"]] <= ctl$boundary_tol) {
    idx <- setdiff(idx, 4L)
    boundary <- c(boundary, "sigma=0")
  }
  if (5L %in% idx && 1 - p[["q"]] <= ctl$boundary_tol) {
    idx <- setdiff(idx, 5L)
    boundary <- c(boundary, "q=1")
  }
  at_limit <- any(abs(p[idx] - fit$lower[idx]) < ctl$boundary_tol |
                  abs(p[idx] - fit$upper[idx]) < ctl$boundary_tol)
  if (at_limit) boundary <- c(boundary, "numerical limit")
  h <- matrix(NA_real_, length(idx), length(idx),
              dimnames = list(names(p)[idx], names(p)[idx]))
  v <- h
  se <- NA_real_
  objective <- function(x) {
    p[idx] <- x
    -sum(.metaZIP_logprob(dat, p, rule))
  }
  reliable <- fit$convergence == 0 && !at_limit && fit$quad_ok
  if (reliable) {
    h0 <- tryCatch(stats::optimHess(p[idx], objective,
                    control = list(ndeps = rep(1e-4, length(idx)))), error = function(e) NULL)
    if (!is.null(h0) && all(is.finite(h0))) {
      h <- (h0 + t(h0)) / 2
      dimnames(h) <- dimnames(v)
      ev <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
      reliable <- min(ev) > max(ev) * 1e-8 && min(ev) > 0
      if (reliable) {
        v <- solve(h)
        se <- sqrt(v["mu", "mu"])
        reliable <- is.finite(se) && se > 0
      }
    } else reliable <- FALSE
  }
  if (!reliable) { se <- NA_real_; v[,] <- NA_real_ }
  fit$se <- se
  fit$vcov <- v
  fit$hessian <- h
  fit$hessian_ok <- reliable
  fit$boundary <- paste(boundary, collapse = "; ")
  fit$se_type <- if (!reliable) "unavailable" else if (length(boundary))
    "Wald; boundary nuisance held fixed" else "Wald; observed information"
  fit
}

.metaZIP_fit <- function(dat, method, para0, conf.level, eps, control, call) {
  if (!is.character(method) || !length(method) || anyNA(method) ||
      any(!method %in% names(.model_indices)) || anyDuplicated(method))
    stop("method must contain unique values from POIF, POIR, ZIPF, ZIPR.", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be strictly between 0 and 1.", call. = FALSE)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0)
    stop("eps must be a positive finite number.", call. = FALSE)
  if (eps != .Machine$double.eps)
    warning("eps is retained for call compatibility and no longer sets parameter bounds; ",
            "use control for numerical settings.", call. = FALSE)
  ctl <- .metaZIP_control(control)
  initial <- .metaZIP_start(dat, para0, ctl)
  warm <- NULL
  for (refinement in 0:ctl$max_refine) {
    rule <- .normal_rule(ctl$nquad)
    fits <- .fit_model_set(dat, method, initial, rule, ctl, warm)
    check_rule <- .normal_rule(2L * ctl$nquad + 1L)
    differences <- vapply(fits, function(f)
      abs(sum(.metaZIP_logprob(dat, f$par, check_rule)) + f$objective), numeric(1))
    if (all(differences <= ctl$quad_tol) || refinement == ctl$max_refine) break
    warm <- fits
    ctl$nquad <- 2L * ctl$nquad + 1L
  }
  n <- if (dat$design == "single") length(dat$y) else length(dat$y1)
  for (model in names(fits)) {
    fits[[model]]$quad_difference <- differences[[model]]
    fits[[model]]$quad_ok <- differences[[model]] <= ctl$quad_tol
    fits[[model]] <- .fit_information(fits[[model]], dat, rule, ctl)
  }
  fits <- fits[method]
  z <- stats::qnorm((1 + conf.level) / 2)
  rows <- lapply(fits, function(f) {
    mu <- f$par[["mu"]]
    c(-f$objective, mu, f$se, exp(mu - z * f$se), exp(mu + z * f$se),
      f$npara, 2 * f$npara + 2 * f$objective,
      log(n) * f$npara + 2 * f$objective, exp(mu))
  })
  out <- as.data.frame(do.call(rbind, rows))
  names(out) <- c("Log Likelihood", "Estimate", "SE",
                  paste0(format(100 * conf.level, trim = TRUE), "% CI lower"),
                  "upper", "npara", "AIC", "BIC", "RR")
  parameters <- lapply(fits, function(f) {
    p <- f$par
    c(alpha = exp(p[["log_alpha"]]), beta = exp(p[["log_alpha"]] - p[["log_rate"]]),
      "tau(mu)" = p[["mu"]], sigma = p[["sigma"]], q = p[["q"]])
  })
  paraEst <- as.data.frame(do.call(rbind, parameters))
  diagnostics <- do.call(rbind, lapply(fits, function(f)
    data.frame(converged = f$convergence == 0, optimizer_code = f$convergence,
               boundary = f$boundary, hessian_ok = f$hessian_ok,
               se_type = f$se_type, quadrature_difference = f$quad_difference,
               quadrature_ok = f$quad_ok, nquad = ctl$nquad,
               message = f$message, stringsAsFactors = FALSE)))
  if (ctl$warn && any(!diagnostics$converged | !diagnostics$quadrature_ok | !diagnostics$hessian_ok))
    warning("Some fits have unavailable uncertainty estimates or numerical diagnostics; ",
            "inspect $diagnostics before interpreting results.", call. = FALSE)
  structure(list(out = out, paraEst = paraEst, diagnostics = diagnostics,
                 fits = fits, design = dat$design, nstudies = n, data = dat,
                 conf.level = conf.level, control = ctl, call = call,
                 estimand = "mean log risk ratio; RR is exp(mean log risk ratio)"),
            class = "metaZIP")
}

metaZIPs <- function(y, x, ni, method = c("POIF", "POIR", "ZIPF", "ZIPR"),
                     para0 = NULL, conf.level = .95, eps = .Machine$double.eps,
                     control = list()) {
  values <- list(y = y, x = x, ni = ni)
  .metaZIP_validate(values, "y", "ni")
  if (!all(x %in% c(0, 1)) || length(unique(x)) != 2L)
    stop("x must contain both groups, coded 0 and 1.", call. = FALSE)
  if (any(y > ni)) stop("Event counts cannot exceed sample sizes.", call. = FALSE)
  .metaZIP_fit(c(values, list(design = "single")), method, para0, conf.level,
               eps, control, match.call())
}

metaZIPr <- function(y1, n1, y0, n0, method = c("POIF", "POIR", "ZIPF", "ZIPR"),
                     para0 = NULL, conf.level = .95, eps = .Machine$double.eps,
                     control = list()) {
  values <- list(y1 = y1, n1 = n1, y0 = y0, n0 = n0)
  .metaZIP_validate(values, c("y1", "y0"), c("n1", "n0"))
  if (any(y1 > n1 | y0 > n0)) stop("Event counts cannot exceed sample sizes.", call. = FALSE)
  .metaZIP_fit(c(values, list(design = "paired")), method, para0, conf.level,
               eps, control, match.call())
}
