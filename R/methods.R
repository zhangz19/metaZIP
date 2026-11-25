print.metaZIP <- function(x, digits = 4, ...) {
  cat("metaZIP:", x$nstudies, if (x$design == "paired") "paired trials" else "single-arm studies", "\n")
  cat("Estimate = mean log risk ratio; RR = exp(Estimate).\n")
  print(x$out, digits = digits, ...)
  boundary <- rownames(x$diagnostics)[nzchar(x$diagnostics$boundary)]
  if (length(boundary)) cat("Boundary fits:", paste(boundary, collapse = ", "),
                            "-- inspect $diagnostics for SE interpretation.\n")
  invisible(x)
}

summary.metaZIP <- function(object, ...) {
  print(object, ...)
  cat("\nNumerical and uncertainty diagnostics:\n")
  print(object$diagnostics)
  invisible(object)
}

coef.metaZIP <- function(object, ...) {
  stats::setNames(object$out$Estimate, rownames(object$out))
}

confint.metaZIP <- function(object, parm, level = object$conf.level,
                            exponentiate = TRUE, ...) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be strictly between 0 and 1.", call. = FALSE)
  selected <- rownames(object$out)
  if (!missing(parm)) {
    if (is.numeric(parm)) {
      if (anyNA(parm) || any(parm %% 1 != 0 | parm < 1 | parm > length(selected)))
        stop("Invalid model index in parm.", call. = FALSE)
      selected <- selected[parm]
    } else if (is.character(parm) && !anyNA(parm) && all(parm %in% selected)) {
      selected <- parm
    } else stop("parm must identify fitted models.", call. = FALSE)
  }
  tab <- object$out[selected, , drop = FALSE]
  z <- stats::qnorm((1 + level) / 2)
  ans <- cbind(tab$Estimate - z * tab$SE, tab$Estimate + z * tab$SE)
  if (exponentiate) ans <- exp(ans)
  dimnames(ans) <- list(selected, paste0(format(100 * c((1 - level) / 2, (1 + level) / 2),
                                                trim = TRUE), "%"))
  ans
}

.select_fit <- function(object, model) {
  if (is.null(model)) {
    if (length(object$fits) != 1L)
      stop("Specify model when the object contains multiple fitted models.", call. = FALSE)
    model <- names(object$fits)[[1L]]
  }
  if (!is.character(model) || length(model) != 1L || is.na(model) || !model %in% names(object$fits))
    stop("model must name one fitted model.", call. = FALSE)
  object$fits[[model]]
}

logLik.metaZIP <- function(object, model = NULL, ...) {
  f <- .select_fit(object, model)
  structure(-f$objective, df = f$npara, nobs = object$nstudies, class = "logLik")
}

vcov.metaZIP <- function(object, model = NULL, ...) {
  .select_fit(object, model)$vcov
}

print.metaZIP_legacy <- function(x, digits = 4, ...) {
  cat("metaZIP historical EUS calculation:", x$nstudies, "single-arm studies\n")
  cat("Estimate = mean log risk ratio; RR = exp(Estimate).\n")
  cat("Legacy spread: inverse-Hessian SE for POIF/ZIPF; heterogeneity SD for POIR/ZIPR.\n")
  cat("Random-model legacy limits are not confidence intervals for the mean effect.\n")
  print(x$out, digits = digits, ...)
  cat("Legacy npara/AIC/BIC retain the original c,d parameter count.\n")
  invisible(x)
}

summary.metaZIP_legacy <- function(object, ...) {
  print(object, ...)
  cat("\nHistorical optimizer diagnostics (no adaptive quadrature check):\n")
  print(object$diagnostics)
  invisible(object)
}

coef.metaZIP_legacy <- function(object, ...) {
  stats::setNames(object$out$Estimate, rownames(object$out))
}

confint.metaZIP_legacy <- function(object, ...) {
  stop("Historical spreads include heterogeneity SDs, not SEs of the mean. ",
       "Use engine = 'current' for mean-effect confidence intervals; ",
       "historical limits are in $out.", call. = FALSE)
}

vcov.metaZIP_legacy <- function(object, ...) {
  stop("The historical engine does not provide a validated covariance matrix. ",
       "Use engine = 'current' for observed-information covariance estimates.",
       call. = FALSE)
}

logLik.metaZIP_legacy <- function(object, ...) {
  stop("Historical likelihoods and parameter counts are in $out. ",
       "Use engine = 'current' for a logLik object with identifiable parameter counts.",
       call. = FALSE)
}
