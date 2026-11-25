# Historical EUS results supplied with the original metaZIPs() console output.
# These checks reproduce that calculation, not every rounded Table 3 entry.
library(metaZIP)
data("eus_fna", package = "metaZIP")

close_to <- function(actual, expected, tolerance) {
  stopifnot(length(actual) == length(expected),
            all(is.finite(actual)),
            max(abs(actual - expected)) <= tolerance)
}
error_message <- function(expr) {
  result <- tryCatch(force(expr), error = function(e) e)
  stopifnot(inherits(result, "error"))
  conditionMessage(result)
}

fit <- with(eus_fna, metaZIPs(y, x, ni, engine = "legacy",
                            control = list(warn = FALSE)))
models <- c("POIF", "POIR", "ZIPF", "ZIPR")
stopifnot(inherits(fit, "metaZIP_legacy"), !inherits(fit, "metaZIP"),
          identical(rownames(fit$out), models),
          !"SE" %in% names(fit$out),
          identical(names(fit$out),
                    c("Log Likelihood", "Estimate", "Legacy spread",
                      "Legacy lower", "Legacy upper", "Legacy npara",
                      "Legacy AIC", "Legacy BIC", "RR", "Spread type")))

close_to(fit$out[["Log Likelihood"]],
         c(-67.886194884, -67.836258377, -66.237842719, -66.179161504), 1e-4)
close_to(coef(fit), c(.365391844, .127773406, .470111992, .260050641), 5e-4)
# The redundant c,d parameterization makes the historical ZIPF Hessian
# sensitive to numerical details; its supplied SE was 0.7721672.
close_to(fit$out[["Legacy spread"]],
         c(.7756458, .8209697, .7721672, .7156950), 1e-3)
close_to(fit$out[["Legacy lower"]],
         c(.3151061, .2273425, .3522874, .3189612), .02)
close_to(fit$out[["Legacy upper"]],
         c(6.590503, 5.679393, 7.268369, 5.273990), .02)
close_to(fit$out[["Legacy AIC"]],
         c(141.7724, 143.6725, 142.4757, 144.3583), 1e-3)
close_to(fit$out[["Legacy BIC"]],
         c(147.5679, 151.3998, 152.1348, 155.9493), 1e-3)
stopifnot(identical(as.numeric(fit$out[["Legacy npara"]]), c(3, 4, 5, 6)),
          identical(fit$out[["Spread type"]],
                    c("Inverse-Hessian SE", "Heterogeneity SD",
                      "Inverse-Hessian SE", "Heterogeneity SD")),
          all(fit$diagnostics$converged))

# An explicit historical label must accompany the SD, and standard inference
# extractors must not reinterpret it as uncertainty in the estimated mean.
close_to(fit$out[c("POIR", "ZIPR"), "Legacy spread"],
         fit$paraEst[c("POIR", "ZIPR"), "sigma"], 1e-12)
for (extractor in list(stats::confint, stats::vcov, stats::logLik)) {
  msg <- error_message(extractor(fit))
  stopifnot(grepl("legacy|historical|current", msg, ignore.case = TRUE))
}

# ZIPR-only calls still fit their required starting model internally and
# return just the requested row to the user.
zipr <- with(eus_fna, metaZIPs(y, x, ni, method = "ZIPR", engine = "legacy",
                             control = list(warn = FALSE)))
stopifnot(identical(rownames(zipr$out), "ZIPR"),
          identical(names(coef(zipr)), "ZIPR"))
close_to(zipr$out[["Log Likelihood"]], -66.179161504, 1e-4)
close_to(coef(zipr), .260050641, 5e-4)

# Unsupported modern controls and impossible starting values should fail
# explicitly instead of silently changing the historical calculation.
stopifnot(grepl("control|nquad", error_message(with(eus_fna,
  metaZIPs(y, x, ni, engine = "legacy", control = list(nquad = 61)))),
  ignore.case = TRUE))
invisible(error_message(with(eus_fna,
  metaZIPs(y, x, ni, engine = "legacy", para0 = c(alpha = -1),
           control = list(warn = FALSE)))))
invisible(error_message(with(eus_fna,
  metaZIPs(y, x, ni, engine = "legacy", para0 = c(q = 1),
           control = list(warn = FALSE)))))
