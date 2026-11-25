library(metaZIP)
data("eus_fna", package = "metaZIP")
data("rosiglitazone_mi", package = "metaZIP")

# These are observed-data regressions and likelihood checks, not a simulation.
stopifnot(nrow(eus_fna) == 51L, sum(eus_fna$y) == 77L,
          sum(eus_fna$ni) == 6172L, nrow(rosiglitazone_mi) == 48L,
          sum(rosiglitazone_mi$y1) == 86L, sum(rosiglitazone_mi$y0) == 72L)
fits <- list(with(eus_fna, metaZIPs(y, x, ni)),
             with(rosiglitazone_mi, metaZIPr(y1, n1, y0, n0)))
for (fit in fits) {
  ll <- setNames(fit$out[["Log Likelihood"]], rownames(fit$out))
  stopifnot(ll[["POIR"]] >= ll[["POIF"]] - 1e-7,
            ll[["ZIPF"]] >= ll[["POIF"]] - 1e-7,
            ll[["ZIPR"]] >= max(ll[c("POIR", "ZIPF")]) - 1e-7,
            identical(as.numeric(fit$out$npara), c(3, 4, 4, 5)),
            all(fit$diagnostics$converged), all(fit$diagnostics$quadrature_ok),
            all(is.finite(fit$out$SE)))
}
# Independently reported Table 3 likelihoods/point estimates survive cleanup.
stopifnot(max(abs(fits[[1]]$out[["Log Likelihood"]] -
                    c(-67.886, -67.836, -66.237, -66.179))) < .001,
          max(abs(coef(fits[[1]]) - c(.36, .13, .47, .26))) < .01)
# Heterogeneity SD must not be substituted for the SE of the pooled mean.
stopifnot(all(abs(fits[[1]]$out[c("POIR", "ZIPR"), "SE"] -
                   fits[[1]]$paraEst[c("POIR", "ZIPR"), "sigma"]) > .1))
# The 48-trial dataset selects the no-inflation/no-random-slope boundary.
stopifnot(all(abs(fits[[2]]$paraEst$q - 1) < 1e-6),
          all(fits[[2]]$paraEst$sigma < 1e-5),
          all(nzchar(fits[[2]]$diagnostics$boundary[-1L])))
# Common data failures must be explicit rather than silently dropping studies.
fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
stopifnot(fails(metaZIPs(c(0, 1), c(0, 1), c(10, 10))),
          fails(metaZIPs(c(1, NA), c(0, 1), c(10, 10))),
          fails(metaZIPr(1, 10, 1:2, 10:11)))
