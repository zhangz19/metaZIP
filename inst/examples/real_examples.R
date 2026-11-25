# Run after installing metaZIP:
# source(system.file("examples", "real_examples.R", package = "metaZIP"))
# Or: Rscript real_examples.R path/to/output_directory
# This script fits the two observed datasets; it does not run a simulation study.

library(metaZIP)
data("eus_fna", package = "metaZIP")
data("rosiglitazone_mi", package = "metaZIP")

fit_eus <- with(eus_fna, metaZIPs(y = y, x = x, ni = ni))
fit_ros <- with(rosiglitazone_mi,
                metaZIPr(y1 = y1, n1 = n1, y0 = y0, n0 = n0))
print(fit_eus)
print(fit_ros)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1L]] else "metaZIP-results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_tables <- function(fit, prefix) {
  for (component in c("out", "paraEst", "diagnostics")) {
    tab <- fit[[component]]
    tab <- cbind(Model = rownames(tab), tab)
    rownames(tab) <- NULL
    write.csv(tab, file.path(output_dir, paste0(prefix, "_", component, ".csv")),
              row.names = FALSE)
  }
}
save_tables(fit_eus, "EUS_FNA")
save_tables(fit_ros, "Rosiglitazone_MI_48")
saveRDS(list(eus_fna = fit_eus, rosiglitazone_mi = fit_ros),
         file.path(output_dir, "real_example_fits.rds"), version = 2)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("\nSaved results to", normalizePath(output_dir), "\n")
