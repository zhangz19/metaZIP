# Reproduce the original EUS-FNA calculation after installing metaZIP:
# source(system.file("examples", "eus_fna.R", package = "metaZIP"))
# This example fits observed data and writes no files.

library(metaZIP)
data("eus_fna", package = "metaZIP")

fit_eus <- with(eus_fna, metaZIPs(y = y, x = x, ni = ni, engine = "legacy"))
print(fit_eus, digits = 6)
print(fit_eus$paraEst, digits = 6)
print(fit_eus$diagnostics)

# The historical random-model spread is sigma, not SE(mu).
# For observed-information SEs and confidence intervals for the mean effect:
# fit_current <- with(eus_fna, metaZIPs(y, x, ni))
# fit_current$out
# confint(fit_current)
