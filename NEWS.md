# metaZIP 0.1.1

* Add `engine = "legacy"` to `metaZIPs()` to reproduce the historical EUS-FNA
  calculation with its original starts, optimizers, and quadrature settings.
* Label historical spreads, interval limits, parameter counts, and information
  criteria explicitly. Random-model spreads are heterogeneity SDs, not SEs
  of the estimated mean effect.
* Return a separate `metaZIP_legacy` class with print, summary, and coefficient
  methods; do not expose historical spreads as mean-effect confidence intervals.
* Add an EUS-only installed example and regression checks against the supplied
  original console results. Document the remaining difference from Table 3's
  ZIPR SE and interval.
* Retain the current engine as the default and preserve the paired-trial API.

# metaZIP 0.1.0

* Package the original `metaZIPs()` and `metaZIPr()` interfaces with the
  POIF, POIR, ZIPF and ZIPR Gamma-baseline models.
* Include the EUS-FNA and 48-trial rosiglitazone MI datasets and real examples.
* Replace nonidentifiable beta mixing parameters `c,d` with their identifiable
  ratio `q`, preserving the observed-data models and legacy starting values.
* Evaluate likelihoods on the log scale, use quadrature centered on the normal
  random effect, and retain exact nested-model boundary solutions.
* Calculate mean log-risk-ratio standard errors from observed information;
  identify intervals that hold boundary nuisance parameters fixed.
* Add input validation, convergence and quadrature diagnostics, information
  criteria, and standard display and extraction methods.
* Retain `eps` for call compatibility, with a warning for nondefault values.
* Exclude the historical exploratory simulation script from package builds.

