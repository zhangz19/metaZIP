# metaZIP

**Meta-analysis of rare events with Poisson mixture models.** Fit independent
single-arm studies or paired two-arm trials using Gamma baseline rates and
fixed or normally distributed log risk ratios. Individual zero-event studies
remain in the analysis.

## Install and reproduce the EUS-FNA example

The package includes the 51-study EUS-FNA dataset. This call runs the original
EUS calculation with its original starting values, optimizers, and quadrature
settings:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("zhangz19/metaZIP")

library(metaZIP)
data(eus_fna)
fit_eus <- with(eus_fna, metaZIPs(y, x, ni, engine = "legacy"))
fit_eus$out
fit_eus$paraEst
fit_eus$diagnostics
```

R 4.1 or later is required. Installation also installs the `gaussquad` and
`statmod` dependencies used by the legacy calculation. No external data file
or separately sourced `util.R` is needed. An installed example is available as:

```r
source(system.file("examples", "eus_fna.R", package = "metaZIP"))
```

### Expected historical EUS results

| Model | Log likelihood | Estimated log RR | Legacy spread | Legacy limits |
| --- | ---: | ---: | ---: | --- |
| POIF | -67.886195 | 0.365392 | 0.775646 | (0.3151, 6.5905) |
| POIR | -67.836258 | 0.127773 | 0.820970 | (0.2273, 5.6794) |
| ZIPF | -66.237843 | 0.470112 | approximately 0.7722 | approximately (0.3523, 7.269) |
| ZIPR | -66.179162 | 0.260051 | 0.715695 | (0.3190, 5.2740) |

These reproduce the historical EUS console output closely. Small numerical
differences can occur across R versions and platforms, particularly for the
ZIPF inverse-Hessian calculation.

`Legacy spread` follows the original function: it is an inverse-Hessian SE
for POIF and ZIPF, and the fitted heterogeneity SD (`sigma`) for POIR and ZIPR.
The `Spread type` column identifies the quantity in each row. `Legacy lower`
and `Legacy upper` reproduce `exp(Estimate +/- z * Legacy spread)`, with
`z = qnorm((1 + conf.level) / 2)`. For the random-effect models, these are
**not confidence limits for the estimated mean effect**. The legacy result
has its own class, and `confint()` and `vcov()` do not interpret it as such.

The `Legacy npara`, `Legacy AIC`, and `Legacy BIC` columns retain the original
counts of 3, 4, 5, and 6 parameters. The two ZIP models count both `c` and `d`,
although their likelihoods identify only `q = c/(c + d)`. These historical
criteria are included for comparison with the original output; `logLik()`
does not attach those redundant counts as model degrees of freedom.

The log likelihoods and effect estimates closely match the original
manuscript's Table 3. Its ZIPR entry of SE 0.64 and interval (0.37, 4.55)
is not reproduced by the supplied original function; the reproducible
historical values are shown above.

## Mean-effect standard errors and confidence intervals

The default `engine = "current"` uses the same four models with identifiable
mixture parameters, checked quadrature, and observed-information uncertainty
for the estimated mean log risk ratio:

```r
fit_current <- with(eus_fna, metaZIPs(y, x, ni))
fit_current$out           # SE and confidence limits for the estimated mean
fit_current$paraEst       # alpha, beta, tau(mu), sigma and q
fit_current$diagnostics
coef(fit_current)         # mean log risk ratios by model
confint(fit_current)      # exponentiated Wald confidence intervals
```

In this output, `sigma` describes between-study variation in log risk ratios;
`SE` measures sampling uncertainty in the estimated mean. When `sigma = 0`
or `q = 1`, the boundary nuisance parameter is held fixed for that
information calculation, as identified in `diagnostics$se_type`. Unreliable
information or numerical checks produce unavailable uncertainty estimates.

## Models and data

| Model | Gamma baseline rates | Log risk ratio | Extra zero component | Identifiable parameters |
| --- | --- | --- | --- | --- |
| POIF | Vary across studies | Common | No | 3 |
| POIR | Vary across studies | Normal random effect | No | 4 |
| ZIPF | Vary across studies | Common | Yes | 4 |
| ZIPR | Vary across studies | Normal random effect | Yes | 5 |

`Estimate` and `coef()` report the mean log risk ratio, **mu**. `RR` is
**exp(mu)**; for random-effect models this is the geometric mean study-level
risk ratio. `q` is the probability of the Poisson count component and
**1 - q** is the structural-zero probability.

`metaZIPs()` treats each row as an independent single-arm study: `x = 1` is
the treatment group and `x = 0` is the reference. For EUS-FNA these are 19G
and 22/25G needles, respectively. Single-arm comparisons rely on
between-study comparability and do not themselves adjust for confounding.
Data provenance and column definitions are in `?eus_fna`.

Select individual models with `method = "ZIPR"` or a character vector of
model names. Starting values and numerical controls are documented in
`?metaZIPs`. The current engine accepts controls such as
`control = list(nquad = 81)`. The legacy engine fixes its original numerical
settings and accepts only `control = list(warn = FALSE)` to suppress
diagnostic warnings; it does not silently apply current-engine controls.

### Paired-trial interface

The existing `metaZIPr()` interface and bundled 48-trial rosiglitazone MI
dataset remain available:

```r
data(rosiglitazone_mi)
fit_mi <- with(rosiglitazone_mi, metaZIPr(y1, n1, y0, n0))
fit_mi$out
fit_mi$diagnostics
```

`metaZIPr()` preserves the shared Gamma baseline within each paired trial;
matching positions in its four input vectors must refer to the same trial.
The two arms have independent count-component indicators with the same `q`.
The data come from [Cai, Parast and Ryan (2010)](https://doi.org/10.1002/sim.3964);
see `?rosiglitazone_mi`. The legacy engine is specific to `metaZIPs()`.

## Scope and citation

Version 0.1.1 adds the historical EUS calculation alongside the existing
current engine. The examples fit observed data; they do not reproduce the
manuscript's simulation study. The historical exploratory script remains
in `tests/codes.R` in the repository and is excluded from package builds.

Run `citation("metaZIP")` for the software citation. Software is provided
under the MIT license, copyright 2022 Zhen Zhang; see [LICENSE.md](LICENSE.md).
Dataset documentation identifies the original sources of the aggregate
study counts.
