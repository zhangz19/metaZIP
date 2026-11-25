# metaZIP

**Meta-analysis of rare events with Poisson mixture models.** Fit independent
single-arm studies or paired two-arm trials using Gamma baseline rates and
fixed or normally distributed log risk ratios. Individual zero-event studies
remain in the analysis.

## Install

Install the prepared source archive from the directory where it is saved:

```r
install.packages("metaZIP_0.1.0.tar.gz", repos = NULL, type = "source")
```

After the `package-real-examples` branch is published to GitHub, install it with:

```r
install.packages("remotes") # once, if needed
remotes::install_github("zhangz19/metaZIP", ref = "package-real-examples")
```

R 4.1 or later is required. The package itself uses only R's standard `stats`
and `utils` packages.

## Run the two real examples

```r
library(metaZIP)

data(eus_fna)             # 51 independent studies: 19G versus 22/25G
fit_eus <- with(eus_fna, metaZIPs(y = y, x = x, ni = ni))
fit_eus$out               # effect estimates, SEs, RR intervals, AIC and BIC
fit_eus$paraEst           # alpha, beta, tau(mu), sigma and q
fit_eus$diagnostics       # optimizer, quadrature and uncertainty checks

data(rosiglitazone_mi)    # 48 paired trials: rosiglitazone versus control
fit_mi <- with(rosiglitazone_mi, metaZIPr(y1, n1, y0, n0))
fit_mi$out
fit_mi$diagnostics

coef(fit_eus)             # mean log risk ratios by model
confint(fit_eus)          # exponentiated Wald confidence intervals
logLik(fit_mi, model = "ZIPR")
```

`metaZIPs()` treats each row as an independent single-arm study: `x = 1` is
the treatment group and `x = 0` is the reference. For EUS-FNA these are 19G
and 22/25G needles, respectively. `metaZIPr()` preserves the shared Gamma
baseline within each paired trial; matching positions in its four input
vectors must refer to the same trial.

Both examples are included with the package; no external data files are
needed. Their provenance and column definitions are in `?eus_fna` and
`?rosiglitazone_mi`. The rosiglitazone example uses the 48-trial MI dataset
from [Cai, Parast and Ryan (2010)](https://doi.org/10.1002/sim.3964).

## Models and output

| Model | Gamma baseline rates | Log risk ratio | Extra zero component | Parameters |
| --- | --- | --- | --- | --- |
| POIF | Vary across studies | Common | No | 3 |
| POIR | Vary across studies | Normal random effect | No | 4 |
| ZIPF | Vary across studies | Common | Yes | 4 |
| ZIPR | Vary across studies | Normal random effect | Yes | 5 |

`Estimate` and `coef()` report the mean log risk ratio, **mu**. `RR` is
**exp(mu)**; for random-effect models this is the geometric mean study-level
risk ratio. `sigma` describes between-study variation in log risk ratios.
`SE` measures sampling uncertainty in the estimated mean log risk ratio.

`q` is the probability of the Poisson count component; **1 - q** is the
structural-zero probability. In paired models, the two arms have independent
component indicators with the same `q`, while retaining the shared Gamma
baseline. The former beta parameters `c` and `d` enter these likelihoods only
through `q = c/(c + d)`, so the package estimates `q` directly.

Confidence intervals use observed-information Wald standard errors. When
`sigma = 0` or `q = 1`, the boundary nuisance parameter is held fixed for
that calculation and the interval is labeled accordingly in
`diagnostics$se_type`. Unreliable information or numerical checks produce
unavailable uncertainty estimates; inspect `diagnostics` alongside the
results. Single-arm comparisons rely on between-study comparability and
do not themselves adjust for confounding.

Select individual models with `method = "ZIPR"` or a character vector of
model names. Starting values and numerical controls are documented in
`?metaZIPs`; for example, `control = list(nquad = 81)` changes the initial
quadrature order. Existing calls using legacy `para0` names remain supported.
The old `eps` argument is retained for call compatibility and no longer sets
parameter bounds.

## Relationship to the original analysis

Version 0.1.0 retains the original four Gamma-baseline models and the
`metaZIPs()`/`metaZIPr()` entry points. It adds an installable package,
bundled real data, stable likelihood evaluation, identifiable parameters,
and numerical diagnostics. The EUS-FNA point estimates and log likelihoods
reproduce Table 3 of the original manuscript; its random-model standard
errors change because the package uses uncertainty in the estimated mean,
rather than the heterogeneity standard deviation.

The historical exploratory script remains in `tests/codes.R` in the
repository and is excluded from package builds. The package contains the
real examples; simulation-study code is outside its scope. Manuscript
revisions will be handled separately.

## Citation and license

Run `citation("metaZIP")` for the software citation. Software is provided
under the MIT license, copyright 2022 Zhen Zhang; see [LICENSE.md](LICENSE.md).
Dataset documentation identifies the original sources of the aggregate
study counts.
