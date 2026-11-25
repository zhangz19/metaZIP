\name{metaZIPs}
\alias{metaZIPs}
\alias{metaZIPr}
\title{Fit Rare-Event Meta-Analysis Models}
\description{
Fit four Poisson mixture models to independent single-arm studies or
paired two-arm trials. Gamma baseline rates account for between-study
variation; optional zero inflation and normally distributed log risk
ratios extend the model.
}
\usage{
metaZIPs(y, x, ni, method = c("POIF", "POIR", "ZIPF", "ZIPR"),
         para0 = NULL, conf.level = .95, eps = .Machine$double.eps,
         control = list(), engine = c("current", "legacy"))

metaZIPr(y1, n1, y0, n0, method = c("POIF", "POIR", "ZIPF", "ZIPR"),
         para0 = NULL, conf.level = .95, eps = .Machine$double.eps,
         control = list())
}
\arguments{
\item{y}{Nonnegative integer event counts, one per single-arm study.}
\item{x}{Treatment indicators, coded 0 for reference and 1 for treatment.
Both groups must be represented.}
\item{ni}{Positive integer sample sizes for the single-arm studies.}
\item{y1, y0}{Nonnegative integer event counts in the treatment and
reference arms, respectively. Corresponding elements belong to the same trial.}
\item{n1, n0}{Positive integer sample sizes in the corresponding arms.}
\item{method}{Character vector selecting one or more models:
\code{"POIF"}, \code{"POIR"}, \code{"ZIPF"}, or \code{"ZIPR"}.}
\item{para0}{Optional named numeric vector or list of starting values.
Recognized parameters are \code{alpha}, \code{beta}, \code{mu},
\code{sigma}, and \code{q}. Legacy aliases \code{alph}, \code{bet},
\code{tau}, \code{tau(mu)}, and \code{sig} are accepted. For the current
engine, positive legacy \code{c} and \code{d} values are converted to
\code{q = c/(c + d)}; an explicitly supplied \code{q} takes precedence.
See \sQuote{Historical EUS calculation} for legacy-engine defaults and
parameter mapping. These are starting values, not fixed parameter values.}
\item{conf.level}{Level between 0 and 1 for confidence limits (current
engine) or the historical limits (legacy engine).}
\item{eps}{For the current engine, a positive compatibility argument that
no longer controls parameter bounds; a nondefault value produces a warning.
For the legacy engine, the original positive lower bound for Gamma, mixture,
and ZIPR heterogeneity parameters; it must be strictly between 0 and 1.}
\item{control}{Named list of numerical settings. See \sQuote{Control settings}.
The legacy engine accepts only \code{warn}.}
\item{engine}{For \code{metaZIPs}, \code{"current"} (default) uses the
identifiable parameterization and observed-information uncertainty for the
mean effect. \code{"legacy"} reproduces the original EUS-FNA calculation
with explicitly labelled historical spreads and limits.}
}
\details{
All input vectors must be finite, nonempty, and of equal length, with no
missing values. Counts cannot exceed sample sizes. Individual zero-event
studies are retained, including double-zero paired trials. Each treatment
group must have at least one event overall for a finite fitted log risk ratio.

In \code{metaZIPs}, study \eqn{i} has rate
\deqn{\lambda_i = \xi_i \exp(\tau_i x_i),}
where \eqn{\xi_i} follows a Gamma distribution with shape \eqn{\alpha}
and rate \eqn{\beta}. The Poisson mean is \eqn{n_i\lambda_i}, so
\eqn{\log(n_i)} is an offset. In \code{metaZIPr}, both arms of a
trial share the same \eqn{\xi_i}; the paired likelihood preserves that
dependence.

POIF and ZIPF use a common log risk ratio \eqn{\tau_i=\mu}.
POIR and ZIPR use \eqn{\tau_i\sim N(\mu,\sigma^2)}, independently of
the Gamma baseline. ZIPF and ZIPR add a structural-zero component:
\eqn{q} is the probability of the Poisson count component and
\eqn{1-q} is the structural-zero probability. The paired models use
independent arm-level component indicators with the same \eqn{q}.
The identifiable parameter counts are 3, 4, 4, and 5 for POIF, POIR,
ZIPF, and ZIPR, respectively. The current engine estimates \eqn{q}
directly because these likelihoods depend on beta mixing parameters
\eqn{c,d} only through \eqn{c/(c+d)}.

The reported \code{Estimate} is \eqn{\mu}; \code{RR} is
\eqn{\exp(\mu)}. In random-effect models this is the geometric mean
of study-specific risk ratios, not their arithmetic mean
\eqn{\exp(\mu+\sigma^2/2)}.

The current engine uses multiple starting values, exact nested-model boundary
candidates, and normal quadrature for the random log risk ratio.
Quadrature is checked against a higher-order rule and refined when needed.
Standard errors use the inverse observed information for \eqn{\mu}.
At \eqn{\sigma=0} or \eqn{q=1}, a boundary nuisance parameter is held
fixed when calculating this information; \code{diagnostics$se_type}
identifies such intervals. Failed numerical checks produce unavailable
standard errors and confidence limits rather than a substituted uncertainty
measure. These are asymptotic Wald intervals.

Single-arm comparisons use between-study information. Their interpretation
requires comparable underlying populations or an appropriately justified
common baseline distribution; the model itself does not adjust for confounding.
}
\section{Historical EUS calculation}{
\code{engine = "legacy"} retains the original single-arm likelihoods,
starting values, and optimization settings. It uses 15-node Hermite
quadrature from \pkg{statmod} for POIR and 21-node Hermite quadrature from
\pkg{gaussquad} for ZIPR and the original POIR fallback. These settings are
fixed; current-engine quadrature and optimization controls are rejected.

With \code{para0 = NULL}, the original sequence starts POIF at
\code{alph = 1, bet = 1, tau = 0}, then initializes the other fits from
the preceding fitted models. Supplied legacy starting values may use
\code{alph}, \code{bet}, \code{tau}, \code{mu}, \code{sig}, \code{c},
and \code{d}, or the modern aliases \code{alpha}, \code{beta}, and
\code{sigma}. Missing supplied starts default to 1 for \code{alph},
\code{bet}, \code{sig}, \code{c}, and \code{d}, and 0 for \code{tau}
and \code{mu}; either effect start supplies the other when only one is
present. A supplied \code{q} must be strictly between 0 and 1 and takes
precedence by setting \code{c = q}, \code{d = 1 - q}.

The original function used an inverse-Hessian SE for POIF and ZIPF but
the fitted heterogeneity SD \eqn{\sigma} for POIR and ZIPR. The legacy
engine therefore labels this quantity \code{Legacy spread}, with a
\code{Spread type} column identifying its meaning. Historical limits are
\eqn{\exp(\widehat\mu \mathbin{\pm} z s)}, where \eqn{s} is the legacy
spread and \eqn{z} is the normal quantile for \code{conf.level}. For the
random-effect models these are not confidence limits for the estimated
mean effect. Use the current engine for its observed-information mean-effect
SEs and confidence intervals.

\code{Legacy npara}, \code{Legacy AIC}, and \code{Legacy BIC} retain
the original counts of 3, 4, 5, and 6 parameters. These include both
\eqn{c} and \eqn{d} for the ZIP models, although only their ratio is
identified. They are supplied for historical comparison and are not
exposed as model degrees of freedom through \code{logLik}.

For the bundled EUS data, the legacy output closely reproduces the
original console results. ZIPR has a legacy spread of approximately
0.715695 and historical limits of approximately (0.318961, 5.273990).
The original manuscript's Table 3 entry of 0.64 and (0.37, 4.55) is not
reproduced by this function. Small platform-dependent differences are
possible, particularly in ZIPF's inverse-Hessian spread. A failed
inverse-Hessian calculation returns an unavailable spread rather than
substituting another quantity.
}
\section{Control settings}{
These controls apply to the current engine. The legacy engine accepts only
\code{warn}, which controls warnings for failed optimization or unavailable
historical spreads; its numerical settings remain fixed.
\describe{
\item{\code{nquad}}{Initial number of normal quadrature nodes;
integer from 15 to 301, default 61.}
\item{\code{maxit}}{Maximum optimizer iterations per start;
positive integer, default 1500.}
\item{\code{reltol}}{Positive relative optimizer tolerance, default \code{1e-10}.}
\item{\code{quad_tol}}{Positive maximum accepted absolute log-likelihood
difference between quadrature rules, default \code{1e-5}.}
\item{\code{max_refine}}{Maximum quadrature refinements;
integer from 0 to 3, default 2. Each refinement changes \eqn{k} nodes
to \eqn{2k+1}.}
\item{\code{sigma_max}}{Positive upper numerical bound for \eqn{\sigma}, default 5.}
\item{\code{alpha_min}, \code{alpha_max}}{Positive lower and upper
numerical bounds for \eqn{\alpha}, defaults \code{1e-6} and \code{1e6}.}
\item{\code{mu_bound}}{Positive absolute numerical bound for \eqn{\mu}, default 20.}
\item{\code{boundary_tol}}{Positive tolerance for identifying boundary
solutions, default \code{1e-5}.}
\item{\code{warn}}{Whether to warn about unsuccessful numerical or
uncertainty checks, default \code{TRUE}. Diagnostics are always returned.}
}}
\value{
An object of class \code{"metaZIP"} for the current engine, containing:
\describe{
\item{\code{out}}{Data frame with one row per requested model: log likelihood,
mean log risk ratio (\code{Estimate}), its standard error (\code{SE}),
exponentiated confidence limits, parameter count (\code{npara}), AIC, BIC,
and \code{RR}. AIC and BIC use the nominal identifiable parameter count,
including when a parameter estimate lies on a boundary.}
\item{\code{paraEst}}{Estimates of Gamma shape \code{alpha}, Gamma rate
\code{beta}, mean log risk ratio \code{tau(mu)}, heterogeneity standard
deviation \code{sigma}, and count-component probability \code{q}.}
\item{\code{diagnostics}}{Convergence code and message, boundary flags,
information-matrix check, type of standard error, quadrature difference,
quadrature check, and final number of nodes.}
\item{\code{fits}}{Individual fitted-model objects used by extraction methods.}
\item{\code{design}, \code{nstudies}, \code{data}}{Study design, number of
independent studies or paired trials, and analyzed inputs.}
\item{\code{conf.level}, \code{control}, \code{call}, \code{estimand}}{
Confidence level, final numerical controls, original call, and estimand label.}
}

The legacy engine returns an object of class \code{"metaZIP_legacy"}:
\describe{
\item{\code{out}}{Log likelihood, \code{Estimate}, \code{Legacy spread},
\code{Legacy lower}, \code{Legacy upper}, \code{Legacy npara},
\code{Legacy AIC}, \code{Legacy BIC}, \code{RR}, and \code{Spread type}.
The spread and limit columns follow the historical definitions above.}
\item{\code{paraEst}}{The original six parameter columns: \code{alpha},
\code{beta}, \code{tau(mu)}, \code{sigma}, \code{c}, and \code{d}.
Parameters absent from a model are \code{NA}.}
\item{\code{diagnostics}}{Optimizer convergence code, convergence indicator,
message, spread type, and spread availability. These are historical-fit
diagnostics, not the current engine's quadrature or information checks.}
\item{\code{fits}}{Raw fitted \code{optim} objects for the requested models.}
}
}
\references{
Cai T, Parast L, Ryan L (2010). Meta-analysis for rare events.
\emph{Statistics in Medicine}, 29(20), 2078--2089.
\doi{10.1002/sim.3964}.
}
\seealso{\code{\link{metaZIP-methods}}, \code{\link{eus_fna}},
\code{\link{rosiglitazone_mi}}}
\examples{
data(eus_fna)
historical <- with(eus_fna, metaZIPs(y, x, ni, engine = "legacy"))
print(historical)
coef(historical)

fit <- with(eus_fna, metaZIPs(y, x, ni, method = "POIF"))
print(fit)
confint(fit)

data(rosiglitazone_mi)
fit_paired <- with(rosiglitazone_mi,
  metaZIPr(y1, n1, y0, n0, method = "POIF"))
summary(fit_paired)
}
