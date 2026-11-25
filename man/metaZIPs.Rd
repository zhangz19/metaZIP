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
         control = list())

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
\code{tau}, \code{tau(mu)}, and \code{sig} are accepted. Positive legacy
\code{c} and \code{d} values are converted to \code{q = c/(c + d)};
an explicitly supplied \code{q} takes precedence. These are starting
values, not fixed parameter values.}
\item{conf.level}{Confidence level between 0 and 1.}
\item{eps}{Positive compatibility argument. It no longer controls
parameter bounds; supplying a nondefault value produces a warning.
Use \code{control} for numerical settings.}
\item{control}{Named list of numerical settings. See \sQuote{Control settings}.}
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
ZIPF, and ZIPR, respectively. Separate beta mixing parameters \eqn{c,d}
are not estimated because these likelihoods depend only on \eqn{c/(c+d)}.

The reported \code{Estimate} is \eqn{\mu}; \code{RR} is
\eqn{\exp(\mu)}. In random-effect models this is the geometric mean
of study-specific risk ratios, not their arithmetic mean
\eqn{\exp(\mu+\sigma^2/2)}.

Estimation uses multiple starting values, exact nested-model boundary
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
\section{Control settings}{
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
An object of class \code{"metaZIP"}, containing:
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
}}
\references{
Cai T, Parast L, Ryan L (2010). Meta-analysis for rare events.
\emph{Statistics in Medicine}, 29(20), 2078--2089.
\doi{10.1002/sim.3964}.
}
\seealso{\code{\link{metaZIP-methods}}, \code{\link{eus_fna}},
\code{\link{rosiglitazone_mi}}}
\examples{
data(eus_fna)
fit <- with(eus_fna, metaZIPs(y, x, ni, method = "POIF"))
print(fit)
confint(fit)

data(rosiglitazone_mi)
fit_paired <- with(rosiglitazone_mi,
  metaZIPr(y1, n1, y0, n0, method = "POIF"))
summary(fit_paired)
}
