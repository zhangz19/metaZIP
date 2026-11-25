\name{metaZIP-methods}
\alias{metaZIP-methods}
\alias{print.metaZIP}
\alias{summary.metaZIP}
\alias{coef.metaZIP}
\alias{confint.metaZIP}
\alias{logLik.metaZIP}
\alias{vcov.metaZIP}
\alias{print.metaZIP_legacy}
\alias{summary.metaZIP_legacy}
\alias{coef.metaZIP_legacy}
\alias{confint.metaZIP_legacy}
\alias{logLik.metaZIP_legacy}
\alias{vcov.metaZIP_legacy}
\title{Display and Extract metaZIP Results}
\description{
Display fitted models and numerical diagnostics or extract effect estimates,
confidence intervals, likelihoods, and covariance matrices.
}
\usage{
\method{print}{metaZIP}(x, digits = 4, ...)
\method{summary}{metaZIP}(object, ...)
\method{coef}{metaZIP}(object, ...)
\method{confint}{metaZIP}(object, parm, level = object$conf.level,
                         exponentiate = TRUE, ...)
\method{logLik}{metaZIP}(object, model = NULL, ...)
\method{vcov}{metaZIP}(object, model = NULL, ...)

\method{print}{metaZIP_legacy}(x, digits = 4, ...)
\method{summary}{metaZIP_legacy}(object, ...)
\method{coef}{metaZIP_legacy}(object, ...)
\method{confint}{metaZIP_legacy}(object, ...)
\method{logLik}{metaZIP_legacy}(object, ...)
\method{vcov}{metaZIP_legacy}(object, ...)
}
\arguments{
\item{x, object}{An object returned by \code{\link{metaZIPs}} or
\code{\link{metaZIPr}}.}
\item{digits}{Number of significant digits used when printing.}
\item{parm}{Optional fitted model names or their numeric indices.
By default, intervals are returned for all fitted models.}
\item{level}{Confidence level between 0 and 1.}
\item{exponentiate}{Whether to return confidence limits on the risk-ratio
scale. Set to \code{FALSE} for the mean log risk ratio.}
\item{model}{One fitted model name. Required for \code{logLik} and
\code{vcov} when more than one model is present.}
\item{...}{Additional arguments passed to printing methods, or unused by
extraction methods.}
}
\details{
\code{coef} returns the mean log risk ratio, not its exponential, for
both current and legacy objects.
For current objects, \code{confint} calculates Wald intervals from the stored estimates and
standard errors; it does not refit a model. The same boundary and numerical
qualifications as those reported in \code{object$diagnostics} apply.

\code{vcov} returns the covariance matrix in the fitted internal
parameterization: \code{log_alpha}, \code{log_rate} (the logarithm of
\eqn{\alpha/\beta}), \code{mu}, and the active \code{sigma} and
\code{q} parameters. Boundary nuisance parameters held fixed for
information-based inference are omitted. Entries are unavailable when the
information calculation fails its numerical checks. This is not the
covariance matrix of all columns of \code{object$paraEst}.

Objects from \code{metaZIPs(engine = "legacy")} have a separate class.
Their print and summary methods label historical spreads and limits:
POIF and ZIPF use inverse-Hessian SEs, whereas POIR and ZIPR use the
heterogeneity SD. The random-model limits are not confidence limits for
the estimated mean. \code{confint} and \code{vcov} for legacy objects
stop with guidance to use the current engine for mean-effect uncertainty.
\code{logLik} also stops: historical parameter counts include redundant
\code{c,d} parameters. The numeric likelihoods and historical criteria
remain available in \code{object$out}.
}
\value{
\code{print} and \code{summary} invisibly return the supplied object.
\code{coef} returns a named numeric vector. \code{confint} returns a
two-column matrix. \code{logLik} returns an object of class
\code{"logLik"} with nominal identifiable parameter count \code{df}
and the number of independent studies or paired trials as \code{nobs}.
\code{vcov} returns a covariance matrix.
For legacy objects, only print, summary, and coefficient extraction return
results; the other three methods deliberately report the unsupported
interpretation as an error.
}
\seealso{\code{\link{metaZIPs}}, \code{\link{metaZIPr}}}
\examples{
data(eus_fna)
historical <- with(eus_fna, metaZIPs(y, x, ni, engine = "legacy"))
print(historical)
coef(historical)

fit <- with(eus_fna, metaZIPs(y, x, ni, method = "POIF"))
coef(fit)
confint(fit)
confint(fit, exponentiate = FALSE)
logLik(fit)
vcov(fit)
}
