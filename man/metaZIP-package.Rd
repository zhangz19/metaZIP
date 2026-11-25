\name{metaZIP-package}
\alias{metaZIP}
\alias{metaZIP-package}
\docType{package}
\title{Meta-Analysis of Rare Events with Poisson Mixture Models}
\description{
Maximum likelihood estimation of Poisson and zero-inflated Poisson
meta-analysis models with Gamma baseline rates and fixed or normally
distributed log risk ratios. The package supports independent single-arm
studies and paired two-arm trials and includes two real-data examples.
}
\details{
Use \code{\link{metaZIPs}} for independent single-arm studies and
\code{\link{metaZIPr}} for paired trials. Both functions fit the POIF,
POIR, ZIPF, and ZIPR models and return effect estimates, confidence
intervals, information criteria, and numerical diagnostics.

Example data are available as \code{\link{eus_fna}} and
\code{\link{rosiglitazone_mi}}. See \code{\link{metaZIP-methods}} for
display and extraction methods.
}
\author{Zhen Zhang \email{zhangquake1@outlook.com}}
\seealso{\url{https://github.com/zhangz19/metaZIP}}
\keyword{package}
