\name{rosiglitazone_mi}
\alias{rosiglitazone_mi}
\docType{data}
\title{Myocardial Infarction Counts in 48 Rosiglitazone Trials}
\description{
Observed myocardial infarction (MI) counts and sample sizes from the 48 paired
trials reported in Table I of Cai, Parast and Ryan (2010).
}
\usage{data(rosiglitazone_mi)}
\format{
A data frame with 48 rows and 5 variables:
\describe{
  \item{study_id}{Study number 1 through 48, retained from Table I of the source.}
  \item{y1}{Number of MI events in the rosiglitazone arm.}
  \item{n1}{Sample size of the rosiglitazone arm.}
  \item{y0}{Number of MI events in the control arm.}
  \item{n0}{Sample size of the control arm.}
}
}
\details{
The rosiglitazone arms include 86 MI events among 16856 participants; the
control arms include 72 MI events among 12962 participants. Ten trials have
zero MI events in both arms. All 48 trials are retained without a continuity
correction. Control arms are not uniformly placebo arms.

This is the original 48-trial dataset used by Cai et al., including six
additional trials excluded from the 42-trial analysis of Nissen and Wolski
(2007) because both arms had zero MI events and zero cardiovascular deaths.
It is distinct from the later 56-trial update. The endpoint in this dataset is
MI; death counts also appearing in the source table are not included.
}
\source{
Cai T, Parast L, Ryan L (2010).
Meta-analysis for rare events.
\emph{Statistics in Medicine}, 29(20), 2078--2089.
Table I, page 2083. \doi{10.1002/sim.3964}.
\url{https://pmc.ncbi.nlm.nih.gov/articles/PMC2932857/}

Table I was verified against a full-text copy of the original article at
\url{https://www.academia.edu/143468083/Meta_analysis_for_rare_events}.
}
\references{
Nissen SE, Wolski K (2007).
Effect of rosiglitazone on the risk of myocardial infarction and death from
cardiovascular causes.
\emph{New England Journal of Medicine}, 356(24), 2457--2471.
\doi{10.1056/NEJMoa072761}.
}
\seealso{\code{\link{metaZIPr}}}
\examples{
data(rosiglitazone_mi)
colSums(rosiglitazone_mi[c("y1", "n1", "y0", "n0")])
sum(with(rosiglitazone_mi, y1 == 0 & y0 == 0))
}
\keyword{datasets}
