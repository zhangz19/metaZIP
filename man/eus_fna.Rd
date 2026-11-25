\name{eus_fna}
\alias{eus_fna}
\docType{data}
\title{Complications in EUS-FNA Single-Arm Studies}
\description{
Observed complication counts and sample sizes for 51 single-arm studies of
endoscopic ultrasound-guided fine-needle aspiration (EUS-FNA), comparing
19-gauge needles with 22- or 25-gauge needles.
}
\usage{data(eus_fna)}
\format{
A data frame with 51 rows and 5 variables:
\describe{
  \item{study_id}{Integer display identifier, assigned sequentially within this
  dataset. These are not original study identifiers.}
  \item{needle}{Needle group: \code{"19G"} or \code{"22/25G"}.}
  \item{x}{Treatment indicator: 1 for 19G and 0 for 22/25G.}
  \item{y}{Number of participants with a reported complication.}
  \item{ni}{Study sample size.}
}
}
\details{
The six 19G studies include 3 events among 215 participants. The 45 22/25G
studies include 74 events among 5957 participants. Four and 27 studies,
respectively, report zero events.

Rows were transcribed from Table 1 of the author-supplied manuscript, reading
each column block from top to bottom and then moving from left to right.
Original study names are not provided in that table. All observed counts,
including zeros, are retained without a continuity correction. Each row is a
separate single-arm study; rows from different needle groups are not paired.
}
\source{
Chawla A, Maiti T, Pradhan V, Zhang Z.
\emph{Meta-analysis of rare events with application to medical data}.
Author-supplied manuscript, Table 1 (\file{ZIPR26.tex}; also present in
\file{ZIPR21.pdf}).
}
\seealso{\code{\link{metaZIPs}}}
\examples{
data(eus_fna)
aggregate(cbind(y, ni) ~ needle, data = eus_fna, FUN = sum)
}
\keyword{datasets}
