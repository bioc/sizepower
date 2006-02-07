\name{power.multi}
\alias{power.multi}
\title{
Power Calculations for Multiple Treatments Design
with an Isolated Treatment Effect in Microarray Studies
}
\description{
  Assume \code{numTrt} treatment conditions are being studied in either a completely randomized or randomized block design. Under the alternative hypothesis H1, one treatment is distinguished from the other \code{numTrt - 1} treatments by exhibiting differential expression for the gene. This computer routine calculates the individual power value for the design. This power value is the expected fraction of truly differentially expressed genes that will be correctly declared as differentially expressed by the tests.
}
\usage{
  power.multi(ER0, G0, numTrt, absMu1, sigma, n)
}
\arguments{
  \item{ER0}{mean number of false positives.}
  \item{G0}{anticipated number of genes in the experiment that are not differentially expressed.}
  \item{numTrt}{total number of treatment conditions.}
  \item{absMu1}{
  the absolute difference in expression between the distinguished treatment and the other treatments on the log-intensity scale.
  }
  \item{sigma}{
 anticipated experimental error standard deviation of the difference in log-expression between treatments.
  }
  \item{n}{the sample size for each group.}
}
\value{
  \item{power}{power.} 
  \item{psi1}{non-centrality parameter.}
}
\references{
  Lee, M.-L. T. (2004).
  Analysis of Microarray Gene Expression Data.
  \emph{Kluwer Academic Publishers, ISBN 0-7923-7087-2.}

  Lee, M.-L. T., Whitmore, G. A. (2002).
  Power and sample size for DNA microarray studies.
  \emph{Statistics in Medicine, \bold{21}:3543-3570.}
}
\note{
  Examples and explainations can be found in \url{
http://www.biostat.harvard.edu/people/faculty/mltlee/pdf/Web-power-isolated050510.pdf}.
}
\author{
  Weiliang Qiu (\email{weiliang.qiu@gmail.com}),
  Mei-Ling Ting Lee (\email{meilinglee@sph.osu.edu}),
  George Alex Whitmore (\email{george.whitmore@mcgill.ca})
}
\seealso{ \code{\link{power.randomized}},
  \code{\link{power.matched}},
  \code{\link{sampleSize.randomized}},
  \code{\link{sampleSize.matched}}
}
\examples{
  power.multi(ER0=2, G0=10000, numTrt=6, absMu1=0.585, sigma=0.3, n=8)
}
\keyword{htest}
\keyword{design}

