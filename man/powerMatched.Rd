\name{power.matched}
\alias{power.matched}
\title{
  Power Calculations for Matched-Pairs Designs in Microarray Studies
}
\description{
  This routine computes the individual power value for a matched-pairs design having \code{n} treatment units and \code{n} matched control units. This power value is the expected fraction of truly differentially expressed genes that will be correctly declared as differentially expressed by the tests.
}
\usage{
  power.matched(ER0, G0, absMu1, sigmad, n)
}
\arguments{
  \item{ER0}{mean number of false positives.}
  \item{G0}{
  anticipated number of genes in the experiment that are not differentially expressed.}
  \item{absMu1}{
absoulte mean difference in log-expression between treatment and control conditions as postulated under the alternative hypothesis H1.
  }
  \item{sigmad}{
anticipated standard deviation of the difference in log-expression between matched treatment and control units. The relation between the standard deviation of the difference (\code{sigmad}) and the experimental error standard deviation (\code{sigma}) is \code{sigmad=sqrt(2)/sigma}.}
  \item{n}{the sample size for each group.}
}
\value{
  \item{power}{power.} 
  \item{psi1}{non-centrality parameter.}
}
\references{
  Lee, M.-L. T. (2004).
  Analysis of Microarray Gene Expression Data.
  \emph{Kluwer Academic Publishers}, ISBN 0-7923-7087-2.

  Lee, M.-L. T., Whitmore, G. A. (2002).
  Power and sample size for DNA microarray studies.
  \emph{Statistics in Medicine}, \bold{21}:3543-3570.
}
\note{
  Examples and explainations can be found in \url{
http://www.biostat.harvard.edu/people/faculty/mltlee/pdf/Web-power-matched050510.pdf}.
}
\author{
  Weiliang Qiu (\email{weiliang.qiu@gmail.com}),
  Mei-Ling Ting Lee (\email{meilinglee@sph.osu.edu}),
  George Alex Whitmore (\email{george.whitmore@mcgill.ca})
}
\seealso{ \code{\link{power.randomized}},
  \code{\link{power.multi}},
  \code{\link{sampleSize.randomized}}, 
  \code{\link{sampleSize.matched}}
}
\examples{
  power.matched(ER0=2, G0=5000, absMu1=1, sigmad=0.4243, n=4)
}
\keyword{htest}
\keyword{design}

