\name{sampleSize.randomized}
\alias{sampleSize.randomized}
\title{Sample Size Calculation for Completely Randomized Treatment-Control Designs in Microarray Studies}
\description{
  For any specified power, this routine computes the required sample size \code{n} for completely randomized designs in which differential expression between \code{n} treatment units and \code{n} control units is of interest. The total number of experimental units for the study is \code{2n}.
}
\usage{
  sampleSize.randomized(ER0, G0, power, absMu1, sigmad)
}
\arguments{
  \item{ER0}{mean number of false positives.}
  \item{G0}{anticipated number of genes in the experiment that are not differentially expressed.}
  \item{power}{
specified power level for an individual gene, which represents the expected proportion of differentially expressed genes that will be declared as such by the tests.
}
  \item{absMu1}{
absolute mean difference in log-expression between treatment and control conditions as postulated under the alternative hypothesis H1.
}
  \item{sigmad}{
  anticipated standard deviation of the difference in log-expression between treatment and control conditions. The relation between the standard deviation of the difference (\code{sigmad}) and the experimental error standard deviation (\code{sigma}) is \code{sigmad=sqrt(2)/sigma}.} 
}
\value{
  \item{n}{sample size for each group.} 
  \item{d}{statistical difference between treatment and control conditions under H1 (i.e. \code{d=absMu1/sigmad}).}
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
  Examples and explainations can be found in \url{http://www.biostat.harvard.edu/people/faculty/mltlee/pdf/Web-sampsize-trt-cont-050511r.pdf}.
}
\author{
  Weiliang Qiu (\email{weiliang.qiu@gmail.com}),
  Mei-Ling Ting Lee (\email{meilinglee@sph.osu.edu}),
  George Alex Whitmore (\email{george.whitmore@mcgill.ca})
}
\seealso{ \code{\link{power.randomized}},
  \code{\link{power.matched}},
  \code{\link{power.multi}},
  \code{\link{sampleSize.matched}}
}
\examples{
  sampleSize.randomized(ER0=1, G0=2000, power=0.9, absMu1=1, sigmad=0.566)
}
\keyword{htest}
\keyword{design}
