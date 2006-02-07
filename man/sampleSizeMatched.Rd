\name{sampleSize.matched}
\alias{sampleSize.matched}
\title{
  Sample Size Calculation for Matched-Pairs Designs in Microarray Studies
}
\description{
  This routine computes the sample size \code{n} required to achieve a specified power level for a matched-pairs design in which differential expression between \code{n} treatment units and \code{n} matched control units is of interest. The total number of experimental units for the study is \code{2n}.}
\usage{
  sampleSize.matched(ER0, G0, power, absMu1, sigmad)
}
\arguments{
  \item{ER0}{mean number of false positives.}
  \item{G0}{anticipated number of genes in the experiment that are not differentially expressed.}
  \item{power}{
specified power level for an individual gene, which represents the expected
proportion of differentially expressed genes that will be declared as such by the tests.
}
  \item{absMu1}{
absolute mean difference in log-expression between treatment and control units as postulated under the alternative hypothesis H1.
}
  \item{sigmad}{
anticipated standard deviation of the difference in log-expression between matched treatment and control units.}
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
  Examples and explainations can be found in \url{
http://www.biostat.harvard.edu/people/faculty/mltlee/pdf/Web-sampsize-matched050510.pdf}.
}
\author{
  Weiliang Qiu (\email{weiliang.qiu@gmail.com}),
  Mei-Ling Ting Lee (\email{meilinglee@sph.osu.edu}),
  George Alex Whitmore (\email{george.whitmore@mcgill.ca})
}
\seealso{ \code{\link{power.randomized}},
  \code{\link{power.matched}}
  \code{\link{power.multi}},
  \code{\link{sampleSize.randomized}} 
}
\examples{
  sampleSize.matched(ER0=1, G0=2000, power=0.9, absMu1=1, sigmad=0.5)
}
\keyword{htest}
\keyword{design}

