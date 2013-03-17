\name{MixtureDistribution}
\alias{qMixtureDistribution}
\alias{rMixtureDistribution}
\alias{pMixtureDistribution}
\alias{dMixtureDistribution}
\title{Mixture Distribution functions}
\description{
  pdf, cdf, inverse cdf, and random deviates of a user defined mixture distribution with equal weights.
}
\usage{
dMixtureDistribution(x, functionList)
pMixtureDistribution(x, functionList)
qMixtureDistribution(p, functionList, xMin, xMax, nPoints=1000, logScale=FALSE)
rMixtureDistribution(n, functionList)
}
\arguments{
  \item{x}{The value at which the cdf or pdf are computed}
  \item{functionList}{A list of functions forming the mixture distribution.  The functions must be a continuous CDF dependent on 'x'.}
  \item{p}{A vector of probabilities.}
  \item{xMin}{The minimum value for which probabilities and quantiles are calculated.}
  \item{xMax}{The maximum value for which probabilities and quantiles are calculated.}
  \item{nPoints}{The number of points at which the mixed CDF is calculated.  Default=1000.}
  \item{logScale}{Are the CDF estimates to be performed on a logarithmic scale?  Default=FALSE}
  \item{n}{The number of random deviates to compute}
}
\value{
  A vector of densities, probabilities, quantiles, or random deviates from the mixture distribution.
}
\author{Rob Carnell}
\seealso{\code{\link{pnorm}}}
\examples{
  qMixtureDistribution(c(0.5), list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"), 0, 4)
  pMixtureDistribution(2, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  pMixtureDistribution(c(0.5,1,2,4,5,6), list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  x1 <- rMixtureDistribution(10000, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  mean(x1)
  rMixtureDistribution(1, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  rMixtureDistribution(3, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  rMixtureDistribution(4, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  integrate(dMixtureDistribution, 0, 4, functionList=list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  integrate(dMixtureDistribution, 0, Inf, functionList=list("pnorm(x, 1, 0.1)","plnorm(x, .2, 0.1)","pnorm(x, 3, 0.1)"))
}
\keyword{distribution}


