#' Distribution Functions for a Mixture Distribution
#'
#' @description pdf, cdf, inverse cdf, and random deviates of a user defined mixture distribution with equal weights.
#'
#' @param x The value at which the cdf or pdf are computed
#' @param p A vector of probabilities.
#' @param functionList A list of functions forming the mixture distribution.  The functions must be a continuous CDF dependent on 'x'.
#' @param xMin The minimum value for which probabilities and quantiles are calculated.
#' @param xMax The maximum value for which probabilities and quantiles are calculated.
#' @param nPoints The number of points at which the mixed CDF is calculated.  Default=1000.
#' @param logScale Are the CDF estimates to be performed on a logarithmic scale?  Default=FALSE
#' @param n The number of random deviates to compute
#'
#' @return A vector of densities, probabilities, quantiles, or random deviates from the mixture distribution.
#'
#' @name mixturedist
#'
#' @examples
#'   qMixtureDistribution(c(0.5), list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)",
#'                                     "pnorm(x, 3, 0.1)"), 0, 4)
#'   pMixtureDistribution(2, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)",
#'                                "pnorm(x, 3, 0.1)"))
#'   pMixtureDistribution(c(0.5,1,2,4,5,6),
#'                        list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)",
#'                             "pnorm(x, 3, 0.1)"))
#'   x1 <- rMixtureDistribution(10000, list("pnorm(x, 1, 0.1)",
#'                                          "pnorm(x, 2, 0.1)",
#'                                          "pnorm(x, 3, 0.1)"))
#'   mean(x1)
#'   rMixtureDistribution(1, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)",
#'                                "pnorm(x, 3, 0.1)"))
#'   rMixtureDistribution(3, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)",
#'                                "pnorm(x, 3, 0.1)"))
#'   rMixtureDistribution(4, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)",
#'                                "pnorm(x, 3, 0.1)"))
#'   integrate(dMixtureDistribution, 0, 4,
#'             functionList = list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)",
#'                                 "pnorm(x, 3, 0.1)"))
#'   integrate(dMixtureDistribution, 0, Inf,
#'             functionList = list("pnorm(x, 1, 0.1)","plnorm(x, .2, 0.1)",
#'                                 "pnorm(x, 3, 0.1)"))

.checkFunctions <- function(functionList)
{
  if (any(sapply(functionList, is.na)))
    stop("functionList contains NA values in MixtureDistribution")
  if (!(length(grep("x", functionList)) == length(functionList)))
    stop("one element of functionList does not contain an x")
  if (!(length(grep("^[:blank:]*p", functionList)) == length(functionList)))
    stop("one element of functionList is not a probability function")
}

#' @rdname mixturedist
#' @importFrom stats approx
#' @importFrom assertthat assert_that
#' @export
qMixtureDistribution <- function(p, functionList, xMin, xMax, nPoints=1000,
                                 logScale=FALSE)
{
  assert_that(all(p >= 0) & all(p <= 1),
              msg = "probabilities must be >= 0 and <= 1")
  assert_that(logScale == FALSE || (logScale == TRUE && xMin > 0),
              msg = "xMin must be >0 when logScale = TRUE")
  assert_that(xMin <= xMax, msg = "xMin must be <= xMax")
  .checkFunctions(functionList)

  if (!logScale)
  {
    x <- seq(xMin, xMax, length = nPoints)
  } else
  {
    x <- 10^seq(log10(xMin), log10(xMax), length = nPoints)
  }

  Z <- sapply(functionList, function(z) eval(parse(text = z)))

  Zmean <- apply(Z, 1, mean)

  if (Zmean[1] > min(p))
    warning("probabilities < ", Zmean[1], " will return ", xMin,
            ": min p is ", min(p))
  if (Zmean[length(Zmean)] < max(p))
    warning("probabilities > ", Zmean[length(Zmean)], " will return ", xMax,
            ": max p is ", max(p))

  ret <- approx(x = Zmean, y = x, xout = p, method = "linear",
                yleft = xMin, yright = xMax, ties = mean)

  return(ret$y)
}

#' @rdname mixturedist
#' @export
rMixtureDistribution <- function(n, functionList)
{
  assert_that(n >= 1, msg = paste("n >= 1 is required, ", n, " was supplied"))
  len <- length(functionList)

  if (n < len)
    warning("There are not enough requested samples to sample once from each function")
  .checkFunctions(functionList)

  nOrig <- n
  n <- ceiling(n / length(functionList))

  functionList2 <- gsub("^[:blank:]*p", "r", functionList)
  functionList2 <- gsub("[(][[:blank:]]*x[[:blank:]]*[,]", "(n,", functionList2)

  Z <- sapply(functionList2, function(z) eval(parse(text = z)))

  Z[1:nOrig]
}

#' @rdname mixturedist
#' @export
pMixtureDistribution <- function(x, functionList)
{
  .checkFunctions(functionList)

  Z <- sapply(functionList, function(z) eval(parse(text = z)))

  if (!is.null(dim(Z)))
  {
    Zmean <- apply(Z, 1, mean)
  } else if (length(Z) > 0)
  {
    Zmean <- mean(Z)
  }

  Zmean
}

#' @rdname mixturedist
#' @export
dMixtureDistribution <- function(x, functionList)
{
  .checkFunctions(functionList)

  functionList2 <- gsub("^[:blank:]*p", "d", functionList)

  Z <- sapply(functionList2, function(z) eval(parse(text = z)))

  if (!is.null(dim(Z)))
  {
    Zmean <- apply(Z, 1, mean)
  } else if (length(Z) > 0)
  {
    Zmean <- mean(Z)
  }

  Zmean
}

