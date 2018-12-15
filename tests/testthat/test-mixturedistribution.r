context("test-mixturedistribution")

test_that("mixture distribution works", {
  expect_equal(2, qMixtureDistribution(c(0.5), list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"), 0, 4))
  expect_equal(0.5, pMixtureDistribution(2, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)")))
  x1 <- pMixtureDistribution(c(0.5,1,2,4,5,6), list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  expect_equal(x1[4:6], c(1,1,1))
  x1 <- rMixtureDistribution(10000, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))
  expect_true(abs(mean(x1) - 2) < 0.1)
  expect_warning(rMixtureDistribution(1, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)")))
  expect_true(length(rMixtureDistribution(3, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))) == 3)
  expect_true(length(rMixtureDistribution(4, list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))) == 4)
  expect_equal(1, integrate(dMixtureDistribution, 0, 4, functionList = list("pnorm(x, 1, 0.1)","pnorm(x, 2, 0.1)","pnorm(x, 3, 0.1)"))$value)
  expect_equal(1, integrate(dMixtureDistribution, 0, Inf, functionList = list("pnorm(x, 1, 0.1)","plnorm(x, .2, 0.1)","pnorm(x, 3, 0.1)"))$value)

  f <- function(x) qMixtureDistribution(x, list("pnorm(x, 1, 1)", "pnorm(x, 2, 1)"), xMin = 3, xMax = 20)
  expect_warning(f(0))
  expect_true(f(1) > 14.0)
  expect_error(f(-1))
  expect_error(f(2))

  expect_error(qMixtureDistribution(0.1, list(NA, "pnorm(x, 2, 1)"), xMin = 0, xMax = 5))
  expect_error(qMixtureDistribution(0.1, list("pnorm(y, 1, 1)", "pnorm(x, 2, 1)"), xMin = 0, xMax = 5))
  expect_error(qMixtureDistribution(0.1, list("pnorm(x, 1, 1)", "qnorm(x, 3, 4)"), xMin = 0, xMax = 5))

  expect_warning(qMixtureDistribution(0.0001, list("pnorm(x, 1, 1)", "pnorm(x, 3, 4)"),
                       xMin = 0, xMax = 10, nPoints = 10))
  expect_warning(qMixtureDistribution(0.9999, list("pnorm(x, 1, 1)", "pnorm(x, 3, 4)"),
                       xMin = 0, xMax = 10, nPoints = 10))
  expect_equal(1.5, qMixtureDistribution(0.5, list("pnorm(x, 1, 1)", "pnorm(x, 2, 1)"), xMin = 1, xMax = 5, logScale = TRUE))
  expect_error(qMixtureDistribution(0.5, list("pnorm(x, 1, 1)", "pnorm(x, 2, 1)"), xMin = 0, xMax = 5, logScale = TRUE))

  expect_error(rMixtureDistribution(-1, list("pnorm(x,1,2)","pnorm(x,1,1)")))
  expect_warning(rMixtureDistribution(1, list("pnorm(x,1,2)","pnorm(x,1,1)")))

  expect_equal(0, dMixtureDistribution(0.0001, list("pnorm(x, 1, 0.001)")))
})
