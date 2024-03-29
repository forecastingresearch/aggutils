library(testthat)

forecasts <- c(rep(.1, 5), rep(.2, 5), rep(.4, 5), rep(.5, 2), .8, .9, 1.0)

test_that("trimmed mean works", {
  # Untrimmed, the mean is .36. Trimmed, it's .31875.
  expect_equal(round(trim(forecasts), 2), .32)
})

test_that("trimmed mean works with p arg", {
  # Same vector but trimming top and bottom 25% (leaving middle 50%).
  expect_equal(trim(forecasts, .25), .3)
})

test_that("P(a) + P(not a) = 1 always", {
  not_forecasts <- 1 - forecasts
  expect_equal(trim(forecasts) + aggutils::trim(not_forecasts), 1)
  expect_equal(neymanAggCalc(forecasts) + neymanAggCalc(not_forecasts), 1)
  expect_equal(hd_trim(forecasts) + hd_trim(not_forecasts), 1)
  # GeoMean doesn't have this property
})

test_that("Highest isn't more than 10x lowest", {
  # Make a vector containing each of the aggregation methods on forecasts.
  # The highest value should be no more than 10x the lowest.
  agg_vec <- c(
    trim(forecasts), hd_trim(forecasts),
    neymanAggCalc(forecasts), geoMeanCalc(forecasts),
    geoMeanOfOddsCalc(forecasts)
  )
  expect_equal(max(agg_vec) / min(agg_vec) < 10, TRUE)
})

test_that("It works when they're all 0", {
  # Trim surely works
  forecasts <- rep(0, 100)
  expect_equal(trim(forecasts), 0)
  # What about geometric mean of odds?
  expect_equal(geoMeanOfOddsCalc(forecasts), 0)
})

test_that("Geo mean of odds of (a, b, c) and of (1-a, 1-b, 1-c) sums to 1", {
  vecGMOD <- geoMeanOfOddsCalc(c(.7, .2, .1))
  recipVecGMOD <- geoMeanOfOddsCalc(c(.3, .8, .9))
  expect_equal(vecGMOD + recipVecGMOD, 1)
})

test_that("Same but directly provide odds", {
  vecGMOD <- geoMeanOfOddsCalc(c(3. / 1., 1. / 2., 5. / 2.), odds = TRUE)
  recipVecGMOD <- geoMeanOfOddsCalc(c(1. / 3., 2. / 1., 2. / 5.), odds = TRUE)
  expect_equal(vecGMOD + recipVecGMOD, 1)
})

test_that("Geom Mean (Probabilities) reciprocal test (sum to less than...)", {
  vec <- geoMeanCalc(c(.7, .2, .1))
  recipVec <- geoMeanCalc(c(.3, .8, .9))
  expect_lt(vec + recipVec, 1)
})
