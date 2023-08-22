library(docstring)
library(stats)

preprocess <- function(x, q = 0) {
  #' Preprocessing function for agg methods
  #'
  #' @description This does the preprocessing steps that all the agg methods
  #' have in common.
  #'
  #' @param x A vector of forecasts
  #' @param q The quantile to use for replacing 0s and 1s (between 0 and 1)
  #' @return A vector of forecasts with 0s are replaced by the qth quantile and
  #' 100s are replaced by the (1 - q)th quantile.
  #'
  #' @importFrom stats quantile
  #'
  #' @note Assumes forecasts are in the range 0 to 100, inclusive.

  x <- sort(x)
  x <- x[!is.na(x) & !is.nan(x)]

  if (q != 0 && length(x) > 1) {
    if (all(x == 0)) {
      # Replace them with 10^-12
      return(rep(10^-12, length(x)))
    } else if (all(x == 100)) {
      return(rep(100 - 10^-12, length(x)))
    }
    x[x == 100] <- as.numeric(quantile(x[x != 100], 1 - q, na.rm = TRUE))
    x[x == 0] <- as.numeric(quantile(x[x != 0], q, na.rm = TRUE))
  }

  return(x)
}

trim <- function(x, p = 0.1) {
  #' Trimmed mean
  #'
  #' @description Trim the top and bottom (p*100)% of forecasts
  #'
  #' @param x Vector of forecasts in 0 to 100 range (%)
  #' @param p The proportion of forecasts to trim from each end (between 0 and
  #' 1)
  #' @return (numeric) The trimmed mean of the vector
  #'
  #' @export

  x <- preprocess(x, q = 0)
  trimN <- round(p * length(x))
  lastRow <- length(x) - trimN
  trimVec <- x[(trimN + 1):lastRow]
  trimmedMean <- mean(trimVec)
  return(trimmedMean)
}

hd_trim <- function(x, p = 0.1) {
  #' Highest-Density Trimmed Mean
  #'
  #' @description From Powell et al. (2022) \doi{10.1037/dec0000191}. You find
  #' the shortest interval containing (1-p) * 100% of the data and take the mean
  #' of the forecasts within that interval.
  #'
  #' @note As p gets bigger this acts like a mode in a similar way to
  #' the symmetrically-trimmed mean acting like a median.
  #'
  #' @param x Vector of forecasts in 0 to 100 range (%)
  #' @param p The proportion of forecasts to trim (between 0 and 1)
  #' @return (numeric) The highest-density trimmed mean of the vector
  #'
  #' @export

  x <- preprocess(x, q = 0)
  if (length(x) == 0) {
    return(NA)
  }
  n_out <- floor(length(x) * p)
  n_in <- length(x) - n_out
  d <- c()
  # Get all the intervals of length n_in
  for (i in 1:(n_out + 1)) {
    d[i] <- x[i + n_in - 1] - x[i]
  }
  # Which of those intervals is the smallest?
  i <- which.min(d)
  # Take the mean that starts at that index and goes for n_in
  mean(x[i:(i + n_in - 1)])
}

soften_mean <- function(x, p = 0.1) {
  #' Soften the mean.
  #'
  #' @description If the mean is > .5, trim the top trim%; if < .5, the bottom
  #' trim%. Return the new mean (i.e. soften the mean).
  #'
  #' @param x Vector of forecasts in 0 to 100 range (%)
  #' @param p The proportion of forecasts to trim from each end (between 0
  #' and 1)
  #' @return (numeric) The softened mean of the vector
  #' @note This goes against usual wisdom of extremizing the mean, but performs
  #' well when the crowd has some overconfident forecasters in it.
  #'
  #' @export

  x <- preprocess(x, q = 0)
  mymean <- mean(x)
  if (mymean > .5) {
    return(mean(x[1:floor(length(x) * (1 - p))]))
  } else {
    return(mean(x[ceiling(length(x) * p):length(x)]))
  }
}

neymanAggCalc <- function(x) {
  #' Neyman Aggregation (Extremized)
  #'
  #' @description Takes the arithmetic mean of the log odds of the forecasts,
  #' then extremizes the mean by a factor d, where d is
  #'
  #' (n*(sqrt((3*n^2) - (3*n) + 1) - 2))/(n^2 - n - 1)
  #'
  #' where n is the number of forecasts.
  #'
  #' @param x Vector of forecasts in 0 to 100 range (%)
  #' @return (numeric) The extremized mean of the vector
  #' @references Neyman, E. and Roughgarden, T. (2021). Are you
  #' smarter than a random expert? The robust aggregation of substitutable
  #' signals: \doi{10.1145/3490486.3538243}. Also Jaime Sevilla's EAF post
  #' ``Principled extremizing of aggregated forecasts."
  #'
  #' @export

  x <- preprocess(x, q = 0.05)  # replace 0% forecasts with 5th percentile forecast
  x <- (x / 100)
  n <- length(x)
  d <- (n * (sqrt((3 * n^2) - (3 * n) + 1) - 2)) / (n^2 - n - 1)
  logodds <- log(x / (1 - x))  # log odds from probabilities
  logodds <- mean(logodds) * d  # arithmetic mean, THEN extremize by factor d
  return(100 * exp(logodds) / (1 + exp(logodds)))  # back to probability
}

geoMeanCalc <- function(x, q = 0.05) {
  #' Geometric Mean
  #'
  #' Calculate the geometric mean of a vector of forecasts. We handle 0s by
  #' replacing them with the qth quantile of the non-zero forecasts.
  #'
  #' @param x Vector of forecasts in 0 to 100 range (%)
  #' @param q The quantile to use for replacing 0s (between 0 and 1)
  #' @return (numeric) The geometric mean of the vector
  #' @note agg(a) + agg(not a) does not sum to 1 for this aggregation method.
  #'
  #' @export

  x <- preprocess(x, q)
  geoMean <- exp(mean(log(x)))
  return(geoMean)
}

geoMeanOfOddsCalc <- function(x, q = 0.05, odds = FALSE) {
  #' Geometric Mean of Odds
  #'
  #' Convert probabilities to odds, and calculate the geometric mean of the
  #' odds. We handle 0s by replacing them with the qth quantile of the non-zero
  #' forecasts, before converting.
  #'
  #' @param x A vector of forecasts (probabilities! unless odds = TRUE)
  #' @param q The quantile to use for replacing 0s (between 0 and 1)
  #' @param odds Whether x is already in odds form (TRUE) or probabilities
  #' @return (numeric) The geometric mean of the odds
  #' @note agg(a) + agg(not a) does not sum to 1 for this aggregation method.
  #'
  #' @export

  if (!odds) {
    x <- preprocess(x, q) / 100
    odds <- x / (1 - x)
  } else {
    odds <- x
  }
  geoMeanOfOdds <- exp(mean(log(odds)))
  # Convert back to probability
  geoMeanOfOdds <- geoMeanOfOdds / (geoMeanOfOdds + 1)
  return(geoMeanOfOdds * 100)
}
