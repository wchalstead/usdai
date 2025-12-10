# Variance testing test statistic
# data nxp matrix of data
# m is integer cross "threshold" value
# Sigma is null hypothesis covariance matrix
#' Title
#'
#' @param data
#' @param m
#' @param Sigma
#'
#' @returns
#' @export
#'
#' @examples
crossUStatVar <- function(data, m, Sigma) {

  data <- as.matrix(data)
  Sigma <- as.matrix(Sigma)
  n <- nrow(data)
  p <- ncol(data)
  pSigma <- ncol(Sigma)

  # Check inputs
  if (m > n) {
    stop("m must be less than the number of rows in data")
  }
  if (m < 2) {
    stop("m must be greater than or equal to 3")
  }
  if (m %% 1 != 0) {
    stop("m must be an integer value")
  }
  if (p != pSigma) {
    stop(paste("Sigma should be of dimension", p, 'x', p))
  }
  if (!isSymmetric(Sigma)) {
    stop("Sigma must be a symmetrix, positive semi-definite matrix")
  }
  if (all(eigen(Sigma)$values >= 0)) {
    stop("Sigma must be a symmetrix, positive semi-definite matrix")
  }

  return(tail(crossUStatVar_cum_c(data, m, Sigma), 1))
}


# Calculate normalized W stat using cumulative stats
#' Title
#'
#' @param data
#' @param m
#' @param Sigma
#'
#' @returns
#' @export
#'
#' @examples
crossWStatVar <- function(data, m, Sigma){
  n <- nrow(data)

  # Calculate cumulative stats
  cumstats <- crossUStatVar_cum(data, m, Sigma)
  theta <- cumstats[n-m]

  # Calculate self-normalizer
  V <- sum((cumstats - (((m+1):n) - m) * cumstats[n-m]/(n - m))^2) / (n - m)
  return(theta^2/V)
}
