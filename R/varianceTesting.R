# Variance testing test statistic
# data nxp matrix of data
# m is integer cross "threshold" value
# Sigma is null hypothesis covariance matrix
#' The Cross U-Statistic for Variance Testing
#'
#' @param data \eqn{n \times p} matrix of data.
#' @param m Cross threshold value. Should be between 3 and \eqn{n}.
#' @param Sigma Symmetric, positive semi-definite matrix to serve as the null hypothesis
#'
#' @returns The Cross U-Statistic for Variance Testing defined as \deqn{
#' \widetilde \theta_{(m+1):k} = \sum_{1 \leq i_1,...,i_{3} \leq m}^{*}\sum_{j = m+1}^k\sum_{l_1,l_2 = 1}^p\frac{1}{4}[(X_{i_1,l_1} - X_{i_2,l_1})(X_{i_1,l_2} - X_{i_2,l_2}) - 2\Sigma_0(l_1,l_2)] \\ \cdot [(X_{i_3,l_1} - X_{j,l_1})(X_{i_3,l_2} - X_{j,l_2}) - 2\Sigma_0(l_1,l_2)].
#' }
#' @export
#'
#' @examples
#' set.seed(12)
#' # Generate Data
#' data <- matrix(rnorm(100 * 2), 100, 2)
#' # Calculate U-Statistic
#' crossUStatVar(data, 50, diag(2))
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
    stop("Sigma must be a symmetric, positive semi-definite matrix")
  }
  if (!all(eigen(Sigma)$values >= 0)) {
    stop("Sigma must be a symmetric, positive semi-definite matrix")
  }

  return(tail(crossUStatVar_cum_c(data, m, Sigma), 1))
}


# Calculate normalized W stat using cumulative stats
#' The Self-Normalized Cross U-Statistic for Variance Testing
#'
#' @param data \eqn{n \times p} matrix of data.
#' @param m Cross threshold value. Should be between 3 and \eqn{n}.
#' @param Sigma Symmetric, positive semi-definite matrix to serve as the null hypothesis
#'
#' @returns The self-normalized cross U-Statistic for variance testing. The Cross U-Statistic for Variance Testing is defined as \deqn{
#' \widetilde \theta_{(m+1):k} = \sum_{1 \leq i_1,...,i_{3} \leq m}^{*}\sum_{j = m+1}^k\sum_{l_1,l_2 = 1}^p\frac{1}{4}[(X_{i_1,l_1} - X_{i_2,l_1})(X_{i_1,l_2} - X_{i_2,l_2}) - 2\Sigma_0(l_1,l_2)] \\ \cdot [(X_{i_3,l_1} - X_{j,l_1})(X_{i_3,l_2} - X_{j,l_2}) - 2\Sigma_0(l_1,l_2)]
#' ,} the self normalizer is given as \deqn{
#' V_{(m+1):n} = \frac{1}{n-m}\sum_{k = m+1}^n\left(\widetilde \theta_{(m+1):k} - \frac{k-m}{n-m}\widetilde \theta_{(m+1):n}\right)^2
#' ,} and the self normalized statistic is given as \deqn{
#' W := \frac{\widetilde \theta_{(m+1):n}^2}{V_{(m+1):n}}
#' .}
#' @export
#'
#' @examples
#' set.seed(12)
#' # Generate Data
#' data <- matrix(rnorm(100 * 2), 100, 2)
#' # Calculate U-Statistic
#' crossWStatVar(data, 50, diag(2))
crossWStatVar <- function(data, m, Sigma){
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
  if (!all(eigen(Sigma)$values >= 0)) {
    stop("Sigma must be a symmetrix, positive semi-definite matrix")
  }


  return(crossWStatVar_c(data, m, Sigma))
}
