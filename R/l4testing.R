# Methods for performing l4-norm testing using R
# <To-Do> Implement in C++


# crossUStatL4 calculates the value $\widetilde \theta_{(m+1):k}$, an estimate of the l4-norm using cross U-statistics, given as
# $$\widetilde \theta_{(m+1):k} = \sum_{1 \leq i_1,...,i_{3} \leq m}^{*}\sum_{j = m+1}^k\left(\sum_{l = 1}^pX_{i_1,l}X_{i_2,l}X_{i_{3},l}X_{j,l}\right).$$
# data is nxp matrix
# m is integer cross "threshold" value
# Works since we can exchange the summations and use Newton's identities for symmetric polynomials
#' The Cross U-Statistic for \eqn{\ell_4} Norm Testing
#'
#' @param data \eqn{n \times p} matrix of data.
#' @param m  Cross threshold value. Should be between 3 and \eqn{n}.
#'
#' @returns The cross U-Statistic for the \eqn{\ell_4}-norm calculated as \deqn{
#' \widetilde \theta_{(m+1):k} = \sum_{1 \leq i_1,...,i_{3} \leq m}^{*}\sum_{j = m+1}^k\left(\sum_{l = 1}^pX_{i_1,l}X_{i_2,l}X_{i_{3},l}X_{j,l}\right)
#' }
#' @export
#'
#' @examples
#' set.seed(12)
#' # Generate Data
#' data <- matrix(rnorm(100 * 2), 100, 2)
#' # Calculate U-Statistic
#' crossUStatL4(data, 50)
crossUStatL4 <- function(data, m) {

  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

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

  result <- crossUStatL4_c(data, m)

  return(result)
}



# Calculate self-normalized value using cumulative stats
# data is nxp matrix
# m is integer cross "threshold" value
#' The Self-Normalized Cross U-Statistic for \eqn{\ell_4} Norm Testing
#'
#' @param data \eqn{n \times p} matrix of data.
#' @param m Cross threshold value. Should be between 3 and \eqn{n}.
#'
#' @returns The self-normalized cross U-Statistic for the \eqn{\ell_4}-norm
#' The cross U-Statistic is given as \deqn{
#' \widetilde \theta_{(m+1):k} = \sum_{1 \leq i_1,...,i_{3} \leq m}^{*}\sum_{j = m+1}^k\left(\sum_{l = 1}^pX_{i_1,l}X_{i_2,l}X_{i_{3},l}X_{j,l}\right)
#' }, the self normalizer is given as \deqn{
#' V_{(m+1):n} = \frac{1}{n-m}\sum_{k = m+1}^n\left(\widetilde \theta_{(m+1):k} - \frac{k-m}{n-m}\widetilde \theta_{(m+1):n}\right)^2.
#' }, and the self normalized statistic is given as \deqn{
#' W := \frac{\widetilde \theta_{(m+1):n}^2}{V_{(m+1):n}}
#' }
#' @export
#'
#' @examples
#' set.seed(12)
#' # Generate Data
#' data <- matrix(rnorm(100 * 2), 100, 2)
#' # Calculate U-Statistic
#' crossWStatL4(data, 50)
crossWStatL4 <- function(data, m){

  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

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

  return(crossWStatL4_c(data, m))
}
