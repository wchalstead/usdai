# Methods for performing l4-norm testing using R
# <To-Do> Implement in C++


# crossUStatL4 calculates the value $\widetilde \theta_{(m+1):k}$, an estimate of the l4-norm using cross U-statistics, given as
# $$\widetilde \theta_{(m+1):k} = \sum_{1 \leq i_1,...,i_{3} \leq m}^{*}\sum_{j = m+1}^k\left(\sum_{l = 1}^pX_{i_1,l}X_{i_2,l}X_{i_{3},l}X_{j,l}\right).$$
# X is nxp matrix
# m is integer cross "threshold" value
# Works since we can exchange the summations and use Newton's identities for symmetric polynomials
#' Title
#'
#' @param X \eqn{n \times p} matrix of data.
#' @param m  Cross threshold value. Should be between 3 and \eqn{n}.
#'
#' @returns The cross U-Statistic for the \eqn{\ell_4}-norm calculated as \deqn{
#' \widetilde \theta_{(m+1):k} = \sum_{1 \leq i_1,...,i_{3} \leq m}^{*}\sum_{j = m+1}^k\left(\sum_{l = 1}^pX_{i_1,l}X_{i_2,l}X_{i_{3},l}X_{j,l}\right)
#' }
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100 * 2), 100, 2)
#' crossUStatL4(X, 50)
crossUStatL4 <- function(X, m) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Check inputs
  if (m > n) {
    stop("m must be less than the number of rows in X")
  }
  if (m < 2) {
    stop("m must be greater than or equal to 3")
  }
  if (m %% 1 != 0) {
    stop("m must be an integer value")
  }

  Xm <- X[1:m, , drop = F]

  # Power sums across the first m rows
  p1 <- colSums(Xm)
  p2 <- colSums(Xm^2)
  p3 <- colSums(Xm^3)

  # e3 (third elementary symmetric polynomial)
  e3 <- (p1^3 - 3*p1*p2 + 2*p3) / 6

  # Sum over last term j = m+1, ..., n
  result <- sum(colSums(X[(m+1):n, , drop = F]) * e3)

  return(result)
}

# Cumulative value of l4 norm estimates as j goes from m+1, ..., n
# X is nxp matrix
# m is integer cross "threshold" value
# Function for calculating cumulative L4 cross stats
# Returns R^(n - m) vector of estimates for each j in (m+1, ..., n)
crossUStatL4_cum <- function(X, m) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Check inputs
  if (m > n) {
    stop("m must be less than the number of rows in X")
  }
  if (m < 2) {
    stop("m must be greater than or equal to 3")
  }
  if (m %% 1 != 0) {
    stop("m must be an integer value")
  }

  # First m rows
  Xm <- X[1:m, , drop = FALSE]

  # Power sums across the first m rows
  p1 <- colSums(Xm)
  p2 <- colSums(Xm^2)
  p3 <- colSums(Xm^3)

  # e3 (third elementary symmetric polynomial)
  e3 <- (p1^3 - 3*p1*p2 + 2*p3) / 6

  # Rows after m
  lastX <- X[(m+1):n, , drop = FALSE]

  # Compute dot products: each row with e3
  dot_vals <- as.vector(lastX %*% e3)

  # Cumulative sum
  cumsum(dot_vals)
}

# Calculate self-normalized value using cumulative stats
# X is nxp matrix
# m is integer cross "threshold" value
crossWStatL4 <- function(X, m){
  n <- nrow(X)

  # Calculate cumulative stats
  cumstats <- crossUStatL4_cum(X,m)

  # Theta estimate is the final value of all the cumulative terms
  theta <- cumstats[n-m]

  # Calculate self-normalizer
  V <- sum((cumstats - (((m+1):n) - m) * cumstats[n-m]/(n - m))^2) / (n - m)

  # Retrun theta^2/v
  return(theta^2/V)
}
