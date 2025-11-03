# Methods for performing l4-norm testing using R
# <To-Do> Implement in C++


# Calculate value $\widetilde \theta_{(m+1):k}$ given as
# $$\widetilde \theta_{(m+1):k} = \sum_{1 \leq i_1,...,i_{3} \leq m}^{*}\sum_{j = m+1}^k\left(\sum_{l = 1}^pX_{i_1,l}X_{i_2,l}X_{i_{3},l}X_{j,l}\right).$$
# crossUStatL4 calculates an estimate of the l4-norm using cross U-statistics
# X is nxp matrix
# m is integer cross "threshold" value
# Works since we can exchange the summations and use Newton's identities for symmetric polynomials
crossUStatL4 <- function(X, m) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

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
