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
  if (all(eigen(Sigma)$values > 0)) {
    stop("Sigma must be a symmetrix, positive semi-definite matrix")
  }

  return(tail(crossUStatVar_cum_c(data, m, Sigma), 1))
}


# Version to give cumulative stats over time
# data nxp matrix of data
# m is integer cross "threshold" value
# Function for calculating cumulative matrix difference norm cross stats
# Returns R^(n - m) vector of estimates for each j in (m+1, ..., n)
crossUStatVar_cum <- function(data, m, Sigma) {
  n <- nrow(data)
  p <- ncol(data)
  nm <- n - m
  q <- p * (p + 1) / 2  # number of unique entries in symmetric matrix

  # helper: flatten symmetric matrix with sqrt(2) scaling
  sym_flat <- function(M) {
    ut <- upper.tri(M)
    M[ut] <- M[ut] * sqrt(2)
    M[upper.tri(M, diag = TRUE)]
  }

  # Step 1: Precompute A_{i1,i2}
  twopairs <- combn(m-1, 2)
  n_two <- ncol(twopairs)
  A <- matrix(0, q, n_two)
  for (i in seq_len(n_two)) {
    diff_vec <- data[twopairs[1, i], ] - data[twopairs[2, i], ]
    M <- tcrossprod(diff_vec) - Sigma
    A[, i] <- sym_flat(M)
  }

  # Step 2: Precompute B_{i3,j}
  i3_seq <- 3:m
  n_i3 <- length(i3_seq)
  B <- array(0, dim = c(q, n_i3, nm))

  for (t in seq_along(i3_seq)) {
    i3 <- i3_seq[t]
    for (j_idx in seq_len(nm)) {
      j <- m + j_idx
      diff_vec <- data[i3, ] - data[j, ]
      M <- tcrossprod(diff_vec) - Sigma
      B[, t, j_idx] <- sym_flat(M)
    }
  }

  # Step 3: Compute cumulative statistic
  results <- numeric(nm)
  cumB <- matrix(0, q, n_i3)


  for (j_idx in seq_len(nm)) {
    cumB <- cumB + B[, , j_idx]
    total_sum <- 0
    for (i in seq_len(n_two)) {
      valid_i3 <- which(i3_seq > twopairs[2, i])
      if (length(valid_i3) > 0) {
        total_sum <- total_sum + sum(A[, i] * rowSums(cumB[, valid_i3, drop = FALSE]))
      }
    }
    results[j_idx] <- total_sum
  }

  return(results/4)
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
