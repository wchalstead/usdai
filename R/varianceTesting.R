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
  n <- nrow(data); p <- ncol(data)
  if (m < 3) stop("m must be >= 3")
  if (n <= m) stop("need n > m")
  Sigma <- (Sigma + t(Sigma)) / 2

  # --- upper-tri indices and weight vector
  ut <- upper.tri(Sigma, diag = TRUE)
  q <- sum(ut)                 # number of stored entries
  # weight: diagonal 1, off-diagonal 2 as in Frobenius inner product
  w_mat <- matrix(2, p, p); diag(w_mat) <- 1
  w_flat_sqrt <- sqrt(w_mat[ut])   # multiply flattened entries by this

  # --- All pairs (i1,i2) with i1 < i2 from 1..m
  twopairs <- combn(m - 1, 2)
  npairs <- ncol(twopairs)
  # compute difference vectors d = x_i1 - x_i2 for each pair
  D12 <- data[twopairs[1, ], , drop = FALSE] - data[twopairs[2, ], , drop = FALSE]  # npairs x p

  # Flattened A matrices (one row per pair), flattened upper-tri and weighted
  # A_pair = tcrossprod(d) - Sigma ; flattened and multiplied by w_flat_sqrt
  A_flat <- matrix(0, nrow = npairs, ncol = q)
  for (r in seq_len(npairs)) {
    mat <- tcrossprod(D12[r, ]) - Sigma
    A_flat[r, ] <- mat[ut] * w_flat_sqrt
  }

  # --- For each i3 in 3..m compute block_sum_i3 := sum_{j=m+1..n} (tcrossprod(x_i3 - x_j) - Sigma)
  i3_vals <- 3:m
  n_i3 <- length(i3_vals)
  n_j <- n - m
  block_flat <- matrix(0, nrow = n_i3, ncol = q)  # one row per i3

  rows_j <- (m+1):n
  for (idx in seq_len(n_i3)) {
    i3 <- i3_vals[idx]
    # D matrix of differences: rows are j in (m+1):n
    D <- sweep(data[rows_j, , drop = FALSE], 2, data[i3, ], "-")  # (n-m) x p
    # sum_j tcrossprod(d_j) = t(D) %*% D  (p x p)
    S_mat <- crossprod(D)        # crossprod(D) == t(D) %*% D
    # subtract Sigma for each j: sum_j Sigma = (n-m) * Sigma
    block_mat <- S_mat - (n_j) * Sigma
    block_flat[idx, ] <- block_mat[ut] * w_flat_sqrt
  }

  # --- Now compute cumulative sums over i3 in decreasing order
  # For each i2 (second index in twopairs), we need sum over i3 > i2 of block_mat
  # Precompute cumulative sums: cum_from_idx[t] = sum_{i3_idx >= t} block_flat[i3_idx, ]
  # We'll create cum_by_i3 such that cum_by_i3[idx] = sum_{i3' >= i3_vals[idx]} block_flat[i3', ]
  cum_rev <- matrix(apply(block_flat, 2, function(col) rev(cumsum(rev(col)))), n_i3, q)
  # cum_rev has same dim n_i3 x q; cum_rev[idx, ] gives sum_{i3' >= i3_vals[idx]} block_flat[i3', ]

  # For each pair (i1,i2) we need sum_{i3 > i2} block_flat
  # For a given i2, find the first index in i3_vals satisfying i3_vals > i2:
  # idx_i3_first <- which(i3_vals > i2)[1]  (if NULL -> sum is zero)
  # We'll group pairs by their i2 for vectorization
  i2_vals <- twopairs[2, ]
  unique_i2 <- sort(unique(i2_vals))

  stat <- 0
  for (t in unique_i2) {
    pair_idx <- which(i2_vals == t)   # indices of pairs whose second index is t
    pos <- which(i3_vals > t)
    if (length(pos) == 0L) next       # no i3 > t -> contribute 0
    first_pos <- pos[1]               # index into cum_rev
    cum_row <- cum_rev[first_pos, ]   # 1 x q
    # sum over pairs in this group of A_flat[pair_idx, ] %*% cum_row
    # vectorized dot products:
    # compute row-wise sum of A_flat[pair_idx, ] * matrix(rep(cum_row,...), nrow=length(pair_idx), byrow=TRUE)
    blockA <- A_flat[pair_idx, , drop = FALSE]
    stat <- stat + sum(blockA * matrix(rep(cum_row, each = nrow(blockA)), nrow = nrow(blockA), byrow = FALSE))
  }

  return(stat / 4)
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
