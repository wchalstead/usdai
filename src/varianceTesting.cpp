#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat crossUStatVar_cum_c(const arma::mat& X,
                              int m,
                              const arma::mat& Sigma){

  // Initialize constants
  int n = X.n_rows; // Number of rows
  int p = X.n_cols; // Dimension of data
  int nm = n - m; // Difference between n and m
  int q = p * (p + 1) / 2; // Number of elements in upper triangular matrix of size p x p
  const double s2 = std::sqrt(2); // Store sqrt(2) as constant

  // Matrix A where each column represents the sum of all pairs {i_1, i_2} for a given i_2
  arma::mat A(q, m - 2, arma::fill::zeros);
  for (int i = 0; i < m - 2; i++) {
    // Calculate elements of the matrix crossprod(difference vetor) - 2sigma one at a time
    // Flatten and add matrix as scaled column in A
    // Exploit symmetry by only keeping upper triangular values
    // Off diagonal elements are multiplied by sqrt(2) to double count the symmetric values when multiplying out later
    int l = 0;
    for (int ii = 0; ii < p; ii++) {
      for (int jj = ii; jj < p; jj++) {
        if (ii == jj) {
          A(l++, i) = accu(X(arma::span(0, i), ii) % X(arma::span(0,i), jj)) - X(i + 1, ii) * accu(X(arma::span(0, i), jj)) -
            X(i + 1, jj) * accu(X(arma::span(0, i), ii)) + (i + 1) * X(i + 1, ii) * X(i + 1, jj) - (2 * (i + 1) * Sigma(ii, jj));
        } else {
          A(l++, i) = s2 * (accu(X(arma::span(0, i), ii) % X(arma::span(0,i), jj)) - X(i + 1, ii) * accu(X(arma::span(0, i), jj)) -
            X(i + 1, jj) * accu(X(arma::span(0, i), ii)) + (i + 1) * X(i + 1, ii) * X(i + 1, jj) - (2 * (i + 1) * Sigma(ii, jj)));
        }
      }
    }
  }

  // Compute cumulative stats
  arma::vec results(nm, arma::fill::zeros);

  // Iterate over each slice in the cube B
  // For each "first two" pair, add up all the columns in the slice whose column index come after the second value of the pair
  // Add up these sums as we iterate to calculate the cumulative value
  arma::mat B(q, m - 2, arma::fill::zeros);
  for (int j = 0; j < nm; j++) {

    // Pre-compute slice of cube-object of the final two terms, B_{i3, j} summed over i3 = i ,..., m (backwards cumsum)
    // Do the last column separately so previous columns can be built cumulatively from them
    arma::vec difference_vec = X.row((m - 2) - 1 + 2).t() - X.row(j + m).t();

    // Last column
    int k = 0;
    for (int ii = 0; ii < p; ii++) {
      for (int jj = ii; jj < p; jj++) {
        if (ii == jj) {
          B(k, (m - 2) - 1) = difference_vec(ii) * difference_vec(jj) - (2 * Sigma(ii, jj));
          k++;
        } else {
          B(k, (m - 2) - 1) =  s2 *  (difference_vec(ii) * difference_vec(jj) - (2 * Sigma(ii, jj)));
          k++;
        }
      }
    }

    // All previous columns
    for (int i = (m - 2) - 2; i > -1; i--){
      // Generate difference matrix
      arma::vec difference_vec = X.row(i + 2).t() - X.row(j + m).t();

      // Flatten and add matrix as scaled column in B plus all future columns
      int k = 0;
      for (int ii = 0; ii < p; ii++) {
        for (int jj = ii; jj < p; jj++) {
          if (ii == jj) {
            B(k, i) = B(k, i + 1) + difference_vec(ii) * difference_vec(jj) - (2 * Sigma(ii, jj));
            k++;
          } else {
            B(k, i) = B(k, i + 1) + s2 * (difference_vec(ii) * difference_vec(jj) - (2 * Sigma(ii, jj)));
            k++;
          }
        }
      }
    }
    results(j) = arma::accu(A % B);
  }

  return arma::cumsum(results)/4;
}

// [[Rcpp::export]]
double crossWStatVar_c(const arma::mat& X,
                       int m,
                       const arma::mat& Sigma){

  // Number of rows
  double n = X.n_rows;

  // Calculate row vector of cumulative stats
  arma::vec cumulatives = crossUStatVar_cum_c(X, m, Sigma);

  // Theta_{(m+1) : n} estimate
  double theta = cumulatives(n - m - 1);

  // Normalizer estimate
  double normalizer = accu(pow(cumulatives - theta * arma::regspace(1, n - m)/(n - m), 2)) / (n - m);


  // W test statistic
  return pow(theta, 2) / normalizer;
}
