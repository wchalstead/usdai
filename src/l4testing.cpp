#include <RcppArmadillo.h>

// [[Rcpp::export]]
double crossUStatL4_c(const arma::mat& X,
                      int m){

  // Get number of rows in X
  double n = X.n_rows;


  // Calculate power sums
  arma::rowvec p1 = arma::sum(X.rows(0, m-1), 0);
  arma::rowvec p2 = arma::sum(pow(X.rows(0, m-1), 2), 0);
  arma::rowvec p3 = arma::sum(pow(X.rows(0, m-1), 3), 0);

  // Calculate product of first three terms
  arma::rowvec e3 = (pow(p1, 3) - 3 * p1 % p2 + 2 * p3) / 6;


  // Dot product with sum of last term
  double result = as_scalar(sum(X.rows(m, n - 1), 0) * e3.t());

  return(result);
}


// Cumulative value of l4 norm estimates as j goes from m+1, ..., n
// data is nxp matrix
// m is integer cross "threshold" value
// Function for calculating cumulative L4 cross stats
// Returns R^(n - m) vector of estimates for each j in (m+1, ..., n)
// [[Rcpp::export]]
arma::vec crossUStatL4_cum_c(const arma::mat& X,
                             int m){
  // Get number of rows in X
  double n = X.n_rows;


  // Calculate power sums
  arma::rowvec p1 = arma::sum(X.rows(0, m-1), 0);
  arma::rowvec p2 = arma::sum(pow(X.rows(0, m-1), 2), 0);
  arma::rowvec p3 = arma::sum(pow(X.rows(0, m-1), 3), 0);

  // Calculate product of first three terms
  arma::rowvec e3 = (pow(p1, 3) - 3 * p1 % p2 + 2 * p3) / 6;

  arma::mat dot_vals = X.rows(m, n-1) * e3.t();

  return cumsum(dot_vals);
}


// Calculate self-normalized value using cumulative stats
// [[Rcpp::export]]
double crossWStatL4_c(const arma::mat& X,
                      int m){

  // Number of rows
  double n = X.n_rows;

  // Calculate col vector of cumulative stats
  arma::vec cumulatives = crossUStatL4_cum_c(X, m);

  // Theta_{(m+1) : n} estimate
  double theta = cumulatives(n - m - 1);

  // Normalizer estimate
  double normalizer = accu(pow(cumulatives - theta * arma::regspace(1, n - m)/(n - m), 2)) / (n - m);


  // W test statistic
  return pow(theta, 2) / normalizer;
}
