#' Approximate Cross W Distribution Quantiles
#'
#' @param p Vector of probabilities
#'
#' @returns Approximate quantile of the distribution based on \eqn{1e6} simulations of the random variable defined as
#' \deqn{
#' W:=\frac{{B(1)^2}}{\int_{0}^1(B(r) - rB(1))^2dr}
#' }
#' @export
#'
#' @examples
qcrossW <- function(p){
  # Check input of p
  if( any(p > 1) | any(p < 0)){
    stop("p must be in [0, 1]")
  }
  if(!is.numeric(p)) {
    stop('p must be numeric')
  }

  # Return vectorized list of quantile values
  results <- sapply(p, \(p) {
    if(p == 1){
      return(Inf)
    } else if (p == 0) {
      return(0)
    } else {
    return(quantile(DATASET, p, names = FALSE))
    }}
  )

  return(results)
}

#' Approximate Cross W Distribution Probabilities
#'
#' @param q vector of quantiles
#'
#' @returns Empirical CDF of the distribution based on \eqn{1e6} simulations of the random variable defined as
#' \deqn{
#' W:=\frac{{B(1)^2}}{\int_{0}^1(B(r) - rB(1))^2dr}
#' }
#' @export
#'
#' @examples
pcrossW <- function(q){
  # Check input
  if(!is.numeric(p)) {
    stop('q must be numeric')
  }

  # Return vectorized list of empirical cdf values
  results <- sapply(q, \(q) {
    mean(q > DATASET)
  }
  )
}

#' Title
#'
#' @param n
#' @param p
#'
#' @returns
#' @export
#'
#' @examples
rcrossW <- function(n, p = 200) {
  # Check inputs
  if(!is.numeric(n) || (n %% 1 != 0) ){
    stop('n must be an integer')
  }
  if(!is.numeric(p) || (p %% 1 != 0) ){
    stop('p must be an integer')
  }

  # Add one to adjust for numerator generation
  p <- p + 1

  # Simulated generation
  Z <- matrix(rnorm(n * p), n, p)
  Z0 <- Z[ , 1]
  lambda <- matrix((pi * 1:(p - 1)), n, (p - 1), byrow = T)
  Q <- rowSums((Z[ , -1] / lambda)^2)
  results <- Z0^2 / Q

  return(results)
}

dcrossW <- function(x) {
  # Check input
  if(!is.numeric(x)) {
    stop('x must be numeric')
  }

  dens.func(x)
}


# Basic density function approximation
dens = density(DATASET, n = 1e6)
dens.func = approxfun(dens)
