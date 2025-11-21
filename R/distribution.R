qcrossW <- function(p){
  # Check input of p
  if( any(p > 1) | any(p < 0)){
    stop("p must be in [0, 1]")
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
