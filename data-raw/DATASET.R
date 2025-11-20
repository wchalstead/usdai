## code to prepare `DATASET` dataset goes here

# DATASET is n realizations of an approximation of the random variable that our self normalizer converges to
# Calculated using the KL expansion of the Brownian bridge
set.seed(9430857)

# Set parameters
n <- 1e6
p <- 201

# Generate random variables and calculate summations
Z <- matrix(rnorm(n * p), n, p)
Z0 <- Z[ , 1]
lambda <- matrix((pi * 1:(p - 1)), n, (p - 1), byrow = T)
Q <- rowSums((Z[ , -1] / lambda)^2) + 0.00051 # Add approximate expected remainder to reduce bias
DATASET <- Z0^2/Q

# Remove extraneous things from memory
rm(lambda, Z, n, Q, Z0, p)

# Save DATASET to be used internally
usethis::use_data(DATASET, overwrite = TRUE, internal = TRUE)
