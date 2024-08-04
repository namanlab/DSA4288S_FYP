# Load necessary library
library(MCMCpack)
library(GGally)
library(tidyverse)
library(MASS)
library(mvtnorm)

F1 <- function(x1){pnorm(x1, mean = 2, sd = 3)}
F2 <- function(x2){pexp(x2, rate = 0.4)}
F1Inv <- function(p1) { qnorm(p1, mean = 2, sd = 3) }
F2Inv <- function(p2) { qexp(p2, rate = 0.4) }


# Function to generate i.i.d. samples from a Gaussian copula
generate_gaussian_copula_samples <- function(n, d, rho_matrix) {
  # Step 1: Generate multivariate normal samples
  mean_vector <- rep(0, d)  # Mean vector for multivariate normal
  mvn_samples <- mvrnorm(n = n, mu = mean_vector, Sigma = rho_matrix)
  print(dim(mvn_samples))
  # Step 2: Transform samples to uniform using the CDF of the standard normal
  uniform_samples <- pnorm(mvn_samples)
  return(uniform_samples)
}
gauss_copula_2 <- function(n, p){
  rho_matrix <- matrix(c(
    1, p,
    p, 1
  ), nrow = 2, ncol = 2)
  smpls = generate_gaussian_copula_samples(n, 2, rho_matrix)
  print(ggpairs(as.data.frame(smpls)))
  return(smpls)
}

n <- 10000 
res = gauss_copula_2(n, 0.8)
x1 = F1Inv(res[,1])
x2 = F2Inv(res[,2])
ggpairs(data.frame(x1 = x1, x2 = x2))


# Function to generate i.i.d. samples from a t-copula
generate_t_copula_samples <- function(n, d, rho_matrix, df) {
  # Step 1: Generate multivariate t samples
  mvn_samples <- rmvt(n = n, sigma = rho_matrix, df = df)
  # Step 2: Transform samples to uniform using the CDF of the t-distribution
  uniform_samples <- pt(mvn_samples, df = df)
  return(uniform_samples)
}
t_copula_2 <- function(n, p){
  rho_matrix <- matrix(c(
    1, p,
    p, 1
  ), nrow = 2, ncol = 2)
  smpls = generate_t_copula_samples(n, 2, rho_matrix, 5)
  print(ggpairs(as.data.frame(smpls)))
  return(smpls)
}

n <- 10000 
res = t_copula_2(n, 0.8)
x1 = F1Inv(res[,1])
x2 = F2Inv(res[,2])
ggpairs(data.frame(x1 = x1, x2 = x2))


