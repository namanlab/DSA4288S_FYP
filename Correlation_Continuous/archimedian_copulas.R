library(MCMCpack)
library(GGally)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(akima)
library(actuar)

set.seed(123)
n <- 10000  # Number of samples
d <- 2     # Dimension (number of variables)


################################################################################
################################ CLAYTON ##################################
################################################################################
generate_clayton_copula_samples <- function(n, d, theta) {
  # Step 1: Sample V ~ F = LS^{-1}(ψ), which corresponds to V ~ Gamma(1/theta, 1) 
  V <- rgamma(n, shape = 1 / theta, rate = 1)
  
  # Step 2: Sample Xi ~ U[0,1] (independent uniform samples)
  Xi <- matrix(runif(n * d), nrow = n, ncol = d)
  
  # Step 3: Calculate the U_i using the transformation U_i = (1 - log(Xi) / V) ^ (-1 / theta)
  U <- (1 - log(Xi) / V) ^ (-1 / theta)
  
  return(U)
}
theta_values <- seq(0.1, 5, by = 0.01)  # Range of theta values
correlations <- sapply(theta_values, function(theta) {
  samples <- generate_clayton_copula_samples(n, d, theta)
  cor(samples[, 1], samples[, 2])  # Calculate the correlation between the two dimensions
})
# Create a data frame for plotting
correlation_df <- data.frame(
  Theta = theta_values,
  Correlation = correlations
)
# Plot the correlation against theta
ggplot(correlation_df, aes(x = Theta, y = Correlation)) +
  geom_line(color = "blue") 
labs(title = "Correlation vs Theta in Clayton Copula",
     x = "Theta",
     y = "Correlation") +
  theme_minimal()

################################################################################
################################ AMH ##################################
################################################################################
generate_amh_copula_samples <- function(n, d, theta) {
  # Step 1: Sample V 
  V <- rgeom(n, 1 - theta)
  
  # Step 2: Sample Xi ~ U[0,1] (independent uniform samples)
  Xi <- matrix(runif(n * d), nrow = n, ncol = d)
  
  # Step 3: Calculate the U_i using the transformation U_i = exp(-(-log(Xi) / V)^(1/theta))
  U <- (1 - theta)/(exp(-log(Xi)/V) - theta)
  
  return(U)
}
theta_values <- seq(0.01, 0.99, by = 0.01)  # Range of theta values
correlations <- sapply(theta_values, function(theta) {
  samples <- generate_amh_copula_samples(n, d, theta)
  cor(samples[, 1], samples[, 2])  # Calculate the correlation between the two dimensions
})
# Create a data frame for plotting
correlation_df <- data.frame(
  Theta = theta_values,
  Correlation = correlations
)
# Plot the correlation against theta
ggplot(correlation_df, aes(x = Theta, y = Correlation)) +
  geom_line(color = "blue") 
labs(title = "Correlation vs Theta in Clayton Copula",
     x = "Theta",
     y = "Correlation") +
  theme_minimal()


################################################################################
################################ FRANK ##################################
################################################################################
generate_frank_copula_samples <- function(n, d, theta) {
  # Step 1: Sample V 
  V <- rlogarithmic(n, 1 - exp(-theta))
  
  # Step 2: Sample Xi ~ U[0,1] (independent uniform samples)
  Xi <- matrix(runif(n * d), nrow = n, ncol = d)
  
  # Step 3: Calculate the U_i using the transformation U_i = exp(-(-log(Xi) / V)^(1/theta))
  t <- -log(Xi/V)
  U <- (-log(exp(-t)*(exp(-theta) - 1) + 1))/theta
  
  return(U)
}
theta_values <- seq(0.1, 10, by = 0.01) # Range of theta values
correlations <- sapply(theta_values, function(theta) {
  samples <- generate_frank_copula_samples(n, d, theta)
  cor(samples[, 1], samples[, 2])  # Calculate the correlation between the two dimensions
})
# Create a data frame for plotting
correlation_df <- data.frame(
  Theta = theta_values,
  Correlation = correlations
)
# Plot the correlation against theta
ggplot(correlation_df, aes(x = Theta, y = Correlation)) +
  geom_line(color = "blue") 
labs(title = "Correlation vs Theta in Clayton Copula",
     x = "Theta",
     y = "Correlation") +
  theme_minimal()
