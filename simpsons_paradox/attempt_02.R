# Load necessary libraries
library(copula)
library(ggplot2)
library(gridExtra)
library(ggExtra)

generate_gaussian_copula_samples <- function(n, d, rho_matrix) {
  # Step 1: Generate multivariate normal samples
  mean_vector <- rep(0, d)  # Mean vector for multivariate normal
  mvn_samples <- mvrnorm(n = n, mu = mean_vector, Sigma = rho_matrix)
  
  # Step 2: Transform samples to uniform using the CDF of the standard normal
  uniform_samples <- pnorm(mvn_samples)
  
  return(uniform_samples)
}
gauss_copula_2 <- function(n, p) {
  rho_matrix <- matrix(c(1, p, p, 1), nrow = 2, ncol = 2)
  smpls <- generate_gaussian_copula_samples(n, 2, rho_matrix)
  return(smpls)
}

genCDFInv_linear <- function(X) {
  ecdf1 <- ecdf(X)
  U <- sort((ecdf1(X)))
  Finv <- sort(X)
  ak2 <- approxfun(U, Finv, method = "linear", rule = 2)
  return(ak2)
}

calculate_tv_distance_empirical <- function(original_data, generated_data) {
  # Create empirical CDFs for original and generated data
  original_ecdf <- ecdf(original_data)
  generated_ecdf <- ecdf(generated_data)
  # Define the grid over which to calculate the distance
  x_values <- sort(unique(c(original_data, generated_data)))
  # Calculate the TV distance
  tv_distance <- max(abs(original_ecdf(x_values) - generated_ecdf(x_values)))
  list(tv_distance = tv_distance)
}

# Step 1: Generate Data
set.seed(123)
n <- 1000

# Generate data for two categories of z
z <- sample(c("A", "B"), prob = c(0.3, 0.7), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 1)
y <- -2*x + rnorm(n, 5, sd = 1)

F1Inv <- Vectorize(genCDFInv_linear(x))
F2Inv <- Vectorize(genCDFInv_linear(y))


t1 = -0.8
t2 = 0.4
p1 = sin(t1*pi/2) 
p2 = sin(t2*pi/2)

# Objective function combining total variation distance and changes in beta coefficients
objective_function <- function(X_prime, Y_prime, X, Y, beta0_orig, beta1_orig, lambda1, lambda2) {
  # Compute total variation distance
  TV_X <- calculate_tv_distance_empirical(X, X_prime)
  TV_Y <- calculate_tv_distance_empirical(Y, Y_prime)
  
  # Fit regression model
  fit <- lm(Y_prime ~ X_prime)
  beta0_new <- coef(fit)[1]
  beta1_new <- coef(fit)[2]
  
  # Compute the changes in regression coefficients
  delta_beta0 <- beta0_new - beta0_orig
  delta_beta1 <- beta1_new - beta1_orig
  
  # Loss function
  loss <- TV_X + TV_Y + lambda1 * delta_beta0^2 + lambda2 * delta_beta1^2
  return(loss)
}


# Simulated Annealing
simulated_annealing <- function(X, Y, lambda1 = 1, lambda2 = 1, max_iter = 1000, initial_temp = 1.0, cooling_rate = 0.99) {
  X_prime <- X
  Y_prime <- Y
  best_X_prime <- X_prime
  best_Y_prime <- Y_prime
  
  # Original regression coefficients
  fit_orig <- lm(Y ~ X)
  beta0_orig <- coef(fit_orig)[1]
  beta1_orig <- coef(fit_orig)[2]
  
  best_loss <- objective_function(X_prime, Y_prime, X, Y, beta0_orig, beta1_orig, lambda1, lambda2)
  
  temp <- initial_temp
  for (i in 1:max_iter) {
    # Generate a new solution (random perturbation)
    new_X_prime <- X_prime + rnorm(length(X_prime), 0, 0.01)
    new_Y_prime <- Y_prime + rnorm(length(Y_prime), 0, 0.01)
    
    # Compute the new loss
    new_loss <- objective_function(new_X_prime, new_Y_prime, X, Y, beta0_orig, beta1_orig, lambda1, lambda2)
    
    # Accept new solution based on probability
    if (new_loss < best_loss || runif(1) < exp((best_loss - new_loss) / temp)) {
      X_prime <- new_X_prime
      Y_prime <- new_Y_prime
      best_loss <- new_loss
      best_X_prime <- X_prime
      best_Y_prime <- Y_prime
    }
    
    # Cool down
    temp <- temp * cooling_rate
    
    # Convergence check
    if (temp < 1e-6) {
      break
    }
  }
  return(list(X_prime = best_X_prime, Y_prime = best_Y_prime))
}

# Example usage
result <- simulated_annealing(x, y)



# Alternating Minimization
alternating_minimization <- function(X, Y, lambda1 = 1, lambda2 = 1, max_iter = 100) {
  X_prime <- X
  Y_prime <- Y
  
  # Original regression coefficients
  fit_orig <- lm(Y ~ X)
  beta0_orig <- coef(fit_orig)[1]
  beta1_orig <- coef(fit_orig)[2]
  
  for (i in 1:max_iter) {
    # Step 1: Minimize total variation distance
    TV_X <- total_variation(X_prime)
    TV_Y <- total_variation(Y_prime)
    
    # Slightly perturb X_prime and Y_prime to minimize TV distance
    X_prime <- X_prime + rnorm(length(X_prime), 0, 0.01)
    Y_prime <- Y_prime + rnorm(length(Y_prime), 0, 0.01)
    
    # Step 2: Adjust regression coefficients
    fit <- lm(Y_prime ~ X_prime)
    beta0_new <- coef(fit)[1]
    beta1_new <- coef(fit)[2]
    
    delta_beta0 <- beta0_new - beta0_orig
    delta_beta1 <- beta1_new - beta1_orig
    
    # Update X_prime and Y_prime based on beta adjustments
    X_prime <- X_prime - lambda1 * delta_beta0
    Y_prime <- Y_prime - lambda2 * delta_beta1
    
    # Convergence check (you can define your own criteria)
    if (i == max_iter) {
      break
    }
  }
  return(list(X_prime = X_prime, Y_prime = Y_prime))
}

# Example usage
result <- alternating_minimization(X, Y)




