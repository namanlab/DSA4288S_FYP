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
  return(tv_distance)
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



# Objective function combining total variation distance, changes in beta coefficients, and R^2 differences across Z categories
objective_function <- function(X_prime, Y_prime, X, Y, Z, beta0_orig, beta1_orig, lambda1, lambda2, lambda3, R2_orig) {
  # Compute total variation distance
  TV_X <- calculate_tv_distance_empirical(X, X_prime)
  TV_Y <- calculate_tv_distance_empirical(Y, Y_prime)
  
  # Fit regression model for the entire dataset
  fit <- lm(Y_prime ~ X_prime)
  beta0_new <- coef(fit)[1]
  beta1_new <- coef(fit)[2]
  
  # Compute the changes in regression coefficients
  delta_beta0 <- beta0_new - beta0_orig
  delta_beta1 <- beta1_new - beta1_orig
  
  # Compute R^2 for each category in Z
  R2_diff <- 0
  categories <- unique(Z)
  for (cat in categories) {
    X_cat <- X[Z == cat]
    Y_cat <- Y[Z == cat]
    X_prime_cat <- X_prime[Z == cat]
    Y_prime_cat <- Y_prime[Z == cat]
    
    # Original R^2 for this category
    fit_orig_cat <- lm(Y_cat ~ X_cat)
    R2_orig_cat <- summary(fit_orig_cat)$r.squared
    
    # New R^2 for this category
    fit_prime_cat <- lm(Y_prime_cat ~ X_prime_cat)
    R2_new_cat <- summary(fit_prime_cat)$r.squared
    
    # Sum of squared differences in R^2 across all categories
    R2_diff <- R2_diff + (R2_new_cat - R2_orig_cat)^2
  }
  
  # Loss function: combine TV distance, coefficient changes, and R^2 differences
  loss_1 <- TV_X + TV_Y + lambda1 * delta_beta0^2 + lambda2 * delta_beta1^2  
  print(loss_1)
  optim_2 <- lambda3 * R2_diff
  print(optim_2)
  return(loss_1 - optim_2)
}

# Simulated Annealing including categorical variable Z and R^2 differences
simulated_annealing <- function(X, Y, Z, lambda1 = 1, lambda2 = 1, lambda3 = 1, max_iter = 1000, initial_temp = 1.0, cooling_rate = 0.99) {
  X_prime <- X
  Y_prime <- Y
  best_X_prime <- X_prime
  best_Y_prime <- Y_prime
  
  # Original regression coefficients
  fit_orig <- lm(Y ~ X)
  beta0_orig <- coef(fit_orig)[1]
  beta1_orig <- coef(fit_orig)[2]
  
  # Original R^2 for each category of Z
  categories <- unique(Z)
  R2_orig <- sapply(categories, function(cat) {
    fit_cat <- lm(Y[Z == cat] ~ X[Z == cat])
    summary(fit_cat)$r.squared
  })
  
  best_loss <- objective_function(X_prime, Y_prime, X, Y, Z, beta0_orig, beta1_orig, lambda1, lambda2, lambda3, R2_orig)
  
  temp <- initial_temp
  for (i in 1:max_iter) {
    # Generate a new solution (random perturbation)
    new_X_prime <- X_prime + rnorm(length(X_prime), 0, 0.01)
    new_Y_prime <- Y_prime + rnorm(length(Y_prime), 0, 0.01)
    
    # Compute the new loss
    new_loss <- objective_function(new_X_prime, new_Y_prime, X, Y, Z, beta0_orig, beta1_orig, lambda1, lambda2, lambda3, R2_orig)
    
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
result <- simulated_annealing(x, y, z)
p1 <- ggplot(data.frame(x = x, xp = result$X_prime)) +
  geom_density(aes(x = x), color = "blue") +
  geom_density(aes(x = xp), color = "red") +
  theme_bw() + labs(main = "x (blue) vs x' (red)")
p2 <- ggplot(data.frame(y = y, yp = result$Y_prime)) +
  geom_density(aes(x = y), color = "blue") +
  geom_density(aes(x = yp), color = "red") +
  theme_bw() + labs(main = "y (blue) vs y' (red)")
p3 <- ggplot(data.frame(x = x, y = y, z = z)) +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() + labs(main = "x vs y")
p4 <- ggplot(data.frame(x = result$X_prime, y = result$Y_prime, z = z)) +
    geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() + labs(main = "x' vs y'")
grid.arrange(p1, p2, p3, p4) 



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




