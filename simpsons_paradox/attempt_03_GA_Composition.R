# Load necessary libraries
library(MCMCpack)
library(stringr)
library(GGally)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(akima)
library(actuar)
library(copula)
library(goftest)
options(scipen = 999)
library(gridExtra)

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


# Objective function combining total variation distance, changes in beta coefficients, and R^2 differences across Z categories
objective_function <- function(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig, 
                               lambda1, lambda2, lambda3, R2_orig, printc = F) {
  # Compute total variation distance
  TV_X <- calculate_tv_distance_empirical(X, X_prime)
  TV_Y <- calculate_tv_distance_empirical(Y, Y_prime)
  
  # Fit regression model for the entire dataset
  fit <- lm(Y_prime ~ X_prime)
  beta0_new <- coef(fit)[1]
  beta1_new <- coef(fit)[2]
  
  # Compute the changes in regression coefficients
  delta_beta0 <- ((beta0_new - beta0_orig)/beta0_orig)^2
  delta_beta1 <- ((beta1_new - beta1_orig)/beta1_orig)^2
  
  # Compute R^2 for each category in Z
  R2_diff <- 0
  categories <- unique(Z)
  for (cat in categories) {
    X_cat <- X[Z == cat]
    Y_cat <- Y[Z == cat]
    X_prime_cat <- X_prime[Z == cat]
    Y_prime_cat <- Y_prime[Z == cat]
    # New R^2 for this category
    rho_curr <- cor(Y_prime_cat, X_prime_cat)
    R2_diff <- R2_diff + (rho_curr - p[cat])^2 # R2_orig_cat is target corr at start
    if (printc){
      cat("rho inter:", cat, "act:", rho_curr, "expected:", p[cat], "sqdiff:", (rho_curr - p[cat])^2, "\n")
    }
  }
  R2_diff <- R2_diff/length(categories)
  cat("TV:", TV_X, TV_Y,"Beta 0:", delta_beta0, "Beta 1:", delta_beta1, "R2:", R2_diff, "\n")
  
  # Loss function: combine TV distance, coefficient changes, and R^2 differences
  loss <- TV_X + TV_Y + lambda1 * delta_beta0^2 + lambda2 * delta_beta1^2  + lambda3 * R2_diff
  return(loss)
}

genetic_algorithm <- function(X, Y, Z, X_st, Y_st, p, 
                              lambda1 = 1, lambda2 = 1, lambda3 = 1, 
                              pop_size = 100, num_generations = 100, 
                              mutation_prob = 0.01, crossover_prob = 0.7) {
  
  # Initialize population
  population <- list()
  for (i in 1:pop_size) {
    X_prime <- X_st + rnorm(length(X_st), 0, 0.05)  # Small random noise
    Y_prime <- Y_st + rnorm(length(Y_st), 0, 0.05)
    population[[i]] <- list(X_prime = X_prime, Y_prime = Y_prime)
  }
  
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
  
  # Fitness function
  fitness_function <- function(individual) {
    X_prime <- individual$X_prime
    Y_prime <- individual$Y_prime
    return(objective_function(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig, lambda1, lambda2, lambda3, R2_orig))
  }
  
  # Selection (Tournament Selection)
  select_parents <- function(population, fitness_scores, k = 5) {
    selected <- sample(1:length(population), k, replace = FALSE)
    parent_idx <- selected[which.min(fitness_scores[selected])]
    return(population[[parent_idx]])
  }
  
  # Crossover (Single Point Crossover)
  crossover <- function(parent1, parent2) {
    if (runif(1) < crossover_prob) {
      crossover_point <- sample(1:length(parent1$X_prime), 1)
      child1_X <- c(parent1$X_prime[1:crossover_point], parent2$X_prime[(crossover_point + 1):length(parent2$X_prime)])
      child1_Y <- c(parent1$Y_prime[1:crossover_point], parent2$Y_prime[(crossover_point + 1):length(parent2$Y_prime)])
      child2_X <- c(parent2$X_prime[1:crossover_point], parent1$X_prime[(crossover_point + 1):length(parent1$X_prime)])
      child2_Y <- c(parent2$Y_prime[1:crossover_point], parent1$Y_prime[(crossover_point + 1):length(parent1$Y_prime)])
      
      return(list(list(X_prime = child1_X, Y_prime = child1_Y), list(X_prime = child2_X, Y_prime = child2_Y)))
    } else {
      return(list(parent1, parent2))
    }
  }
  
  # Mutation
  mutate <- function(individual, mutation_prob) {
    if (runif(1) < mutation_prob) {
      individual$X_prime <- individual$X_prime + rnorm(length(individual$X_prime), 0, 0.1)
      individual$Y_prime <- individual$Y_prime + rnorm(length(individual$Y_prime), 0, 0.1)
    }
    return(individual)
  }
  
  mutate_with_copula <- function(individual, Z, p_corrs, F1Inv, F2Inv, mutation_prob) {
    unique_cats <- unique(Z)
    
    for (cat in unique_cats) {
      # Get samples for this category
      x_samples <- individual$X_prime[Z == cat]
      y_samples <- individual$Y_prime[Z == cat]
      
      # Check if mutation should be applied
      if (runif(1) < mutation_prob) {
        # Generate copula samples with the desired correlation
        n_samples <- sum(Z == cat)  # Number of samples for current category
        p_corr <- p_corrs[cat]  # Desired correlation for the current category
        
        # Generate copula samples
        copula_samples <- gauss_copula_2(n_samples, p_corr)
        
        # Transform copula samples using the inverse CDF functions
        x_transformed <- F1Inv(copula_samples[, 1])
        y_transformed <- F2Inv(copula_samples[, 2])
        
        # Replace original samples with transformed ones
        individual$X_prime[Z == cat] <- x_transformed
        individual$Y_prime[Z == cat] <- y_transformed
      }
    }
    
    return(individual)
  }
  
  # Genetic Algorithm
  for (gen in 1:num_generations) {
    # Calculate fitness scores
    fitness_scores <- sapply(population, fitness_function)
    
    # Create a new population using selection, crossover, and mutation
    new_population <- list()
    for (i in seq(1, pop_size, by = 2)) {
      # Selection
      parent1 <- select_parents(population, fitness_scores)
      parent2 <- select_parents(population, fitness_scores)
      
      # Crossover
      offspring <- crossover(parent1, parent2)
      
      # Mutation
      offspring[[1]] <- mutate_with_copula(offspring[[1]], mutation_prob)
      offspring[[2]] <- mutate_with_copula(offspring[[2]], mutation_prob)
      
      # Add offspring to the new population
      new_population[[i]] <- offspring[[1]]
      new_population[[i + 1]] <- offspring[[2]]
    }
    
    # Replace old population with the new one
    population <- new_population
  }
  
  # Select the best solution from the final population
  final_fitness_scores <- sapply(population, fitness_function)
  best_individual_idx <- which.min(final_fitness_scores)
  best_solution <- population[[best_individual_idx]]
  print(objective_function(best_solution$X_prime, best_solution$Y_prime, X, Y, Z, p, 
                           beta0_orig, beta1_orig, lambda1, lambda2, lambda3, R2_orig,
                           printc = T))
  return(best_solution)
}

set.seed(123)
n <- 1000
z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 5) + 5*rbeta(n, 5, 3)
y <- 2*x + rnorm(n, 5, sd = 4)
t = c(-0.8, -0.8, -0.8)
p = sin(t*pi/2) 
names(p) <- unique(z)

# Run the genetic algorithm
best_solution <- genetic_algorithm(X = x, Y = y, Z = z, X_st = x, Y_st = y, p = p,
                                   lambda3 = 10)

p1 <- ggplot(data.frame(x = x, xp = best_solution$X_prime[1:1000])) +
  geom_density(aes(x = x), color = "blue") +
  geom_density(aes(x = xp), color = "red") +
  theme_bw() + labs(main = "x (blue) vs x' (red)")
p2 <- ggplot(data.frame(y = y, yp = best_solution$Y_prime[1:1000])) +
  geom_density(aes(x = y), color = "blue") +
  geom_density(aes(x = yp), color = "red") +
  theme_bw() + labs(main = "y (blue) vs y' (red)")
p3 <- ggplot(data.frame(x = x, y = y, z = z)) +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() + labs(main = "x vs y")
p4 <- ggplot(data.frame(x = best_solution$X_prime[1:1000], y = best_solution$Y_prime[1:1000], z = z)) +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() + labs(main = "x' vs y'")
grid.arrange(p1, p2, p3, p4) 








################################################################################
############################# NEW RULES: COPULA ############################
################################################################################






genetic_algorithm_corr <- function(X, Y, Z, X_st, Y_st, p, 
                                   lambda1 = 1, lambda2 = 1, lambda3 = 1, 
                                   pop_size = 100, num_generations = 100, 
                                   mutation_prob = 0.1, crossover_prob = 0.7,
                                   sel_parents_fn = select_parents_roulette,
                                   crossover_fn = crossover_trad,
                                   mutate_fn = mutate_corr) {
  
  # Initialize population
  population <- list()
  for (i in 1:pop_size) {
    X_prime <- X_st + rnorm(length(X_st), 0, 0.1) 
    Y_prime <- Y_st + rnorm(length(Y_st), 0, 0.1)
    population[[i]] <- list(X_prime = X_prime, Y_prime = Y_prime)
  }
  
  # Original regression coefficients
  fit_orig <- lm(Y ~ X)
  beta0_orig <- coef(fit_orig)[1]
  beta1_orig <- coef(fit_orig)[2]
  
  orig_corr = cor(X, Y)
  
  # Generate inverse CDF functions for x and y samples
  F1Inv <- Vectorize(genCDFInv_linear(X))
  F2Inv <- Vectorize(genCDFInv_linear(Y))
  
  # Fitness function modified to account for category-specific correlations
  fitness_function <- function(individual) {
    X_prime <- individual$X_prime
    Y_prime <- individual$Y_prime
    return(objective_function(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig, lambda1, lambda2, lambda3, R2_orig))
  }
  
  
  # Main Genetic Algorithm loop
  for (gen in 1:num_generations) {
    fitness_scores <- sapply(population, fitness_function)
    
    new_population <- list()
    for (i in seq(1, pop_size, by = 2)) {
      parent1 <- sel_parents_fn(population, fitness_scores)
      parent2 <- sel_parents_fn(population, fitness_scores)
      
      if (identical(crossover_fn, crossover_corr_adjust)){
        offspring <- crossover_corr_adjust(parent1, parent2, orig_corr, crossover_prob)
      } else {
        offspring <- crossover_fn(parent1, parent2, Z, crossover_prob)
      }
      
      if (identical(mutate_fn, mutate_with_copula)){
        offspring[[1]] <- mutate_with_copula(offspring[[1]], Z, p, F1Inv, F2Inv, mutation_prob)
        offspring[[2]] <- mutate_with_copula(offspring[[2]], Z, p, F1Inv, F2Inv, mutation_prob)
      } else if (identical(mutate_fn, mutate_corr_adjust)){
        offspring[[1]] <- mutate_corr_adjust(offspring[[1]], orig_corr, mutation_prob)
        offspring[[2]] <- mutate_corr_adjust(offspring[[2]], orig_corr, mutation_prob)
      } else {
        offspring[[1]] <- mutate_fn(offspring[[1]], Z, mutation_prob)
        offspring[[2]] <- mutate_fn(offspring[[2]], Z, mutation_prob)
      }
      
      new_population[[i]] <- offspring[[1]]
      new_population[[i + 1]] <- offspring[[2]]
    }
    
    population <- new_population
  }
  
  final_fitness_scores <- sapply(population, fitness_function)
  best_individual_idx <- which.min(final_fitness_scores)
  best_solution <- population[[best_individual_idx]]
  print(objective_function(best_solution$X_prime, best_solution$Y_prime, X, Y, Z, p, 
                           beta0_orig, beta1_orig, lambda1, lambda2, lambda3, R2_orig,
                           printc = T))
  return(best_solution)
}



set.seed(123)
n <- 1000
z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 5) + 5*rbeta(n, 5, 3)
y <- 2*x + rnorm(n, 5, sd = 4)
t = c(-0.8, 0.8, -0.8)
p = sin(t*pi/2) 
names(p) <- unique(z)

# Run the genetic algorithm
best_solution <- genetic_algorithm_corr(X = x, Y = y, Z = z, X_st = x, Y_st = y, p = p,
                                   mutation_prob = 0.1, lambda3 = 10,
                                   num_generations = 100,
                                   mutate_fn = mutate_with_copula)

p1 <- ggplot(data.frame(x = x, xp = best_solution$X_prime[1:1000])) +
  geom_density(aes(x = x), color = "blue") +
  geom_density(aes(x = xp), color = "red") +
  theme_bw() + labs(main = "x (blue) vs x' (red)")
p2 <- ggplot(data.frame(y = y, yp = best_solution$Y_prime[1:1000])) +
  geom_density(aes(x = y), color = "blue") +
  geom_density(aes(x = yp), color = "red") +
  theme_bw() + labs(main = "y (blue) vs y' (red)")
p3 <- ggplot(data.frame(x = x, y = y, z = z)) +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() + labs(main = "x vs y")
p4 <- ggplot(data.frame(x = best_solution$X_prime[1:1000], y = best_solution$Y_prime[1:1000], z = z)) +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() + labs(main = "x' vs y'")
grid.arrange(p1, p2, p3, p4) 








################################################################################
################################### COMBINE ####################################
################################################################################



objective_function <- function(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig, 
                               lambda1, lambda2, lambda3, lambda4, R2_orig, printc = F) {
  # Compute total variation distance
  TV_X <- calculate_tv_distance_empirical(X, X_prime)
  TV_Y <- calculate_tv_distance_empirical(Y, Y_prime)
  
  # Fit regression model for the entire dataset
  fit <- lm(Y_prime ~ X_prime)
  beta0_new <- coef(fit)[1]
  beta1_new <- coef(fit)[2]
  
  # Compute the changes in regression coefficients
  delta_beta0 <- ((beta0_new - beta0_orig)/beta0_orig)^2
  delta_beta1 <- ((beta1_new - beta1_orig)/beta1_orig)^2
  
  # Compute R^2 for each category in Z
  R2_diff <- 0
  categories <- unique(Z)
  for (cat in categories) {
    X_cat <- X[Z == cat]
    Y_cat <- Y[Z == cat]
    X_prime_cat <- X_prime[Z == cat]
    Y_prime_cat <- Y_prime[Z == cat]
    
    # New R^2 for this category
    rho_curr <- cor(Y_prime_cat, X_prime_cat)
    R2_diff <- R2_diff + (rho_curr - p[cat])^2 # R2_orig_cat is target corr at start
    if (printc) {
      cat("rho inter:", cat, "act:", rho_curr, "expected:", p[cat], "sqdiff:", (rho_curr - p[cat])^2, "\n")
    }
  }
  R2_diff <- R2_diff / length(categories)
  
  # Compute inter-cluster (centroid) distances
  cluster_centroids <- lapply(categories, function(cat) {
    X_cat_prime <- X_prime[Z == cat]
    Y_cat_prime <- Y_prime[Z == cat]
    return(c(mean(X_cat_prime), mean(Y_cat_prime))) # Centroid for (X, Y) in category 'cat'
  })
  
  inter_cluster_dist <- 0
  num_pairs <- 0
  for (i in 1:(length(categories) - 1)) {
    for (j in (i + 1):length(categories)) {
      centroid_i <- cluster_centroids[[i]]
      centroid_j <- cluster_centroids[[j]]
      
      # Compute Euclidean distance between centroids
      dist_ij <- sqrt(sum((centroid_i - centroid_j)^2))
      inter_cluster_dist <- inter_cluster_dist + dist_ij  
      num_pairs <- num_pairs + 1
    }
  }
  
  # Average inverse of inter-cluster distances
  inter_cluster_dist <- inter_cluster_dist / num_pairs
  
  # Use quantiles to compute a robust min and max distance (to handle outliers)
  quantile_low <- 0.05  # 5th percentile for minimum
  quantile_high <- 0.95 # 95th percentile for maximum
  
  # Quantile-based min and max distances for normalization
  range_X <- quantile(X_prime, probs = c(quantile_low, quantile_high))
  range_Y <- quantile(Y_prime, probs = c(quantile_low, quantile_high))
  
  # Min and max based on quantile ranges
  min_dist <- 0 # Since min_dist is still 0 (clusters perfectly overlap)
  max_dist <- sqrt((range_X[2] - range_X[1])^2 + (range_Y[2] - range_Y[1])^2) # Max possible distance based on quantiles
  
  # Perform quantile-based min-max normalization
  normalized_inter_cluster_dist <- (inter_cluster_dist - min_dist) / (max_dist - min_dist)
  
  # Ensure the normalized value is between 0 and 1
  inter_cluster_dist <- max(0, min(normalized_inter_cluster_dist, 1))
  
  # Display intermediate metrics
  cat("TV:", TV_X, TV_Y, "Beta 0:", delta_beta0, "Beta 1:", delta_beta1, "R2:", R2_diff, "Cluster Dist:", inter_cluster_dist, "\n")
  
  # Loss function: combine TV distance, coefficient changes, R^2 differences, and cluster separation penalty
  loss <- TV_X + TV_Y + lambda1 * delta_beta0^2 + lambda2 * delta_beta1^2 + lambda3 * R2_diff + lambda4 * inter_cluster_dist
  return(loss)
}

genetic_algorithm_corr <- function(X, Y, Z, X_st, Y_st, p, 
                                   lambda1 = 1, lambda2 = 1, lambda3 = 1, lambda4 = 1,
                                   pop_size = 100, num_generations = 100, 
                                   mutation_prob = 0.1, crossover_prob = 0.7,
                                   sel_parents_fn = select_parents_roulette,
                                   crossover_fn = crossover_trad,
                                   mutate_fn = mutate_corr) {
  
  # Initialize population
  population <- list()
  for (i in 1:pop_size) {
    X_prime <- X_st + rnorm(length(X_st), 0, 0.1) 
    Y_prime <- Y_st + rnorm(length(Y_st), 0, 0.1)
    population[[i]] <- list(X_prime = X_prime, Y_prime = Y_prime)
  }
  
  # Original regression coefficients
  fit_orig <- lm(Y ~ X)
  beta0_orig <- coef(fit_orig)[1]
  beta1_orig <- coef(fit_orig)[2]
  
  orig_corr = cor(X, Y)
  
  # Generate inverse CDF functions for x and y samples
  F1Inv <- Vectorize(genCDFInv_linear(X))
  F2Inv <- Vectorize(genCDFInv_linear(Y))
  
  # Fitness function modified to account for category-specific correlations
  fitness_function <- function(individual) {
    X_prime <- individual$X_prime
    Y_prime <- individual$Y_prime
    return(objective_function(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig, lambda1, lambda2, lambda3, lambda4, R2_orig))
  }
  
  
  # Main Genetic Algorithm loop
  for (gen in 1:num_generations) {
    fitness_scores <- sapply(population, fitness_function)
    
    new_population <- list()
    for (i in seq(1, pop_size, by = 2)) {
      parent1 <- sel_parents_fn(population, fitness_scores)
      parent2 <- sel_parents_fn(population, fitness_scores)
      
      if (identical(crossover_fn, crossover_corr_adjust)){
        offspring <- crossover_corr_adjust(parent1, parent2, orig_corr, crossover_prob)
      } else {
        offspring <- crossover_fn(parent1, parent2, Z, crossover_prob)
      }
      
      if (identical(mutate_fn, mutate_with_copula)){
        offspring[[1]] <- mutate_with_copula(offspring[[1]], Z, p, F1Inv, F2Inv, mutation_prob)
        offspring[[2]] <- mutate_with_copula(offspring[[2]], Z, p, F1Inv, F2Inv, mutation_prob)
      } else if (identical(mutate_fn, mutate_corr_adjust)){
        offspring[[1]] <- mutate_corr_adjust(offspring[[1]], orig_corr, mutation_prob)
        offspring[[2]] <- mutate_corr_adjust(offspring[[2]], orig_corr, mutation_prob)
      } else {
        offspring[[1]] <- mutate_fn(offspring[[1]], Z, mutation_prob)
        offspring[[2]] <- mutate_fn(offspring[[2]], Z, mutation_prob)
      }
      
      new_population[[i]] <- offspring[[1]]
      new_population[[i + 1]] <- offspring[[2]]
    }
    
    population <- new_population
  }
  
  final_fitness_scores <- sapply(population, fitness_function)
  best_individual_idx <- which.min(final_fitness_scores)
  best_solution <- population[[best_individual_idx]]
  print(objective_function(best_solution$X_prime, best_solution$Y_prime, X, Y, Z, p, 
                           beta0_orig, beta1_orig, lambda1, lambda2, lambda3, lambda4, R2_orig,
                           printc = T))
  return(best_solution)
}


modify_data <- function(x, y, z, 
                        corr_vector,   # Vector of correlations
                        FInvFunc = genCDFInv_linear,      # Inverse CDF function to transform samples
                        lambda1 = 1, lambda2 = 1, lambda3 = 1, lambda4 = 1, 
                        pop_size = 100, num_generations = 100, 
                        mutation_prob = 0.1, crossover_prob = 0.7,
                        sel_parents_fn = select_parents_roulette,
                        crossover_fn = crossover_trad,
                        mutate_fn = mutate_cor) {     
  
  cat("Starting transformation process...\n")
  
  # Transform correlation vector for Gaussian copula
  p_corr_vec <- sin(corr_vector * pi / 2) 
  
  # Get unique levels of z and sort data by z
  z_levels <- sort(unique(z))
  names(p_corr_vec) <- z_levels
  z <- sort(z)
  
  # Calculate quantiles for x and y based on z distribution
  x_quantiles <- quantile(x, probs = cumsum(table(z) / length(z)))
  y_quantiles <- quantile(y, probs = cumsum(table(z) / length(z)))
  
  df_composition <- data.frame()  # Initialize result dataframe
  
  # Iterate over each category of z
  for (i in seq_along(z_levels)) {
    z_cat <- z_levels[i]
    prob_val <- table(z)[i] / length(z)   # Proportion of category in z
    
    cat("Processing category:", z_cat, "with proportion:", prob_val, "\n")
    
    # Generate samples using the Gaussian copula
    p_corr <- p_corr_vec[i]
    n_samples <- sum(z == z_cat)  # Number of samples for current category
    copula_samples <- gauss_copula_2(n_samples, p_corr)  # Generate copula samples
    
    # Filter x and y samples for the current quantile range
    x_range <- if (i == 1) {
      which(x <= x_quantiles[i])
    } else {
      which(x > x_quantiles[i - 1] & x <= x_quantiles[i])
    }
    y_range <- if (i == 1) {
      which(y <= y_quantiles[i])
    } else {
      which(y > y_quantiles[i - 1] & y <= y_quantiles[i])
    }
    
    x_samples <- x[x_range]
    y_samples <- y[y_range]
    
    # Generate inverse CDF functions for x and y samples
    F1Inv <- Vectorize(FInvFunc(x_samples))
    F2Inv <- Vectorize(FInvFunc(y_samples))
    
    # Transform copula samples using the inverse CDF functions
    x_transformed <- F1Inv(copula_samples[, 1])
    y_transformed <- F2Inv(copula_samples[, 2])
    
    # Append transformed data to the result dataframe
    df_composition <- rbind(df_composition, data.frame(x = x_transformed, y = y_transformed, z = z_cat))
  }
  
  cat("Piecewise copula transformation complete.\n")
  
  # Apply simulated annealing to further optimize the samples
  cat("Starting Genetic Algorithm optimization...\n")
  res_gen <- genetic_algorithm_corr(X = x, Y = y, Z = z, 
                                    X_st = df_composition$x, Y_st = df_composition$y, 
                                    p = p_corr_vec, 
                                    lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, 
                                    lambda4 = lambda4,
                                    pop_size = pop_size, num_generations = num_generations, 
                                    mutation_prob = mutation_prob, crossover_prob = crossover_prob,
                                    sel_parents_fn = sel_parents_fn,
                                    crossover_fn = crossover_fn,
                                    mutate_fn = mutate_fn)
  
  # Create a final dataframe with the optimized x and y values
  df_gen <- data.frame(x_optimized = res_gen$X_prime, 
                             y_optimized = res_gen$Y_prime, 
                             z = z)
  
  cat("GA complete.\n")
  
  # Combine original, copula-transformed, and annealed data into a single dataframe
  df_final <- data.frame(x_original = x, y_original = y, z = z, 
                         x_transformed = df_composition$x, 
                         y_transformed = df_composition$y, 
                         x_optimized = df_gen$x_optimized, 
                         y_optimized = df_gen$y_optimized)
  
  # Plot original, copula-transformed, and simulated annealing results
  p_orig <- ggplot(data.frame(x = x, y = y, z = z)) +
    geom_point(aes(x = x, y = y, color = z)) +
    theme_bw() + guides(color = "none") +
    labs(title = "Original Data")
  p1 <- ggMarginal(p_orig, type = "density")
  
  p_copula <- ggplot(df_composition) +
    geom_point(aes(x = x, y = y, color = z)) +
    theme_bw() + guides(color = "none") +
    labs(title = "After Piecewise Copulas")
  p2 <- ggMarginal(p_copula, type = "density")
  
  p_gen <- ggplot(df_gen) +
    geom_point(aes(x = x_optimized, y = y_optimized, color = z)) +
    theme_bw() + guides(color = "none") +
    labs(title = "After Genetic Algorithm")
  p3 <- ggMarginal(p_gen, type = "density")
  
  # Display all plots side by side
  print(grid.arrange(p1, p2, p3, nrow = 1))
  
  # Calculate total variation (TV) distance before and after annealing
  tv_initial <- calculate_tv_distance_empirical(x, df_composition$x) +
    calculate_tv_distance_empirical(y, df_composition$y)
  tv_final <- calculate_tv_distance_empirical(x, df_gen$x_optimized) +
    calculate_tv_distance_empirical(y, df_gen$y_optimized)
  
  # Print TV distances
  cat("Total Variation Before GA: ", tv_initial, "\n")
  cat("Total Variation After GA: ", tv_final, "\n")
  
  # Return the combined dataframe containing all the transformations
  return(df_final)
}


# Generate data for two categories of z
set.seed(123)
n <- 100
z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 5) + 5*rbeta(n, 5, 3)
y <- 2*x + rnorm(n, 5, sd = 4)
t = c(0.9, 0.7, 0.5) 
res = modify_data(x, y, z, t, sel_parents_fn = select_parents_tournament,
                  mutation_prob = 0.5, 
                  crossover_fn = crossover_trad_z,
                  mutate_fn = mutate_corr_adjust, lambda4 = 10,
                  num_generations = 100)
cor(res$y_original, res$x_original)
cor(res$y_transformed, res$x_transformed)
cor(res$y_optimized, res$x_optimized)


