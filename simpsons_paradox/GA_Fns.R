
################################################################################

# Selection (Tournament Selection)
select_parents_tournament <- function(population, fitness_scores, k = 5) {
  selected <- sample(1:length(population), k, replace = FALSE)
  parent_idx <- selected[which.min(fitness_scores[selected])]
  return(population[[parent_idx]])
}

select_parents_roulette <- function(population, fitness_scores) {
  # Calculate probabilities (inverse of fitness scores, assuming lower is better)
  probabilities <- 1 / fitness_scores
  probabilities <- probabilities / sum(probabilities)  # Normalize to sum to 1
  
  # Select a parent based on the probabilities
  parent_idx <- sample(1:length(population), 1, prob = probabilities)
  
  return(population[[parent_idx]])
}


################################################################################
crossover_trad <- function(parent1, parent2, Z, crossover_prob) {
  if (runif(1) < crossover_prob) {
    child1_X <- parent1$X_prime
    child1_Y <- parent1$Y_prime
    child2_X <- parent2$X_prime
    child2_Y <- parent2$Y_prime
    
    # Get the unique categories from Z
    categories <- unique(Z)
    
    # Pick two or more categories to swap samples between
    num_categories_to_swap <- sample(2:length(categories), 1) # Pick at least two categories
    selected_categories <- sample(categories, num_categories_to_swap)
    
    # For each pair of selected categories, swap samples
    for (i in 1:(num_categories_to_swap - 1)) {
      for (j in (i + 1):num_categories_to_swap) {
        cat1 <- selected_categories[i]
        cat2 <- selected_categories[j]
        
        # Indices for each category
        idx_cat1 <- which(Z == cat1)
        idx_cat2 <- which(Z == cat2)
        
        # Number of samples to swap (ensuring it's the same for both categories)
        num_samples_to_swap <- min(length(idx_cat1), length(idx_cat2))
        swap_idx_cat1 <- sample(idx_cat1, num_samples_to_swap)
        swap_idx_cat2 <- sample(idx_cat2, num_samples_to_swap)
        
        # Swap samples between categories cat1 and cat2
        child1_X[swap_idx_cat1] <- parent2$X_prime[swap_idx_cat2]
        child1_Y[swap_idx_cat1] <- parent2$Y_prime[swap_idx_cat2]
        child1_X[swap_idx_cat2] <- parent2$X_prime[swap_idx_cat1]
        child1_Y[swap_idx_cat2] <- parent2$Y_prime[swap_idx_cat1]
        
        child2_X[swap_idx_cat1] <- parent1$X_prime[swap_idx_cat2]
        child2_Y[swap_idx_cat1] <- parent1$Y_prime[swap_idx_cat2]
        child2_X[swap_idx_cat2] <- parent1$X_prime[swap_idx_cat1]
        child2_Y[swap_idx_cat2] <- parent1$Y_prime[swap_idx_cat1]
      }
    }
    
    return(list(list(X_prime = child1_X, Y_prime = child1_Y), list(X_prime = child2_X, Y_prime = child2_Y)))
  } else {
    return(list(parent1, parent2))
  }
}


# Crossover: Modify to consider correlation adjustment per category
crossover_trad_z <- function(parent1, parent2, Z, crossover_prob) {
  if (runif(1) < crossover_prob) {
    crossover_point <- sample(1:length(parent1$X_prime), 1)
    
    # Category-wise crossover: Ensure correlation is adjusted category-wise
    child1_X <- parent1$X_prime
    child1_Y <- parent1$Y_prime
    child2_X <- parent2$X_prime
    child2_Y <- parent2$Y_prime
    
    for (cat in unique(Z)) {
      idx <- which(Z == cat)
      
      # Combine values from both parents for this category
      crossover_idx <- sample(idx, round(length(idx) * 0.5))
      
      child1_X[crossover_idx] <- parent2$X_prime[crossover_idx] + rnorm(length(crossover_idx), 0, 0.1)
      child1_Y[crossover_idx] <- parent2$Y_prime[crossover_idx] + rnorm(length(crossover_idx), 0, 0.1)
      
      child2_X[crossover_idx] <- parent1$X_prime[crossover_idx] + rnorm(length(crossover_idx), 0, 0.1)
      child2_Y[crossover_idx] <- parent1$Y_prime[crossover_idx] + rnorm(length(crossover_idx), 0, 0.1)
    }
    
    return(list(list(X_prime = child1_X, Y_prime = child1_Y), list(X_prime = child2_X, Y_prime = child2_Y)))
  } else {
    return(list(parent1, parent2))
  }
}


# Crossover with Correlation Adjustment
crossover_corr_adjust <- function(parent1, parent2, target_corr, crossover_prob, corr_adjust_factor = 0.1) {
  # Calculate initial correlation between parent1 and parent2
  initial_corr <- cor(parent1$X_prime, parent1$Y_prime)
  
  if (runif(1) < crossover_prob) {
    # Blend the values from parent1 and parent2 with a small perturbation to adjust correlation
    blend_ratio <- runif(1)  # Random blend ratio between 0 and 1
    child1_X <- blend_ratio * parent1$X_prime + (1 - blend_ratio) * parent2$X_prime
    child1_Y <- blend_ratio * parent1$Y_prime + (1 - blend_ratio) * parent2$Y_prime
    
    child2_X <- blend_ratio * parent2$X_prime + (1 - blend_ratio) * parent1$X_prime
    child2_Y <- blend_ratio * parent2$Y_prime + (1 - blend_ratio) * parent1$Y_prime
    
    # Adjust correlation towards target
    new_corr <- cor(child1_X, child1_Y)
    corr_diff <- target_corr - new_corr
    adjustment <- corr_diff * corr_adjust_factor
    
    # Apply the adjustment to slightly modify X or Y to achieve target correlation
    child1_X <- child1_X + adjustment * rnorm(length(child1_X), 0, 0.1)
    child1_Y <- child1_Y + adjustment * rnorm(length(child1_Y), 0, 0.1)
    
    # Adjust correlation towards target
    new_corr <- cor(child2_X, child2_Y)
    corr_diff <- target_corr - new_corr
    adjustment <- corr_diff * corr_adjust_factor
    
    child2_X <- child2_X + adjustment * rnorm(length(child2_X), 0, 0.1)
    child2_Y <- child2_Y + adjustment * rnorm(length(child2_Y), 0, 0.1)
    
    return(list(list(X_prime = child1_X, Y_prime = child1_Y), list(X_prime = child2_X, Y_prime = child2_Y)))
  } else {
    return(list(parent1, parent2))
  }
}


################################################################################


# Mutation: Modify to adjust correlation toward target for each category
mutate_corr <- function(individual, Z, mutation_prob) {
  if (runif(1) < mutation_prob) {
    for (cat in unique(Z)) {
      idx <- which(Z == cat)
      curr_corr <- cor(individual$X_prime[idx], individual$Y_prime[idx])
      target_corr <- p[cat]
      
      # Adjust mutation based on correlation difference
      corr_diff <- target_corr - curr_corr
      adjustment <- rnorm(length(idx), mean = 0, sd = abs(corr_diff) * 0.5)
      
      # Mutate more in the direction of reducing correlation difference
      individual$X_prime[idx] <- individual$X_prime[idx] + adjustment
      individual$Y_prime[idx] <- individual$Y_prime[idx] + adjustment
    }
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


# Mutation with Correlation Adjustment
mutate_corr_adjust <- function(individual, target_corr, mutation_prob, corr_adjust_factor = 0.1) {
  if (runif(1) < mutation_prob) {
    # Calculate current correlation
    current_corr <- cor(individual$X_prime, individual$Y_prime)
    corr_diff <- target_corr - current_corr
    adjustment <- corr_diff * corr_adjust_factor
    
    # Apply distortion to X_prime and Y_prime based on correlation difference
    individual$X_prime <- individual$X_prime + adjustment * rnorm(length(individual$X_prime), 0, 0.1)
    individual$Y_prime <- individual$Y_prime + adjustment * rnorm(length(individual$Y_prime), 0, 0.1)
    
    # Introduce small boundary distortions by slightly shifting values near category boundaries
    category_shifts <- rbinom(length(individual$X_prime), 1, 0.1)  # Randomly decide which values to slightly shift
    shift_amounts <- runif(sum(category_shifts), -0.1, 0.1)  # Small random shifts
    
    individual$X_prime[category_shifts == 1] <- individual$X_prime[category_shifts == 1] + shift_amounts
    individual$Y_prime[category_shifts == 1] <- individual$Y_prime[category_shifts == 1] + shift_amounts
  }
  
  return(individual)
}


# select_parents_tournament
# select_parents_roulette
# crossover_trad
# crossover_trad_z
# crossover_corr_adjust
# mutate_corr
# mutate_with_copula
# mutate_corr_adjust
