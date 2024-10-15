library(tidyverse)
library(stringr)
library(gridExtra)
# install.packages("genodds")
library(genodds)
# install.packages("DescTools")
library(DescTools) # concord and siscord pairs

# Function to calculate mutual information
MutInf <- function(table) {
  rowsums <- rowSums(table)
  colsums <- colSums(table)
  ent_x <- sum(-rowsums/sum(rowsums) * log2(rowsums/sum(rowsums)), na.rm = TRUE)
  ent_y <- sum(-colsums/sum(colsums) * log2(colsums/sum(colsums)), na.rm = TRUE)
  return(ent_x + ent_y - entropy_pair(table))
}

# Function to calculate entropy of a pair
entropy_pair <- function(table) {
  total <- sum(table)
  probs <- table / total
  return(sum(-probs * log2(probs), na.rm = TRUE))
}

assoc_val <- function(tab){
  tab[1, 1]*tab[2, 2] - tab[1, 2]*tab[2, 1]
}

log_odds <- function(tab){
  log(tab[1, 1]/tab[1, 2]*(tab[1, 2] + tab[2, 2])/(tab[1, 1] + tab[2, 1]))
}

# Bar plot comparing log odds before and after transformation
plot_log_odds <- function(before_log_odds, after_log_odds, names_matrices) {
  # Add "Overall" to the list of matrices
  names_matrices <- c(names_matrices, "Overall")
  
  # Calculate overall log odds (sum of all matrices)
  overall_before <- log_odds_general(Reduce("+", matrices))
  overall_after <- log_odds_general(Reduce("+", result[[1]]))
  
  # Add overall log odds to the before and after values
  before_log_odds <- c(before_log_odds, overall_before)
  after_log_odds <- c(after_log_odds, overall_after)
  
  # Create a dataframe for plotting
  data <- data.frame(matrix_name = c(names_matrices, names_matrices),
                     log_odds = c(before_log_odds, after_log_odds),
                     type = c(rep("Before", length(names_matrices)), rep("After", length(names_matrices))))
  
  # Create the bar plot with before and after values
  data %>%
    mutate(type = factor(type, levels = c("Before", "After"), ordered = TRUE)) %>%
    ggplot(aes(x = matrix_name, y = log_odds, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    xlab("Matrix") +
    ylab("Log Odds") +
    ggtitle("Log Odds Before and After Transformation")
}

# Define the generalized log_odds function for n*m matrix
log_odds_general <- function(tab) {
  res = ConDisPairs(tab)[c("C","D")]
  concordant <- res$C
  discordant <- res$D
  if (discordant == 0) return(Inf)  # Handle case where discordant is 0 to avoid division by 0
  return(log(concordant / discordant))
}

log_odds <- log_odds_general



################################################################################
################################ Admissions Data ###############################
################################################################################

tm <- matrix(c(1198, 557, 1493, 1278), ncol = 2, byrow = T)
log_odds(tm)

ta <- matrix(c(512, 89, 313, 19), ncol = 2, byrow = T)
log_odds(ta)

tb <- matrix(c(353, 17, 207, 8), ncol = 2, byrow = T)
log_odds(tb)

tc <- matrix(c(120, 202, 205, 391), ncol = 2, byrow = T)
log_odds(tc)

td <- matrix(c(138, 131, 279, 244), ncol = 2, byrow = T)
log_odds(td)

te <- matrix(c(53, 94, 138, 299), ncol = 2, byrow = T)
log_odds(te)

tf <- matrix(c(22, 24, 351, 317), ncol = 2, byrow = T)
log_odds(tf)



# Function to calculate log odds by varying the matrix
calculate_log_odds <- function(mat, e) {
  log_odds_values <- numeric(length(e))
  
  for (i in seq_along(e)) {
    modified_mat_plus <- mat - matrix(c(e[i], -e[i], -e[i], e[i]), ncol = 2)
    log_odds_values[i] <- log_odds(modified_mat_plus)
  }
  
  return(log_odds_values)
}

# Define the matrices
matrices <- list(
  ta = matrix(c(512, 89, 313, 19), ncol = 2, byrow = TRUE),
  tb = matrix(c(353, 17, 207, 8), ncol = 2, byrow = TRUE),
  tc = matrix(c(120, 202, 205, 391), ncol = 2, byrow = TRUE),
  td = matrix(c(138, 131, 279, 244), ncol = 2, byrow = TRUE),
  te = matrix(c(53, 94, 138, 299), ncol = 2, byrow = TRUE),
  tf = matrix(c(22, 24, 351, 317), ncol = 2, byrow = TRUE),
  tm = matrix(c(1198, 557, 1493, 1278), ncol = 2, byrow = TRUE)
)

# Define the range for e
e_values <- seq(0, 20, by = 1)  # Adjust range as needed

# Initialize an empty data frame to store all results for single plot
all_data <- data.frame()

# Calculate log odds for each matrix and store in all_data for combined plot
for (name in names(matrices)) {
  if (name == "tm") {
    log_odds_values <- calculate_log_odds(matrices[[name]], 6 * e_values)
  } else {
    log_odds_values <- calculate_log_odds(matrices[[name]], e_values)
  }
  # Append data to all_data with matrix name
  df <- data.frame(e = e_values, log_odds = log_odds_values, matrix_name = name)
  all_data <- rbind(all_data, df)
}

# Create individual plots and store in a list for grid.arrange
plots <- list()
for (name in names(matrices)) {
  p <- ggplot(all_data %>% filter(matrix_name == name), aes(x = e, y = log_odds)) +
    geom_line() +
    ggtitle(paste("Log Odds for", name)) +
    xlab("e") +
    ylab("Log Odds") +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  plots[[name]] <- p
}

# Arrange all individual plots in a grid
grid.arrange(grobs = plots, ncol = 3)

# Plot all log odds on a single graph with legend
ggplot(all_data, aes(x = e, y = log_odds, color = matrix_name)) +
  geom_line() +
  ggtitle("Log Odds for All Matrices") +
  xlab("e") +
  ylab("Log Odds") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = rainbow(length(matrices)))

################################################################################
############################ Adusting Function (2*2) ###########################
################################################################################

# Define the log_odds function
log_odds <- function(tab){
  log(tab[1, 1]/tab[1, 2]*(tab[1, 2] + tab[2, 2])/(tab[1, 1] + tab[2, 1]))
}

# Function to augment a matrix by adding (phi, -phi; -phi, phi)
augment_matrix <- function(mat, phi) {
  augmented_mat <- mat + matrix(c(phi, -phi, -phi, phi), ncol = 2)
  return(augmented_mat)
}


# Define the adjust_e function
adjust_e <- function(matrices) {
  
  # Step 1: Calculate ei for each matrix
  ei_values <- sapply(matrices, function(ai) {
    floor((ai[1, 2] * ai[2, 1] - ai[1, 1] * ai[2, 2]) / sum(ai))
  })
  cat("ei_values", ei_values, "\n")
  
  # Step 2: Calculate em by summing all matrices element-wise
  summed_matrix <- Reduce("+", matrices)
  em <- ceiling((summed_matrix[1, 2] * summed_matrix[2, 1] - summed_matrix[1, 1] * summed_matrix[2, 2]) / sum(summed_matrix))
  cat("em", em, "\n")
  
  # Step 3: Check if sum(ei) >= em
  sum_ei <- sum(ei_values)
  cat("sum_ei", sum_ei, "\n")
  
  if (sum_ei < em) {
    cat("Impossible\n")
    return(NULL)
  }
  
  # Step 4: Calculate delta and phi
  delta <- sum_ei - em
  phi <- ceiling(delta / 2)
  cat("phi", phi, "\n")
  
  # Step 5: Augment matrices ensuring no negative or zero entries
  n <- length(matrices)
  augmented_matrices <- matrices  # Initialize with original matrices
  
  # Augment all matrices by phi
  if (phi >= n) {
    remaining_budget <- phi
    while(remaining_budget > 0){
      flag_check = 0
      for (i in 1:n) {
        min_val <- min(diag(augmented_matrices[[i]]))
        budget_used <- min(floor(phi/n), min_val - 1)
        augmented_matrices[[i]] <- augment_matrix(matrices[[i]], -budget_used)
        remaining_budget <- remaining_budget - budget_used
        if (budget_used == 0){flag_check = flag_check + 1}
      }
      if (flag_check == n){break}
    }
  } else {
    # Augment top phi matrices with lowest ei by (1, -1; -1, 1)
    sorted_indices <- order(ei_values)
    for (i in 1:phi) {
      idx <- sorted_indices[i]
      augmented_matrices[[idx]] <- augment_matrix(matrices[[idx]], -1)
    }
  }
  print(augmented_matrices)
  
  # Step 6: Create bar plot comparing log odds before and after transformation
  
  # Prepare data for plotting
  data <- data.frame(matrix_name = rep(c(names(matrices), "Overall"), each = 2),
                     log_odds = numeric(2 * (length(matrices) + 1)),
                     type = rep(c("Before", "After"), length(matrices) + 1))
  
  for (i in seq_along(matrices)) {
    before_log_odds <- log_odds(matrices[[i]])
    after_log_odds <- log_odds(augmented_matrices[[i]])
    
    data$log_odds[2 * i - 1] <- before_log_odds
    data$log_odds[2 * i] <- after_log_odds
  }
  
  # Log odds for the overall summed matrix
  overall_before <- log_odds(summed_matrix)
  overall_after <- log_odds(Reduce("+", augmented_matrices))
  
  data$log_odds[2 * length(matrices) + 1] <- overall_before
  data$log_odds[2 * length(matrices) + 2] <- overall_after
  
  # Create the bar plot
  data %>% 
    mutate(type = factor(type, levels = c("Before", "After"), ordered = TRUE)) %>%
    ggplot(aes(x = matrix_name, y = log_odds, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    xlab("Matrix") +
    ylab("Log Odds") +
    ggtitle("Log Odds Before and After Transformation")
}

# Example matrices
matrices <- list(
  ta = matrix(c(512, 89, 313, 19), ncol = 2, byrow = TRUE),
  tb = matrix(c(353, 17, 207, 8), ncol = 2, byrow = TRUE),
  tc = matrix(c(120, 202, 205, 391), ncol = 2, byrow = TRUE),
  td = matrix(c(138, 131, 279, 244), ncol = 2, byrow = TRUE),
  te = matrix(c(53, 94, 138, 299), ncol = 2, byrow = TRUE),
  tf = matrix(c(22, 24, 351, 317), ncol = 2, byrow = TRUE)
)

adjust_e(matrices)




################################################################################
############################ Adusting Function (n*m) ###########################
################################################################################

############################ Generalized Log-Odds Adjustment ##################

# Define the generalized log_odds function for n*m matrix
log_odds_general <- function(tab) {
  res = ConDisPairs(tab)[c("C","D")]
  concordant <- res$C
  discordant <- res$D
  if (discordant == 0) return(Inf)  # Handle case where discordant is 0 to avoid division by 0
  return(log(concordant / discordant))
}

# Define the generalized log_odds function for n*m matrix
# O(n^2)
log_odds_general <- function(tab) {
  res = ConDisPairs(tab)[c("C","D")]
  concordant <- res$C
  discordant <- res$D
  if (discordant == 0) return(Inf)  # Handle case where discordant is 0 to avoid division by 0
  return(log(concordant / discordant))
}

# Softmax function to compute probabilities based on log-odds
softmax <- function(x) {
  exp_x <- exp(x - max(x))  # Subtract max for numerical stability
  return(exp_x / sum(exp_x))
}

# Main function to adjust matrices based on generalized log-odds
adjust_matrices <- function(matrices, max_n = 1000, temp = 10) {
  
  # Define a function to sum all matrices
  sum_matrices <- function(mats) Reduce("+", mats)
  
  # Augment a random 2x2 block in the n*m matrix by (-D, D; D, -D) or (+D, -D; -D, +D)
  augment_matrix_random_block <- function(table, delta) {
    nrows <- nrow(table)
    ncols <- ncol(table)
    iters = 0
    repeat {
      iters = iters + 1
      RR <- sample(1:nrows, 2)
      CC <- sample(1:ncols, 2)
      if (delta < 0 && table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
      if (delta > 0 && table[RR[1], CC[2]] > 0 && table[RR[2], CC[1]] > 0) break
      if (iters > 100){delta = 0}
    }
    table[RR[1], CC[1]] <- table[RR[1], CC[1]] + delta
    table[RR[2], CC[2]] <- table[RR[2], CC[2]] + delta
    table[RR[1], CC[2]] <- table[RR[1], CC[2]] - delta
    table[RR[2], CC[1]] <- table[RR[2], CC[1]] - delta
    
    return(table)
  }
  
  # Initialize the current and best matrices and evaluations
  curr <- matrices
  curr_eval_sum <- log_odds_general(sum_matrices(curr))
  
  # Store history for analysis later
  history <- data.frame(iteration = 0, overall_log_odds = curr_eval_sum)
  
  # Simulated annealing loop
  for (n in 1:max_n) {
    
    # Compute log-odds for all matrices
    log_odds_values <- sapply(curr, log_odds_general)
    
    # Decide order of selection based on current overall evaluation
    if (curr_eval_sum >= 0) {
      selection_probs <- softmax(log_odds_values)  # Higher log-odds, higher probability
    } else {
      selection_probs <- softmax(-log_odds_values)  # Lower log-odds, higher probability
    }
    print(selection_probs)
    
    # Select a matrix based on softmax probabilities
    idx <- sample(1:length(curr), 1, prob = selection_probs)
    
    # current evals
    curr_eval_indiv <- log_odds_general(curr[[idx]])
    
    # Decide the direction of the adjustment based on current log-odds
    delta <- ifelse(curr_eval_sum >= 0, -1, 1)
    
    # Augment the selected matrix
    modified_matrix <- augment_matrix_random_block(curr[[idx]], delta)
    
    # Calculate new overall log-odds after modification
    new_matrices <- curr
    new_matrices[[idx]] <- modified_matrix
    
    # new evals
    new_eval_sum <- log_odds_general(sum_matrices(new_matrices))
    new_eval_indiv <- log_odds_general(new_matrices[[idx]])
    
    # Accept modification if overall log-odds improves or passes Metropolis criteria
    temp <- 0.9 * temp 
    metropolis <- ifelse(curr_eval_sum >= 0, exp((curr_eval_sum - new_eval_sum) / temp),
                         exp((new_eval_sum - curr_eval_sum)/temp))
    
    # Conditions to accept/reject the change
    if (delta == -1 && new_eval_sum > 0 && runif(1) < metropolis) {
      curr <- new_matrices
      curr_eval_sum <- new_eval_sum
    } else if (delta == 1 && new_eval_indiv < 0 && runif(1) < metropolis) {
      curr <- new_matrices
      curr_eval_sum <- new_eval_sum
    }
    
    # Track the progress
    history <- rbind(history, data.frame(iteration = n, overall_log_odds = curr_eval_sum))
  }
  
  # Return the final matrices and log odds history
  return(list(curr, curr_eval_sum, history))
}

# Bar plot comparing log odds before and after transformation
plot_log_odds <- function(before_log_odds, after_log_odds, names_matrices) {
  # Add "Overall" to the list of matrices
  names_matrices <- c(names_matrices, "Overall")
  
  # Calculate overall log odds (sum of all matrices)
  overall_before <- log_odds_general(Reduce("+", matrices))
  overall_after <- log_odds_general(Reduce("+", result[[1]]))
  
  # Add overall log odds to the before and after values
  before_log_odds <- c(before_log_odds, overall_before)
  after_log_odds <- c(after_log_odds, overall_after)
  
  # Create a dataframe for plotting
  data <- data.frame(matrix_name = c(names_matrices, names_matrices),
                     log_odds = c(before_log_odds, after_log_odds),
                     type = c(rep("Before", length(names_matrices)), rep("After", length(names_matrices))))
  
  # Create the bar plot with before and after values
  data %>%
    mutate(type = factor(type, levels = c("Before", "After"), ordered = TRUE)) %>%
    ggplot(aes(x = matrix_name, y = log_odds, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    xlab("Matrix") +
    ylab("Log Odds") +
    ggtitle("Log Odds Before and After Transformation")
}

set.seed(42)
# Example matrix: 1
matrices <- list(
  ta = matrix(c(512, 89, 313, 19), ncol = 2, byrow = TRUE),
  tb = matrix(c(353, 17, 207, 8), ncol = 2, byrow = TRUE),
  tc = matrix(c(120, 202, 205, 391), ncol = 2, byrow = TRUE),
  td = matrix(c(138, 131, 279, 244), ncol = 2, byrow = TRUE),
  te = matrix(c(53, 94, 138, 299), ncol = 2, byrow = TRUE),
  tf = matrix(c(22, 24, 351, 317), ncol = 2, byrow = TRUE)
)

# Run the generalized log-odds adjustment on matrices
result <- adjust_matrices(matrices, max_n = 200)

# Extract log odds before and after
before_log_odds <- sapply(matrices, log_odds_general)
after_log_odds <- sapply(result[[1]], log_odds_general)

# Plot log odds before and after
plot_log_odds(before_log_odds, after_log_odds, names(matrices))

# Plot the trajectory of log odds during the adjustment process
result[[3]] %>%
  ggplot(aes(x = iteration, y = overall_log_odds)) +
  geom_line() +
  theme_minimal() +
  xlab("Iteration") +
  ylab("Overall Log Odds") +
  ggtitle("Log Odds Trajectory During Adjustment")


set.seed(42)
# Example matrix: 2
matrices <- list(
  ta = matrix(c(512, 89, 313, 19, 45, 67), ncol = 3, byrow = TRUE),
  tb = matrix(c(353, 17, 207, 8, 33, 25), ncol = 3, byrow = TRUE),
  tc = matrix(c(120, 202, 205, 391, 55, 23), ncol = 3, byrow = TRUE),
  td = matrix(c(138, 131, 279, 244, 111, 99), ncol = 3, byrow = TRUE),
  te = matrix(c(53, 94, 138, 299, 87, 77), ncol = 3, byrow = TRUE),
  tf = matrix(c(22, 24, 351, 317, 11, 33), ncol = 3, byrow = TRUE),
  tg = matrix(c(22, 24, 351, 317, 11, 33), ncol = 3, byrow = TRUE)
)

# Run the generalized log-odds adjustment on matrices
result <- adjust_matrices(matrices, max_n = 1500)

# Extract log odds before and after
before_log_odds <- sapply(matrices, log_odds_general)
after_log_odds <- sapply(result[[1]], log_odds_general)

# Plot log odds before and after
plot_log_odds(before_log_odds, after_log_odds, names(matrices))

# Plot the trajectory of log odds during the adjustment process
result[[3]] %>%
  ggplot(aes(x = iteration, y = overall_log_odds)) +
  geom_line() +
  theme_minimal() +
  xlab("Iteration") +
  ylab("Overall Log Odds") +
  ggtitle("Log Odds Trajectory During Adjustment")

# 1: other metric
# indiv specification
# check possibility for 2*2



