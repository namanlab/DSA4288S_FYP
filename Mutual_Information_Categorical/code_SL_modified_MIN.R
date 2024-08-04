library(dplyr)

# Function to modify the dataframe and calculate new mutual information
new_mutual_cols <- function(df, max_n = 500, min_inf = 0.5) {
  df <- df
  # Record label names
  y_name <- colnames(df)[1]
  x_name <- colnames(df)[2]
  # Get contingency table
  table <- table(df[, c(y_name, x_name)])
  x_index <- rownames(table)
  y_index <- colnames(table)
  # Convert table to matrix
  table <- as.matrix(table)
  print("Old Table:")
  print(table)
  # Record old row & column sums
  rowsum_old <- rowSums(table)
  colsum_old <- colSums(table)
  old_mut <- MutInf(table)
  
  # Perform simulated annealing
  result <- simul_anneal(table, max_n = max_n, min_inf = min_inf)
  table <- result[[1]]
  new_mut <- result[[2]]
  # New row & column sums
  rowsum_new <- rowSums(table)
  colsum_new <- colSums(table)
  
  if (!identical(rowsum_old, rowsum_new) || !identical(colsum_old, colsum_new)) {
    stop("Marginals not the same")
  }
  
  # Convert to dataframe to return
  df <- as.table(table)
  # Reset field names
  rownames(df) <- x_index
  colnames(df) <- y_index
  
  # Unpivot frequency table
  df = as.data.frame(table)
  df <- df[rep(row.names(df), df$Freq),]
  df <- df[order(row.names(df)), ] 
  
  # Return Dataframe
  # Print old and new mutual information values
  print(paste("Old Mutual Information:", old_mut))
  print("New Table:")
  print(table)
  print(paste("New Mutual Information:", new_mut))
  
  return(list(df, old_mut, new_mut, result[[3]]))
}

# Function to calculate mutual information
MutInf <- function(table) {
  # Get the entropy of the marginals
  rowsums <- rowSums(table)
  colsums <- colSums(table)
  ent_x <- sum(-rowsums/sum(rowsums) * log2(rowsums/sum(rowsums)), na.rm = T)
  ent_y <- sum(-colsums/sum(colsums) * log2(colsums/sum(colsums)), na.rm = T)
  # Get mutual information
  return(ent_x + ent_y - entropy_pair(table))
}

# Function to get row and column sums
get_row_col_sums <- function(x) {
  return(list(rowSums(x), colSums(x)))
}

# Function to calculate entropy of a pair
entropy_pair <- function(table) {
  total <- sum(table)
  probs <- table / total
  return(sum(-probs * log2(probs), na.rm = T))
}

# Function to generate a new number
gen_number <- function(x) {
  table <- x
  nrows <- nrow(table)
  ncols <- ncol(table)
  RR <- sample(1:nrows, 2)
  CC <- sample(1:ncols, 2)
  org_values <- c(table[RR[1], CC[1]], table[RR[2], CC[2]])
  
  while (org_values[1] < 1 || org_values[2] < 1) {
    RR <- sample(1:nrows, 2)
    CC <- sample(1:ncols, 2)
    org_values <- c(table[RR[1], CC[1]], table[RR[2], CC[2]])
  }
  
  S <- table[RR[1], CC[1]] + table[RR[2], CC[2]] + table[RR[1], CC[2]] + table[RR[2], CC[1]]
  delta <- round((table[RR[1], CC[1]]*table[RR[2], CC[2]] - table[RR[1], CC[2]]*table[RR[2], CC[1]])/S)
  
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] - delta
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] - delta
  
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] + delta
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] + delta
  
  return(table)
}

# Function to perform simulated annealing
simul_anneal <- function(initial_table, obj = MutInf, gen_fn = gen_number, max_n = 5000, min_inf = 0.5, temp = 10) {
  # Copy of initial point
  best <- initial_table
  # Evaluate Initial point
  best_eval <- obj(best)
  # Current Working Solution
  curr <- best
  curr_eval <- best_eval
  n <- 0
  # Run the algorithm while the mutual information isn't sufficient and the max number of iterations hasn't been reached
  while (n < max_n) {
    # Take a step
    cand <- gen_fn(curr)
    # Evaluate new candidate point
    cand_eval <- obj(cand)
    # Check for new best solution
    if (cand_eval < best_eval) {
      best <- cand
      best_eval <- cand_eval
    }
    # Difference between candidate and current point evaluation
    diff <- - cand_eval + curr_eval
    # Calculate temperature
    t <- temp / (n + 1)
    # Metropolis acceptance criterion
    metropolis <- exp(diff / t)
    if (diff > 0 || runif(1) < metropolis) {
      curr <- cand
      curr_eval <- cand_eval
    }
    n <- n + 1
  }
  return(list(best, best_eval, n))
}



# Example usage
set.seed(42)
df <- data.frame(
  x = sample(LETTERS[1:20], 1000, replace = TRUE),
  y = sample(LETTERS[6:25], 1000, replace = TRUE)
)


# df <- data.frame(
#   x = sample(c(rep("b", 40), rep("a", 35)), 75),
#   y = sample(c(rep("b", 15), rep("a", 60)), 75)
# )

result <- new_mutual_cols(df)
print(result[[4]])

