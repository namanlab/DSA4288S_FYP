library(dplyr)

# Function to modify the dataframe and calculate new mutual information using Sinkhorn algorithm
new_mutual_cols <- function(df, max_n = 500, tolerance = 1e-5) {
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
  
  # Perform Sinkhorn algorithm
  result <- sinkhorn_algorithm(table, max_iter = max_n, tolerance = tolerance)
  table <- result[[1]]
  new_mut <- result[[2]]
  # New row & column sums
  rowsum_new <- rowSums(table)
  colsum_new <- colSums(table)
  
  print(sum(abs(rowsum_old - rowsum_new))  +  sum(abs(colsum_old - colsum_new)))
  print(sum(table))
  
  
  
  # Return Dataframe
  # Print old and new mutual information values
  print(paste("Old Mutual Information:", old_mut))
  print("New Table:")
  print(table)
  print(paste("New Mutual Information:", new_mut))
  
  return(list(old_mut, new_mut, result[[3]]))
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

sinkhorn_algorithm <- function(initial_table, max_iter = 500, tolerance = 1e-5) {
  S <- sum(initial_table)
  S_r <- rowSums(initial_table)
  S_c <- colSums(initial_table)
  n <- nrow(initial_table)
  m <- ncol(initial_table)
  
  # Initialize alpha and beta
  alpha <- rep(1, m)
  beta <- rep(1, n)
  
  for (iter in 1:max_iter) {
    
    ones_mat <- matrix(1, nrow = m, ncol = n)
    
    # Update alpha
    alpha <- (S_c / (ones_mat %*% beta))
    # Update beta
    beta <- (S_r / (t(ones_mat) %*% alpha))
    
    # Convergence check
    updated_table <- t(diag(as.vector(alpha)) %*% ones_mat %*% diag(as.vector(beta)))
    new_mut <- MutInf(updated_table)
    
    # Check for convergence
    if (max(abs(rowSums(updated_table) - S_r)) < tolerance && 
        max(abs(colSums(updated_table) - S_c)) < tolerance) {
      break
    }
  }
  
  updated_table <- t(diag(as.vector(alpha)) %*% ones_mat %*% diag(as.vector(beta)))
  return(list(updated_table, new_mut, iter))
}


# Example usage
set.seed(42)
df <- data.frame(
  x = sample(str_c("Categ", 1:100), 10000, replace = TRUE),
  y = sample(str_c("Categ", 10:100), 10000, replace = TRUE)
)


# df <- data.frame(
#   x = sample(c(rep("b", 60), rep("a", 15)), 75),
#   y = sample(c(rep("b", 25), rep("a", 50)), 75)
# )

# Find range of mutual information
res <- new_mutual_cols(df, max_n = 10000)