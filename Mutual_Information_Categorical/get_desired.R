library(tidyverse)
library(stringr)

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

# Function to generate a new number for maximizing mutual information
gen_number_max <- function(x, tot_sum) {
  table <- x
  nrows <- nrow(table)
  ncols <- ncol(table)
  repeat {
    RR <- sample(1:nrows, 2)
    CC <- sample(1:ncols, 2)
    if (table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
  }
  
  delta1 <- min(table[RR[1], CC[1]], table[RR[2], CC[2]])
  delta2 <- -min(table[RR[1], CC[2]], table[RR[2], CC[1]])
  
  table1 <- table
  table2 <- table
  
  table1[RR[1], CC[1]] <- table1[RR[1], CC[1]] - delta1
  table1[RR[2], CC[2]] <- table1[RR[2], CC[2]] - delta1
  table1[RR[1], CC[2]] <- table1[RR[1], CC[2]] + delta1
  table1[RR[2], CC[1]] <- table1[RR[2], CC[1]] + delta1
  
  table2[RR[1], CC[1]] <- table2[RR[1], CC[1]] - delta2
  table2[RR[2], CC[2]] <- table2[RR[2], CC[2]] - delta2
  table2[RR[1], CC[2]] <- table2[RR[1], CC[2]] + delta2
  table2[RR[2], CC[1]] <- table2[RR[2], CC[1]] + delta2
  
  mutinf1 <- MutInf(table1)
  mutinf2 <- MutInf(table2)
  
  if (mutinf1 > mutinf2) {
    return(table1)
  } else {
    return(table2)
  }
}

# Function to generate a new number for minimizing mutual information
gen_number_min <- function(x, tot_sum) {
  table <- x
  nrows <- nrow(table)
  ncols <- ncol(table)
  repeat {
    RR <- sample(1:nrows, 2)
    CC <- sample(1:ncols, 2)
    if (table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
  }
  
  S <- (table[RR[1], CC[1]] + table[RR[2], CC[2]] + table[RR[1], CC[2]] + table[RR[2], CC[1]])*tot_sum
  delta <- round((table[RR[1], CC[1]]*tot_sum^2*table[RR[2], CC[2]] - table[RR[1], CC[2]]*tot_sum^2*table[RR[2], CC[1]]) / S)
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] - delta/tot_sum
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] - delta/tot_sum
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] + delta/tot_sum
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] + delta/tot_sum
  
  return(table)
}

# Function to generate a new number for stepwise modification
gen_number_1 <- function(x, tot_sum) {
  table <- x
  nrows <- nrow(table)
  ncols <- ncol(table)
  repeat {
    RR <- sample(1:nrows, 2)
    CC <- sample(1:ncols, 2)
    if (table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
  }
  
  delta <- 1/tot_sum
  
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] - delta
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] - delta
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] + delta
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] + delta
  
  return(table)
}

# Function to perform simulated annealing with stopping condition for target entropy
simul_anneal <- function(initial_table, obj, gen_fn, tot_sum, target, max_n = 5000, temp = 10, maxim = T, readj = F) {
  best <- initial_table
  best_eval <- obj(best)
  curr <- best
  curr_eval <- best_eval
  n <- 0
  mutual_info_history <- data.frame(iteration = 0, mutual_info = curr_eval, 
                                    type = ifelse(readj, "Readjusting", ifelse(maxim, "Maximizing", "Minimizing")))
  
  while (n < max_n) {
    if ((maxim && curr_eval >= target) || (!maxim && curr_eval <= target)) {
      print(paste("Target entropy reached:", curr_eval))
      break
    }
    
    cand <- gen_fn(curr, tot_sum)
    cand_eval <- obj(cand)
    if (maxim){
      if (cand_eval > best_eval) {
        best <- cand
        best_eval <- cand_eval
      }
      diff <- cand_eval - curr_eval
    } else {
      if (cand_eval < best_eval) {
        best <- cand
        best_eval <- cand_eval
      }
      diff <- -cand_eval + curr_eval
    }
    
    temp <- 0.9*temp #/ (n + 1)
    metropolis <- exp(diff / temp)
    if (diff > 0 || runif(1) < metropolis) {
      curr <- cand
      curr_eval <- cand_eval
    }
    n <- n + 1
    mutual_info_history <- rbind(mutual_info_history, data.frame(iteration = n, mutual_info = curr_eval, 
                                                                 type =  ifelse(readj, "Readjusting", ifelse(maxim, "Maximizing", "Minimizing"))))
  }
  return(list(best, best_eval, n, mutual_info_history))
}

# Sinkhorn algorithm function
sinkhorn_algorithm <- function(initial_table, obj, max_iter = 500, tolerance = 1e-5) {
  S <- sum(initial_table)
  S_r <- rowSums(initial_table)
  S_c <- colSums(initial_table)
  n <- nrow(initial_table)
  m <- ncol(initial_table)
  
  # Initialize alpha and beta
  alpha <- rep(1, m)
  beta <- rep(1, n)
  
  curr_eval <- obj(initial_table)
  mutual_info_history <- data.frame(iteration = 0, mutual_info = curr_eval, 
                                    type = "Minimizing")
  
  for (iter in 1:max_iter) {
    
    ones_mat <- matrix(1, nrow = m, ncol = n)
    
    # Update alpha
    alpha <- (S_c / (ones_mat %*% beta))
    # Update beta
    beta <- (S_r / (t(ones_mat) %*% alpha))
    
    # Convergence check
    updated_table <- t(diag(as.vector(alpha)) %*% ones_mat %*% diag(as.vector(beta)))
    curr_eval <- obj(updated_table)
    
    mutual_info_history <- rbind(mutual_info_history, data.frame(iteration = iter, mutual_info = curr_eval, 
                                                                 type = "Minimizing"))
  }
  updated_table <- t(diag(as.vector(alpha)) %*% ones_mat %*% diag(as.vector(beta)))
  new_mut <- obj(updated_table)
  return(list(updated_table, new_mut, iter, mutual_info_history))
}

# Function to solve entropy by adjusting mutual information
solve_entropy <- function(df, target, max_n = 10000, epsilon = 0.001) {
  y_name <- colnames(df)[1]
  x_name <- colnames(df)[2]
  table <- as.matrix(table(df[, c(y_name, x_name)]))
  table_sum <- sum(table)
  table <- table/table_sum # normalize
  old_mut <- MutInf(table)
  
  # Step 1: Find range of entropy values
  
  # Max Entropy (Maximizing mutual information)
  result_max <- simul_anneal(table, obj = MutInf, gen_fn = gen_number_max, 
                             tot_sum = table_sum, target = Inf, max_n = max_n, maxim = TRUE)
  new_mut_max <- result_max[[2]]
  
  # Min Entropy (Minimizing mutual information)
  result_min_sk <- sinkhorn_algorithm(table, obj = MutInf, max_iter = max_n)
  new_mut_min <- result_min_sk[[2]]
  
  # Print the max and min entropy values
  print(paste("Current Entropy:", old_mut))
  print(paste("Min Entropy:", new_mut_min))
  print(paste("Max Entropy:", new_mut_max))
  
  # Step 2: Check if the target entropy is within the range
  if (target < new_mut_min || target > new_mut_max) {
    print("Target entropy is out of range. Please choose a value between the min and max entropy.")
    return(NULL)
  }
  
  # Step 3: Adjust mutual information to reach the target entropy
  if (target > old_mut) {
    gen_fn <- gen_number_max
    result <- simul_anneal(table, obj = MutInf, gen_fn = gen_fn, 
                           tot_sum = table_sum, target = target, max_n = max_n, maxim = TRUE)
    final_hist = result[[4]]
    final_table <- result[[1]]
    final_mut <- result[[2]]
    if (result[[2]] - target > epsilon) {
      print("Target exceeded. Re-adjusting with gen_number_1 to decrease.")
      result_sub <- simul_anneal(result[[1]], obj = MutInf, gen_fn = gen_number_1, 
                             tot_sum = table_sum, target = target, max_n = max_n, 
                             maxim = FALSE, readj = T)
      result_sub[[4]]$iteration = result_sub[[4]]$iteration + max(final_hist$iteration)
      final_hist = rbind(result_sub[[4]], final_hist)
      final_table <- result_sub[[1]]
      final_mut <- result_sub[[2]]
    }
    
  } else {
    gen_fn <- gen_number_min
    result <- simul_anneal(table, obj = MutInf, gen_fn = gen_fn, 
                           tot_sum = table_sum, target = target, max_n = max_n, maxim = FALSE)
    final_hist = result[[4]]
    final_table <- result[[1]]
    final_mut <- result[[2]]
    if (target - result[[2]] > epsilon ) {
      print("Target crossed. Re-adjusting with gen_number_1 to increase.")
      result_sub <- simul_anneal(result[[1]], obj = MutInf, gen_fn = gen_number_1, 
                             tot_sum = table_sum, target = target, max_n = max_n, 
                             maxim = TRUE, readj = T)
      result_sub[[4]]$iteration = result_sub[[4]]$iteration + max(final_hist$iteration)
      final_hist = rbind(result_sub[[4]], final_hist)
      final_table <- result_sub[[1]]
      final_mut <- result_sub[[2]]
    }
  }
  
  print(paste("Final Mutual Information:", final_mut))
  return(list(final_table = final_table, history = final_hist,
              max_mut = new_mut_max, min_mut = new_mut_min))
}

# Increasing
set.seed(33)
df <- data.frame(
  x = sample(str_c("Categ", 1:4), 10000, replace = TRUE),
  y = sample(str_c("Categ", 10:4), 10000, replace = TRUE)
)

target_entropy <- 1.4 # Set your target entropy here
res <- solve_entropy(df, target_entropy)
final_table <- res$final_table
a = res$history
res$history %>%
  ggplot(aes(x = iteration, y = mutual_info, color = type)) +
  geom_line() +
  labs(title = "Mutual Information Over Iterations", x = "Iteration", y = "Mutual Information") +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 60)) +
  xlim(c(0, 1.5 * max(res$history$iteration))) + 
  ylim(c(res$min_mut - 1, res$max_mut + 1)) +
  geom_hline(yintercept = res$min_mut, linetype = "dashed", color = "black") +
  geom_hline(yintercept = res$max_mut, linetype = "dashed", color = "black") +
  geom_point(data = subset(res$history, iteration == max(res$history$iteration)), 
             aes(x = iteration, y = mutual_info), 
             color = "black", size = 1) +
  annotate("text", x = max(res$history$iteration), y = res$min_mut, 
           label = paste0("Min Mutual Info: ", round(res$min_mut, 3)), 
           vjust = 1.6, hjust = 1, color = "black", size = 3) +
  annotate("text", x = max(res$history$iteration), y = res$max_mut, 
           label = paste0("Max Mutual Info: ", round(res$max_mut, 3)), 
           vjust = -1, hjust = 1, color = "black", size = 3) +
  scale_x_log10()


# Decreasing


get_entropy <- function(df) {
  y_name <- colnames(df)[1]
  x_name <- colnames(df)[2]
  table <- as.matrix(table(df[, c(y_name, x_name)]))
  table_sum <- sum(table)
  table <- table/table_sum # normalize
  old_mut <- MutInf(table)
  old_mut
}
set.seed(42)
df <- data.frame(
  x = sample(str_c("Categ", 1:4), 10000, replace = TRUE)
) %>% mutate(y = 
               ifelse(str_ends(x, "1"), sample(str_c("Categ", 10:8), 1, replace = TRUE),
                      ifelse(str_ends(x, "2"), sample(str_c("Categ", 8:6), 1, replace = TRUE),
                             sample(str_c("Categ", 5:4), 1, replace = TRUE))
                      )
             )
get_entropy(df)
target_entropy <- 0.2 # Set your target entropy here
res <- solve_entropy(df, target_entropy)
final_table <- res$final_table
a = res$history
res$history %>%
  ggplot(aes(x = iteration, y = mutual_info, color = type)) +
  geom_line() +
  labs(title = "Mutual Information Over Iterations", x = "Iteration", y = "Mutual Information") +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 60)) +
  xlim(c(0, 1.5 * max(res$history$iteration))) + 
  ylim(c(res$min_mut - 1, res$max_mut + 1)) +
  geom_hline(yintercept = res$min_mut, linetype = "dashed", color = "black") +
  geom_hline(yintercept = res$max_mut, linetype = "dashed", color = "black") +
  geom_point(data = subset(res$history, iteration == max(res$history$iteration)), 
             aes(x = iteration, y = mutual_info), 
             color = "black", size = 1) +
  annotate("text", x = max(res$history$iteration), y = res$min_mut, 
           label = paste0("Min Mutual Info: ", round(res$min_mut, 3)), 
           vjust = 1.6, hjust = 1, color = "black", size = 3) +
  annotate("text", x = max(res$history$iteration), y = res$max_mut, 
           label = paste0("Max Mutual Info: ", round(res$max_mut, 3)), 
           vjust = -1, hjust = 1, color = "black", size = 3) +
  scale_x_log10()


