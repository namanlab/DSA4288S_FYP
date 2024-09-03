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

# Function to perform simulated annealing
simul_anneal <- function(initial_table, obj, gen_fn, tot_sum, max_n = 5000, temp = 10, maxim = T) {
  best <- initial_table
  best_eval <- obj(best)
  curr <- best
  curr_eval <- best_eval
  n <- 0
  mutual_info_history <- data.frame(iteration = 0, mutual_info = curr_eval, 
                                    type = if (maxim) "Maximizing" else "Minimizing")
  
  while (n < max_n) {
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

    t <- temp / (n + 1)
    metropolis <- exp(diff / t)
    if (diff > 0 || runif(1) < metropolis) {
      curr <- cand
      curr_eval <- cand_eval
    }
    n <- n + 1
    mutual_info_history <- rbind(mutual_info_history, data.frame(iteration = n, mutual_info = curr_eval, type = if (maxim) "Maximizing" else "Minimizing"))
  }
  return(list(best, best_eval, n, mutual_info_history))
}

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
    # updated_table <-  round_vals(updated_table, initial_table)
    curr_eval <- obj(updated_table)
    
    mutual_info_history <- rbind(mutual_info_history, data.frame(iteration = iter, mutual_info = curr_eval, 
                                                                 type = "Minimizing"))
    
  }
  updated_table <- t(diag(as.vector(alpha)) %*% ones_mat %*% diag(as.vector(beta)))
  # updated_table <-  round_vals(updated_table, initial_table)
  new_mut <- obj(updated_table)
  return(list(updated_table, new_mut, iter, mutual_info_history))
}

round_vals <- function(final, orig){
  rounded_final <<- final
  table <- orig
  nrows <- nrow(table)
  ncols <- ncol(table)
  for (i in 1:(nrows - 1)){
    for (j in 1:(ncols - 1)){
      delta <- table[i, j] - round(final[i, j])
      table[i, j] <- table[i, j] - delta
      table[i+1, j+1] <- table[i+1, j+1]- delta
      table[i, j+1] <- table[i, j+1] + delta
      table[i+1, j] <- table[i+1, j] + delta
    }
  }
  return(table)
}


# Function to modify the dataframe and calculate new mutual information
new_mutual_cols <- function(df, max_n = 10000, min_inf = 0.5) {
  y_name <- colnames(df)[1]
  x_name <- colnames(df)[2]
  table <- as.matrix(table(df[, c(y_name, x_name)]))
  table_sum <- sum(table)
  table <- table/table_sum # normalize
  x_index <- rownames(table)
  y_index <- colnames(table)
  rowsum_old <- rowSums(table)
  colsum_old <- colSums(table)
  old_mut <- MutInf(table)
  
  # MAX 
  gen_fn <- gen_number_max
  result_max <- simul_anneal(table, obj = MutInf, gen_fn = gen_number_max, 
                             tot_sum = table_sum, max_n = max_n, maxim = T)
  table_max <- result_max[[1]]
  new_mut_max <- result_max[[2]]
  rowsum_new_max <- rowSums(table_max)
  colsum_new_max <- colSums(table_max)
  print(sum(abs(rowsum_old - rowsum_new_max)) + sum(abs(colsum_old - colsum_new_max)))
  max_mutual_info_history <- result_max[[4]]
  
  # MIN 
  gen_fn <- gen_number_min
  result_min <- simul_anneal(table, obj = MutInf, gen_fn = gen_number_min, 
                             tot_sum = table_sum, max_n = max_n, maxim = F)
  table_min <- result_min[[1]]
  new_mut_min <- result_min[[2]]
  rowsum_new_min <- rowSums(table_min)
  colsum_new_min <- colSums(table_min)
  print(sum(abs(rowsum_old - rowsum_new_min)) + sum(abs(colsum_old - colsum_new_min)))
  min_mutual_info_history <- result_min[[4]]
  
  # MAX using gen_number_1
  gen_fn <- gen_number_1
  result_max1 <- simul_anneal(table, obj = MutInf, gen_fn = gen_number_1, 
                              tot_sum = table_sum, max_n = max_n, maxim = TRUE)
  table_max1 <- result_max1[[1]]
  new_mut_max1 <- result_max1[[2]]
  rowsum_new_max1 <- rowSums(table_max1)
  colsum_new_max1 <- colSums(table_max1)
  print(sum(abs(rowsum_old - rowsum_new_max1)) + sum(abs(colsum_old - colsum_new_max1)))
  max1_mutual_info_history <- result_max1[[4]]
  
  # MIN using gen_number_1
  gen_fn <- gen_number_1
  result_min1 <- simul_anneal(table, obj = MutInf, gen_fn = gen_number_1, 
                              tot_sum = table_sum, max_n = max_n, maxim = FALSE)
  table_min1 <- result_min1[[1]]
  new_mut_min1 <- result_min1[[2]]
  rowsum_new_min1 <- rowSums(table_min1)
  colsum_new_min1 <- colSums(table_min1)
  print(sum(abs(rowsum_old - rowsum_new_min1)) + sum(abs(colsum_old - colsum_new_min1)))
  min1_mutual_info_history <- result_min1[[4]]
  
  # SINKHORN 
  result_min_sk <- sinkhorn_algorithm(table, obj = MutInf, max_iter = max_n)
  table_min_sk <- result_min_sk[[1]]
  new_mut_min_sk <- result_min_sk[[2]]
  rowsum_new_min_sk <- rowSums(table_min_sk)
  colsum_new_min_sk <- colSums(table_min_sk)
  print(sum(abs(rowsum_old - rowsum_new_min_sk))  +  sum(abs(colsum_old - colsum_new_min_sk)))
  min_mutual_info_history_sk <- result_min_sk[[4]]
  
  
  mutual_info_history <- rbind(max_mutual_info_history, min_mutual_info_history)
  mutual_info_history$method = "Improved"
  mutual1_info_history <- rbind(max1_mutual_info_history, min1_mutual_info_history)
  mutual1_info_history$method = "Stepwise"
  min_mutual_info_history_sk$method = "Sinkhorn"
  df_stat = rbind(mutual_info_history, mutual1_info_history, min_mutual_info_history_sk)
  
  
  print("Old Table:")
  print(table)
  print(paste("Old Mutual Information:", old_mut))
  # print("New Max Table:")
  # print(table_max)
  print(paste("Max Mutual Information:", new_mut_max))
  # print("New Min Table:")
  # print(table_min)
  print(paste("Min Mutual Information:", new_mut_min))
  # print("New Sinkhorn Min Table:")
  # print(table_min_sk)
  print(paste("Min Sinkhorn Mutual Information:", new_mut_min_sk))
  
  return(list(old_mut = old_mut, new_mut_max = new_mut_max, new_mut_min = new_mut_min, 
              table = table, table_min = table_min, table_max = table_max,table_min_sk = table_min_sk,
              df_stat = df_stat))
}

# Example usage: LARGE
set.seed(42)
df <- data.frame(
  x = sample(str_c("Categ", 1:50), 10000, replace = TRUE),
  y = sample(str_c("Categ", 11:60), 10000, replace = TRUE)
)
res <- new_mutual_cols(df, max_n = 50000)
res$df_stat  %>% 
  ggplot(aes(x = iteration, y = mutual_info, color = method)) +
  geom_line() +
  labs(title = "Mutual Information Over Iterations", x = "Iteration", y = "Mutual Information") +
  theme_minimal() +
  facet_wrap(~type, scales = "free_y") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 60))

# Example usage: MEDIUM
set.seed(9)
df <- data.frame(
  x = sample(str_c("Categ", 1:10), 10000, replace = TRUE)
) %>% mutate(y = 
               ifelse(str_ends(x, "5"), sample(str_c("Categ", 11:20), replace = TRUE),
                             sample(str_c("Categ", 11:21), replace = TRUE))
               )
length(unique(df$y))
res <- new_mutual_cols(df, max_n = 50000)
res$df_stat  %>% 
  mutate(iteration = ifelse(iteration == 0, 1, iteration)) %>%
  ggplot(aes(x = iteration, y = mutual_info, color = method)) +
  geom_line() +
  labs(title = "Mutual Information Over Iterations", x = "Iteration", y = "Mutual Information") +
  theme_minimal() +
  facet_wrap(~type, scales = "free_y") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 60)) + 
  scale_x_log10()

# Example usage: SMALL
set.seed(42)
df <- data.frame(
  x = sample(str_c("Categ", 1:3), 10000, replace = TRUE)
) %>% mutate(y = 
               ifelse(str_ends(x, "1"), sample(str_c("Categ", 1:2), replace = TRUE),
                      sample(str_c("Categ", 1:3), replace = TRUE))
)
res <- new_mutual_cols(df, max_n = 50000)
res$df_stat  %>% 
  mutate(iteration = ifelse(iteration == 0, 1, iteration)) %>%
  ggplot(aes(x = iteration, y = mutual_info, color = method)) +
  geom_line() +
  labs(title = "Mutual Information Over Iterations", x = "Iteration", y = "Mutual Information") +
  theme_minimal() +
  facet_wrap(~type, scales = "free_y") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 60)) + 
  scale_x_log10()

