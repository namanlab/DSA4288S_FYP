library(tidyverse)
library(stringr)
library(gridExtra)
# install.packages("genodds")
library(genodds)
# install.packages("DescTools")
library(DescTools) # concord and siscord pairs

################################################################################
############################# Manual Specification #############################
################################################################################


log_odds_dc <- function(tab) {
  res <- ConDisPairs(tab)[c("C", "D")]
  concordant <- res$C + 1
  discordant <- res$D + 1
  if (discordant == 0) return(Inf)  # Handle division by zero
  return(log(concordant / discordant))
}

log_odds_weighted <- function(mat) {
  n <- nrow(mat)
  m <- ncol(mat)
  total_weighted_odds <- 0
  total_weight <- 0
  for (i in 1:(n - 1)) {
    for (j in 1:(m - 1)) {
      sub_mat <- mat[i:(i + 1), j:(j + 1)]  
      
      a <- sub_mat[1, 1]
      b <- sub_mat[1, 2]
      c <- sub_mat[2, 1]
      d <- sub_mat[2, 2]
      
      # Calculate odds ratio, handling zero division
      if (b * c == 0) {
        odds_ratio <- (a * d) + 1
      } else {
        odds_ratio <- (a * d) / (b * c) 
      }
      
      # calculate weight (total observations in the 2x2 submatrix)
      weight <- a + b + c + d
      
      # accumulate weighted sum of odds ratios
      total_weighted_odds <- total_weighted_odds + (weight * odds_ratio)
      total_weight <- total_weight + weight
    }
  }
  
  # Calculate final weighted average of log-odds ratios
  weighted_log_odds <- log((total_weighted_odds + 1)/(total_weight + 1))
  return(weighted_log_odds)
}





############################# SL Based Procedure #############################

plot_log_odds <- function(matrices, result, names_matrices, log_odds_general = log_odds_dc) {
  before_log_odds <- sapply(matrices, log_odds_general)
  after_log_odds <- sapply(result[[1]], log_odds_general)
  names_matrices <- c(names_matrices, "Overall")
  overall_before <- log_odds_general(Reduce("+", matrices))
  overall_after <- log_odds_general(Reduce("+", result[[1]]))
  before_log_odds <- c(before_log_odds, overall_before)
  after_log_odds <- c(after_log_odds, overall_after)
  data <- data.frame(matrix_name = c(names_matrices, names_matrices),
                     log_odds = c(before_log_odds, after_log_odds),
                     type = c(rep("Before", length(names_matrices)), rep("After", length(names_matrices))))
  data %>%
    mutate(type = factor(type, levels = c("Before", "After"), ordered = TRUE)) %>%
    ggplot(aes(x = matrix_name, y = log_odds, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    xlab("Matrix") +
    ylab("Log Odds") +
    ggtitle("Log Odds Before and After Transformation")
}


softmax <- function(x) {
  # If any value is Inf, assign it a probability of 1, others 0
  if (any(is.infinite(x) & x > 0)) {
    probs <- rep(0, length(x))
    probs[which(x == Inf)] <- 1
    return(probs)
  }
  # Handle the case where all values are -Inf
  if (all(x == -Inf)) {
    return(rep(0, length(x)))
  }
  # Standard softmax calculation with numerical stability
  exp_x <- exp(x - max(x, na.rm = TRUE))  # Stability adjustment
  res = exp_x / sum(exp_x, na.rm = TRUE)
  return(res/sum(res))
}

augment_matrix_random_block <- function(table, delta) {
  nrows <- nrow(table)
  ncols <- ncol(table)
  iters <- 0
  repeat {
    iters <- iters + 1
    RR <- sort(sample(1:nrows, 2))
    CC <- sort(sample(1:ncols, 2))
    if (delta < 0 && table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
    if (delta > 0 && table[RR[1], CC[2]] > 0 && table[RR[2], CC[1]] > 0) break
    if (iters > 100) {
      delta <- 0
      break
    }
  }
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] + delta
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] + delta
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] - delta
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] - delta
  return(table)
}


adjust_matrices <- function(matrices, 
                            manual_vec, target_overall, 
                            margin, margin_overall, 
                            max_n = 1000, max_n_overall = 1000, temp = 10,
                            log_odds_general = log_odds_dc) {
  
  sum_matrices <- function(mats) Reduce("+", mats)
  curr <- matrices
  curr_eval_sum <- log_odds_general(sum_matrices(curr))
  history <- data.frame(iteration = 0, overall_log_odds = curr_eval_sum)
  for (n in 1:max_n) {
    log_odds_values <- sapply(curr, log_odds_general)
    selection_probs <- sapply(1:length(curr), function(i) {
      cur_log_odds <- log_odds_values[i]
      target <- manual_vec[i]
      if (cur_log_odds < -margin && target < 0) return(-Inf)  # No modification
      if (cur_log_odds > margin && target > 0) return(-Inf)  # No modification
      return(abs(cur_log_odds - target))
    })
    # print(selection_probs)
    if (all(is.infinite(selection_probs))) {selection_probs <- rep(0, length(selection_probs))}
    selection_probs <- softmax(selection_probs)  # Normalize to softmax probabilities
    idx <- sample(1:length(curr), 1, prob = selection_probs)
    cur_log_odds <- log_odds_general(curr[[idx]])
    target <- manual_vec[idx]
    delta <- 1
    if (cur_log_odds - target*margin > 0){delta = -1}
    modified_matrix <- augment_matrix_random_block(curr[[idx]], delta)
    new_matrices <- curr
    new_matrices[[idx]] <- modified_matrix
    new_eval_sum <- log_odds_general(sum_matrices(new_matrices))
    temp <- 0.9 * temp  # Annealing step
    metropolis <- ifelse(delta <= 0, exp((curr_eval_sum - new_eval_sum) / temp),
                         exp((new_eval_sum - curr_eval_sum)/temp))
    # Metropolis criteria for acceptance and constraints on log-odds limits
    if (runif(1) < metropolis) {
      curr <- new_matrices
      curr_eval_sum <- new_eval_sum
    }
    history <- rbind(history, data.frame(iteration = n, overall_log_odds = curr_eval_sum))
  }
  
  ## NOW ADJUST OVERALL
  curr_eval_sum <- log_odds_general(sum_matrices(curr))
  if (curr_eval_sum < margin_overall*target_overall && target_overall > 0){
    # increase
    for (n in 1:max_n_overall) {
      log_odds_values <- sapply(curr, log_odds_general)
      selection_probs <- sapply(1:length(curr), function(i) {
        cur_log_odds <- log_odds_values[i]
        target <- manual_vec[i]
        if (target > 0){return(1000*exp(-cur_log_odds))}
        else {return(-margin - cur_log_odds)}
      })
      if (all(is.infinite(selection_probs))) {selection_probs <- rep(0, length(selection_probs))}
      selection_probs <- softmax(selection_probs)  # Normalize to softmax probabilities
      idx <- sample(1:length(curr), 1, prob = selection_probs)
      delta <- 1
      modified_matrix <- augment_matrix_random_block(curr[[idx]], delta)
      new_matrices <- curr
      new_matrices[[idx]] <- modified_matrix
      new_eval_sum <- log_odds_general(sum_matrices(new_matrices))
      temp <- 0.9 * temp  # Annealing step
      metropolis <- ifelse(delta <= 0, exp((curr_eval_sum - new_eval_sum) / temp),
                           exp((new_eval_sum - curr_eval_sum)/temp))
      # CHECK IF MODIFIED MATRIX CONSTRAINTS ARE NOT VIOLATED!
      cur_log_odds <- log_odds_general(modified_matrix)
      target <- manual_vec[idx]
      CHECK_F = 0
      if (cur_log_odds < -margin && target < 0){CHECK_F = 1} 
      if (cur_log_odds > margin && target > 0){CHECK_F = 1} 
      # Metropolis criteria for acceptance and constraints on log-odds limits
      if (runif(1) < metropolis && CHECK_F == 1) {
        curr <- new_matrices
        curr_eval_sum <- new_eval_sum
      }
      history <- rbind(history, data.frame(iteration = n + max_n, overall_log_odds = curr_eval_sum))
    }
  } else if (curr_eval_sum > margin_overall*target_overall && target_overall < 0){
    # decrease
    for (n in 1:max_n_overall) {
      log_odds_values <- sapply(curr, log_odds_general)
      selection_probs <- sapply(1:length(curr), function(i) {
        cur_log_odds <- log_odds_values[i]
        target <- manual_vec[i]
        if (target < 0){return(1000*exp(cur_log_odds))}
        else {return(cur_log_odds - margin)}
      })
      if (all(is.infinite(selection_probs))) {selection_probs <- rep(0, length(selection_probs))}
      selection_probs <- softmax(selection_probs)  # Normalize to softmax probabilities
      idx <- sample(1:length(curr), 1, prob = selection_probs)
      delta <- -1
      modified_matrix <- augment_matrix_random_block(curr[[idx]], delta)
      new_matrices <- curr
      new_matrices[[idx]] <- modified_matrix
      new_eval_sum <- log_odds_general(sum_matrices(new_matrices))
      temp <- 0.9 * temp  # Annealing step
      metropolis <- ifelse(delta <= 0, exp((curr_eval_sum - new_eval_sum) / temp),
                           exp((new_eval_sum - curr_eval_sum)/temp))
      # CHECK IF MODIFIED MATRIX CONSTRAINTS ARE NOT VIOLATED!
      cur_log_odds <- log_odds_general(modified_matrix)
      target <- manual_vec[idx]
      CHECK_F = 0
      if (cur_log_odds < -margin && target < 0){CHECK_F = 1} 
      if (cur_log_odds > margin && target > 0){CHECK_F = 1} 
      # Metropolis criteria for acceptance and constraints on log-odds limits
      if (runif(1) < metropolis && CHECK_F == 1) {
        curr <- new_matrices
        curr_eval_sum <- new_eval_sum
      }
      history <- rbind(history, data.frame(iteration = n + max_n, overall_log_odds = curr_eval_sum))
    }
  }
  # Return the final matrices and log odds history
  return(list(curr, curr_eval_sum, history))
}


################################## EXAMPLE 1 ###################################

set.seed(42)
matrices <- list(
  ta = matrix(c(512, 89, 313, 19), ncol = 2, byrow = TRUE),
  tb = matrix(c(353, 17, 207, 8), ncol = 2, byrow = TRUE),
  tc = matrix(c(120, 202, 205, 391), ncol = 2, byrow = TRUE),
  td = matrix(c(138, 131, 279, 244), ncol = 2, byrow = TRUE),
  te = matrix(c(53, 94, 138, 299), ncol = 2, byrow = TRUE),
  tf = matrix(c(22, 24, 351, 317), ncol = 2, byrow = TRUE)
)
result <- adjust_matrices(matrices, 
                          manual_vec = c(-1, -1, -1, -1, -1, -1), 
                          target_overall = +1, 
                          margin = 0.2,  margin_overall = 0.2,
                          max_n = 200, max_n_overall = 100)
plot_log_odds(matrices, result, names(matrices))
result[[3]] %>%
  ggplot(aes(x = iteration, y = overall_log_odds)) +
  geom_line() +
  theme_minimal() +
  xlab("Iteration") +
  ylab("Overall Log Odds") +
  ggtitle("Log Odds Trajectory During Adjustment")
# USING WEIGHTED LOG ODDS:
# result <- adjust_matrices(matrices, 
#                           manual_vec = c(-1, -1, -1, -1, -1, -1), 
#                           target_overall = +1, 
#                           margin = 0.2,  margin_overall = 0.2,
#                           max_n = 200, max_n_overall = 100,
#                           log_odds_general = log_odds_weighted)
# plot_log_odds(matrices, result, names(matrices),
#               log_odds_general = log_odds_weighted)
# result[[3]] %>%
#   ggplot(aes(x = iteration, y = overall_log_odds)) +
#   geom_line() +
#   theme_minimal() +
#   xlab("Iteration") +
#   ylab("Overall Log Odds") +
#   ggtitle("Log Odds Trajectory During Adjustment")



################################## EXAMPLE 2 ###################################

set.seed(42)
matrices <- list(
  ta = matrix(c(512, 89, 313, 19, 45, 67), ncol = 3, byrow = TRUE),
  tb = matrix(c(353, 17, 207, 8, 33, 25), ncol = 3, byrow = TRUE),
  tc = matrix(c(120, 202, 205, 391, 55, 23), ncol = 3, byrow = TRUE),
  td = matrix(c(138, 131, 279, 244, 111, 99), ncol = 3, byrow = TRUE),
  te = matrix(c(53, 94, 138, 299, 87, 77), ncol = 3, byrow = TRUE),
  tf = matrix(c(22, 24, 351, 317, 11, 33), ncol = 3, byrow = TRUE),
  tg = matrix(c(22, 24, 351, 317, 11, 33), ncol = 3, byrow = TRUE)
)
result <- adjust_matrices(matrices, 
                          manual_vec = c(-1, -1, -1, 1, -1, -1, -1), 
                          target_overall = +1, 
                          margin = 0.1,  margin_overall = 0.1,
                          max_n = 2000, max_n_overall = 1000)
plot_log_odds(matrices, result, names(matrices))
result[[3]] %>%
  ggplot(aes(x = iteration, y = overall_log_odds)) +
  geom_line() +
  theme_minimal() +
  xlab("Iteration") +
  ylab("Overall Log Odds") +
  ggtitle("Log Odds Trajectory During Adjustment")
# USING WEIGHTED LOG ODDS:
# result <- adjust_matrices(matrices, 
#                           manual_vec = c(-1, -1, -1, -1, -1, -1, -1), 
#                           target_overall = +1, 
#                           margin = 0.1,  margin_overall = 0.1,
#                           max_n = 2000, max_n_overall = 2000,
#                           log_odds_general = log_odds_weighted)
# plot_log_odds(matrices, result, names(matrices),
#               log_odds_general = log_odds_weighted)
# result[[3]] %>%
#   ggplot(aes(x = iteration, y = overall_log_odds)) +
#   geom_line() +
#   theme_minimal() +
#   xlab("Iteration") +
#   ylab("Overall Log Odds") +
#   ggtitle("Log Odds Trajectory During Adjustment")


############################# EXAMPLE 3: steroids :) ############################


generate_random_matrices <- function(p, n, m, min_val = 1, max_val = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)  # Set seed for reproducibility
  # use lapply to generate a list of 'p' matrices
  matrices <- lapply(1:p, function(i) {
    x = matrix(sample(min_val:max_val, n * m, replace = TRUE), ncol = m, byrow = TRUE)
    if (i > p/4){
      for (j in round(m/2):m){
        x[, j] = round(x[, j]/10)
      }
    }
    x
  })
  # Name each matrix as "t1", "t2", ..., "tp"
  names(matrices) <- paste0("t", 1:p)
  return(matrices)
}

set.seed(21) 
matrices <- generate_random_matrices(p = 8, n = 4, m = 4)
result <- adjust_matrices(matrices, 
                          manual_vec = c(-1, -1, -1, -1, -1, -1, -1, -1), 
                          target_overall = +1, 
                          margin = 0.05,  margin_overall = 0.1,
                          max_n = 2000, max_n_overall = 1000)
plot_log_odds(matrices, result, names(matrices))
result[[3]] %>%
  ggplot(aes(x = iteration, y = overall_log_odds)) +
  geom_line() +
  theme_minimal() +
  xlab("Iteration") +
  ylab("Overall Log Odds") +
  ggtitle("Log Odds Trajectory During Adjustment")
# USING WEIGHTED LOG ODDS:
# result <- adjust_matrices(matrices, 
#                           manual_vec = c(-1, -1, -1, -1, -1, -1, -1, -1), 
#                           target_overall = +1, 
#                           margin = 0.1,  margin_overall = 0.1,
#                           max_n = 2000, max_n_overall = 1000,
#                           log_odds_general = log_odds_weighted)
# plot_log_odds(matrices, result, names(matrices),
#               log_odds_general = log_odds_weighted)
# result[[3]] %>%
#   ggplot(aes(x = iteration, y = overall_log_odds)) +
#   geom_line() +
#   theme_minimal() +
#   xlab("Iteration") +
#   ylab("Overall Log Odds") +
#   ggtitle("Log Odds Trajectory During Adjustment")













