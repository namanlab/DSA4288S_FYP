# Load necessary packages
library(MASS)        # For polr() - cumulative logit model
library(nnet)        # For multinom() - multinomial logit model
library(ordinal)     # For clm() - adjacent category logit model
library(tidyverse)   # For data wrangling
library(vcd)         # For generating contingency tables
library(stringr)
library(gridExtra)
library(genodds)
library(DescTools) #


log_odds_dc <- function(tab) {
  res <- ConDisPairs(tab)[c("C", "D")]
  concordant <- res$C + 1/6
  discordant <- res$D + 1/6
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
        odds_ratio <- (a * d) + 1/6
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
  weighted_log_odds <- log((total_weighted_odds + 1/6)/(total_weight + 1/6))
  return(weighted_log_odds)
}

generate_random_mat <- function(n, m, min_val = 1, max_val = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)  # Set seed for reproducibility
  matrix(sample(min_val:max_val, n * m, replace = TRUE), ncol = m, byrow = TRUE)
}



################################################################################
############################# EX 1: TRACK RELATION #############################
################################################################################

# function to generate a random matrix
generate_random_mat <- function(n, m, min_val = 1, max_val = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)  # Set seed for reproducibility
  matrix(sample(min_val:max_val, n * m, replace = TRUE), ncol = m, byrow = TRUE)
}

iterations <- 1000
results <- data.frame(iter = integer(), m = integer(), 
                      log_odds_dc = double(), log_odds_weighted = double())

set.seed(123)
for (i in 1:iterations) {
  m <- sample(3:10, 1)  # not 2*2 because 2*2 will be smae values
  mat <- generate_random_mat(m, m)  
  log_dc <- log_odds_dc(mat) 
  log_weighted <- log_odds_weighted(mat)
  results <- rbind(results, data.frame(iter = i, m = m, log_odds_dc = log_dc, log_odds_weighted = log_weighted))
}
ggplot(results, aes(x = log_odds_dc, y = log_odds_weighted)) +
  geom_point(aes(color = factor(m)), alpha = 0.7, size = 2) +
  labs(
    title = "Comparison of Agresti and Weighted Log-Odds Ratios",
    x = "Agresti Log-Odds",
    y = "Log-Odds (Weighted)",
    color = "Matrix Size (m)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right"
  ) +
  scale_colour_viridis_d() + geom_smooth()



################################################################################
############################# EX2: EFFECT OF AUGMENTATION #############################
################################################################################

# Does augmentation/strengtheing diagonals affect association as intended?

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
    if (iters > 100){
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
set.seed(42)
iterations <- 300
results <- data.frame(iter = integer(), dim_mat = integer(),
                      log_odds_dc = double(), log_odds_weighted = double())
dimensions_exp <- c(2, 3, 7, 12)
for (dim_c in dimensions_exp){
  sample_mat = generate_random_mat(dim_c, dim_c)
  for (i in 1:iterations) {
    sample_mat <<- augment_matrix_random_block(sample_mat, 1)
    log_dc <- log_odds_dc(sample_mat) 
    log_weighted <- log_odds_weighted(sample_mat)
    results <- rbind(results, data.frame(iter = i, dim_mat = dim_c, log_odds_dc = log_dc, log_odds_weighted = log_weighted))
  }
}
plot_dimension <- function(dim) {
  results %>% filter(dim_mat == dim) %>%
    pivot_longer(log_odds_dc:log_odds_weighted, names_to = "metric", values_to = "log_odds") %>%
    mutate(metric = ifelse(metric == "log_odds_dc", "Agresti Log Odds", "Weighted Log Odds")) %>%
    ggplot(aes(x = iter, y = log_odds)) +
    geom_line() + labs(title = paste("Dimension:", dim), x = "Iteration", y = "Log-Odds") +
    facet_wrap(~metric, scales = "free_y") + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
}

plots <- lapply(dimensions_exp, plot_dimension)
grid.arrange(grobs = plots, ncol = 2)




################################################################################
############################# EX 3: MODEL FITTINGS #############################
################################################################################


fit_models <- function(table) {
  # Convert table to data frame format
  df <- as.data.frame(as.table(table))
  colnames(df) <- c("Row", "Col", "Freq")
  
  # Model 1: Cumulative Logit (Proportional Odds)
  model_cumlogit <- polr(as.factor(Row) ~ as.factor(Col), weights = Freq, data = df, method = "logistic")
  coef_cumlogit <- coef(model_cumlogit)
  
  # Model 2: Multinomial Logit
  model_multinom <- multinom(as.factor(Row) ~ as.factor(Col), weights = Freq, data = df)
  coef_multinom <- coef(model_multinom)
  
  # Model 3: Adjacent Category Logit Model
  model_adjcat <- clm(as.factor(Row) ~ as.factor(Col), weights = Freq, data = df)
  coef_adjcat <- coef(model_adjcat)
  
  # Return all coefficients as a named list
  return(list(
    cumlogit = coef_cumlogit,
    multinom = coef_multinom,
    adjcat = coef_adjcat
  ))
}

r2_summary <- tibble(size = integer(), model = character(), r2_nc = double(), 
                     r2adj_nc = double(), r2_w = double(), r2adj_w = double())
set.seed(123)
n_iter = 500

# 5. Loop through matrix sizes
for (size in 3:8) {
  
  # Initialize coefficient tibbles for this size
  cumlogit_tibble <<- list()
  multinom_tibble <<- list()
  adjcat_tibble <<- list()
  log_odds_nc_vec <<- rep(0, n_iter)
  log_odds_w_vec <<- rep(0, n_iter)
  
  # Run iterations for the current size
  for (i in 1:n_iter) {
    table <- generate_random_mat(size, size)
    log_odds_nc <- log_odds_dc(table)
    log_odds_nc_vec[i] = log_odds_nc
    log_odds_w <- log_odds_weighted(table)
    log_odds_w_vec[i] = log_odds_w
    model_coefs <- fit_models(table)
    cumlogit_tibble[[i]] <- as_tibble_row(model_coefs$cumlogit, .name_repair = "unique")
    multinom_tibble[[i]] <- as_tibble_row(model_coefs$multinom, .name_repair = "unique")
    adjcat_tibble[[i]] <- as_tibble_row(model_coefs$adjcat, .name_repair = "unique")
  }
  
  # Combine all iterations into one data frame for each model
  cumlogit_df <<- bind_rows(cumlogit_tibble)
  multinom_df <<- bind_rows(multinom_tibble) 
  adjcat_df <<- bind_rows(adjcat_tibble)
  
  # Helper function to calculate R² and adjusted R²
  calculate_r2 <- function(df, log_odds_col) {
    fit <- lm(log_odds_col ~ ., data = df)
    list(r2 = summary(fit)$r.squared, adj_r2 = summary(fit)$adj.r.squared)
  }
  
  # Calculate R² and adjusted R² for each model
  cumlogit_r2nc <- calculate_r2(cumlogit_df , log_odds_nc_vec)
  multinom_r2nc <- calculate_r2(multinom_df, log_odds_nc_vec)
  adjcat_r2nc <- calculate_r2(adjcat_df, log_odds_nc_vec)
  
  cumlogit_r2w <- calculate_r2(cumlogit_df , log_odds_w_vec)
  multinom_r2w <- calculate_r2(multinom_df, log_odds_w_vec)
  adjcat_r2w <- calculate_r2(adjcat_df, log_odds_w_vec)
  
  # Store R² results for each model
  r2_summary <- r2_summary %>%
    add_row(size = size, model = "cumlogit", r2_nc = cumlogit_r2nc$r2, r2adj_nc = cumlogit_r2nc$adj_r2,
            r2_w = cumlogit_r2w$r2, r2adj_w = cumlogit_r2w$adj_r2) %>%
    add_row(size = size, model = "multinom", r2_nc = multinom_r2nc$r2, r2adj_nc = multinom_r2nc$adj_r2,
            r2_w = multinom_r2w$r2, r2adj_w = multinom_r2w$adj_r2) %>%
    add_row(size = size, model = "adjcat", r2_nc = adjcat_r2nc$r2, r2adj_nc = adjcat_r2nc$adj_r2,
            r2_w = adjcat_r2w$r2, r2adj_w = adjcat_r2w$adj_r2)
}

# 6. Plot the R² and adjusted R² metrics
r2_long <- r2_summary %>%
  pivot_longer(cols = starts_with("r2"), names_to = c("metric", "log_odds_type"), names_sep = "_") %>%
  mutate(metric = ifelse(metric == "r2", "R2", "Adj R2"),
         model = ifelse(model == "cumlogit", "Cumulative Logit", 
                        ifelse(model == "multinom", "Multinomial Logit", "Adjacent Category Logit")),
         log_odds_type = ifelse(log_odds_type == "w", "Local Weighted", "Agresti Generalized"))

ggplot(r2_long, aes(x = as.factor(size), y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(log_odds_type~ model) +
  labs(
    title = "R² and Adjusted R² Comparison across Matrix Sizes",
    x = "Matrix Size",
    y = "R² / Adjusted R²",
    fill = "Metric"
  ) +
  theme_bw() +
  scale_fill_manual(values = c('red', "steelblue"))
