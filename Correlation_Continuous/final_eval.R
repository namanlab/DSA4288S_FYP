# Load necessary libraries
library(MCMCpack)
library(GGally)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(akima)
library(actuar)
library(copula)
library(goftest)
options(scipen = 999)


################################################################################
################################ COPULAS FUNC ##################################
################################################################################


# GAUSSIAN
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

# STUDENTS'S T
generate_t_copula_samples <- function(n, d, rho_matrix, df) {
  # Step 1: Generate multivariate t samples
  mvn_samples <- rmvt(n = n, sigma = rho_matrix, df = df)
  
  # Step 2: Transform samples to uniform using the CDF of the t-distribution
  uniform_samples <- pt(mvn_samples, df = df)
  
  return(uniform_samples)
}
t_copula_2 <- function(n, p) {
  rho_matrix <- matrix(c(1, p, p, 1), nrow = 2, ncol = 2)
  smpls <- generate_t_copula_samples(n, 2, rho_matrix, 5)
  return(smpls)
}


# CLAYTON
generate_clayton_copula_samples <- function(n, d, theta) {
  rCopula(n = n, copula = claytonCopula(param = theta, dim = d))
}
clayton_copula_2 <- function(n, p){
  theta <- 2*p/(1 - p) # kendall tau is theta/(2 + theta)
  generate_clayton_copula_samples(n, 2, theta)
}

# GUMBEL
generate_gumbel_copula_samples <- function(n, d, theta) {
  rCopula(n = n, copula = gumbelCopula(param = theta, dim = d))
}
gumbel_copula_2 <- function(n, p){
  if (p == 1){return(NA)}
  if (p < 0){return(NA)}
  theta <- 1/(1 - p)# kendall tau is 1 - 1/theta
  generate_clayton_copula_samples(n, 2, theta)
}

# AMH
generate_amh_copula_samples <- function(n, d, theta) {
  rCopula(n = n, copula = amhCopula(param = theta, dim = d))
}
tau_function_amh  <- function(theta) {
  (3 * theta - 2) / (3 * theta) - (2 * (1 - theta)^2 * log(1 - theta)) / (3 * theta^2)
}
inverse_tau_function_amh <- function(tau_value, tol = 1e-6) {
  # Use uniroot to find the theta that corresponds to the given tau
  result <- uniroot(function(theta) tau_function_amh(theta) - tau_value, 
                    lower = -1 + tol, upper = 1 - tol)
  return(result$root)
}
amh_copula_2 <- function(n, p){
  if (p <= -0.1817){return(NA)}
  if (p >= 0.33333){return(NA)}
  theta <- inverse_tau_function_amh(p)
  generate_amh_copula_samples(n, 2, theta)
}



################################################################################
################################ INV CDF FUNC ##################################
################################################################################

# Inverse CDF generator using quantiles
genCDFInv_quantile <- function(data) {
  return(function(p) {quantile(data, p)})
}

# Inverse CDF generator using linear interpolation
genCDFInv_linear <- function(X) {
  ecdf1 <- ecdf(X)
  U <- sort(unique(ecdf1(X)))
  Finv <- sort(X)
  ak2 <- approxfun(U, Finv, method = "linear", rule = 2)
  return(ak2)
}

# Inverse CDF generator using Akima spline interpolation
genCDFInv_akima <- function(X) {
  ecdf1 <- ecdf(X)
  U <- sort(unique(ecdf1(X)))
  Finv <- sort(X)
  ak2 <- function(x) {aspline(U, Finv, x)}
  return(ak2)
}

# Inverse CDF generator using polynomial regression
genCDFInv_poly <- function(data, degree = 10) {
  x <- sort(unique(data))
  y <- ecdf(data)(x)
  poly_fit <- lm(x ~ poly(y, degree = degree, raw = TRUE))
  cdf_poly <- function(y) predict(poly_fit, newdata = data.frame(y = y))
  return(cdf_poly)
}

# ADD BAYES


################################################################################
################################# EVAL FUNC ####################################
################################################################################


# 1. Function to Evaluate Performance of Inverse CDF Estimators
evaluate_inv_cdf <- function(true_gen_fn, est_inv_cdf, N_test) {
  x_test <- true_gen_fn(N_test)
  x_vals <- seq(0, 1, length.out = N_test)
  samples <- quantile(x_test, x_vals)
  estimated <- est_inv_cdf(x_vals)
  mse <- mean((samples - estimated)^2)
  mae <- mean(abs(samples - estimated))
  list(MSE = mse, MAE = mae)
}


# 2. Functions to calculate Total Variation (TV) distance
calculate_tv_distance_knownfn <- function(true_gen_fn, generated_data, N_test) {
  original_data <- true_gen_fn(N_test)
  return(calculate_tv_distance_empirical(original_data, generated_data))
}
calculate_tv_distance_empirical <- function(original_data, generated_data) {
  # Create empirical CDFs for original and generated data
  original_ecdf <- ecdf(original_data)
  generated_ecdf <- ecdf(generated_data)
  # Define the grid over which to calculate the distance
  x_values <- sort(unique(c(original_data, generated_data)))
  # Calculate the TV distance
  tv_distance <- max(abs(original_ecdf(x_values) - generated_ecdf(x_values)))
  list(tv_distance = tv_distance)
}


# 3. Functions to Perform Kolmogorov-Smirnov Test
perform_ks_test_knownfn <- function(true_gen_fn, generated_data, N_test) {
  original_data <- true_gen_fn(N_test)
  return(perform_ks_test_empirical(original_data, generated_data))
}
perform_ks_test_empirical <- function(original_data, generated_data) {
  ks_test_result <- ks.test(original_data, generated_data)
  list(ks_test_result = ks_test_result)
}

# 4. Functions to Perform Goodness of Fit: Cramer-von Mises Test
perform_cvm_test_knownfn <- function(true_gen_fn, generated_data, N_test) {
  x_test <- true_gen_fn(N_test)
  Fn = ecdf(x_test)
  cvm_test <- cvm.test(generated_data, Fn)
  list(CramerVonMises = cvm_test)
}
perform_cvm_test_empirical <- function(original_data, generated_data) {
  Fn = ecdf(original_data)
  cvm_test <- cvm.test(generated_data, Fn)
  list(CramerVonMises = cvm_test)
}

# 5. Function to Evaluate Correlation Performance
evaluate_performance_corr <- function(x1, x2, target_corr_kendall) {
  kendall_tau <- cor(x1, x2, method = "kendall")
  err_kendall <- kendall_tau - target_corr_kendall
  return(list(err_kendall = err_kendall))
}



################################################################################
################################################################################
################################ FINAL EVAL FN #################################
################################################################################
################################################################################


################################################################################
############################# KNOWN GENERATING FN ##############################
################################################################################

evaluate_copulas_and_inverse_cdfs_parametric <- function(target_corr_kendall,
                                              distribution_1, distribution_2, copula_type, 
                                              inv_cdf_type, n = 1000, n_test = 100000) {
  
  # Generate original data
  d1 <- distribution_1(n)
  d2 <- distribution_2(n)
  
  # Generate copula samples
  if (copula_type == "gaussian") {
    target_rho = sin(target_corr_kendall*pi/2) # Greiner's equality 
    res <- gauss_copula_2(n, target_rho) 
  } else if (copula_type == "t") {
    target_rho = sin(target_corr_kendall*pi/2) # Greiner's equality 
    res <- t_copula_2(n, target_rho)
  } else if (copula_type == "clayton") {
    res <- clayton_copula_2(n, target_corr_kendall)
  } else if (copula_type == "gumbel") {
    res <- gumbel_copula_2(n, target_corr_kendall)
  } else if (copula_type == "amh") {
    res <- amh_copula_2(n, target_corr_kendall)
  } else {
    stop("Unsupported copula type")
  }
  
  if (all(is.na(res))){
    results <- list(
      inv_cdf_eval_1 = list(NA),
      inv_cdf_eval_2 = list(NA),
      tv_result_1 = list(NA),
      tv_result_2 = list(NA),
      ks_test_result_1 = list(NA),
      ks_test_result_2 = list(NA),
      cvm_test_result_1 = list(NA),
      cvm_test_result_2 = list(NA),
      corr_change_test = list(NA)
    )
    return(results)
  }
  
  # Apply chosen inverse CDF
  inv_cdf_d1 <- switch(inv_cdf_type,
                       "quantile" = genCDFInv_quantile(d1),
                       "linear" = genCDFInv_linear(d1),
                       "akima" = genCDFInv_akima(d1),
                       "poly" = genCDFInv_poly(d1))
  
  inv_cdf_d2 <- switch(inv_cdf_type,
                       "quantile" = genCDFInv_quantile(d2),
                       "linear" = genCDFInv_linear(d2),
                       "akima" = genCDFInv_akima(d2),
                       "poly" = genCDFInv_poly(d2))
  
  F1Inv <- Vectorize(inv_cdf_d1)
  F2Inv <- Vectorize(inv_cdf_d2)
  
  x1 <- F1Inv(res[,1])
  x2 <- F2Inv(res[,2])
  
  # Evaluate performance of inverse CDFs
  inv_cdf_eval_1 <- evaluate_inv_cdf(distribution_1, inv_cdf_d1, n_test)
  inv_cdf_eval_2 <- evaluate_inv_cdf(distribution_2, inv_cdf_d2, n_test)
  
  # Calculate TV Distance
  tv_result_1 <- calculate_tv_distance_knownfn(distribution_1, x1, n_test)
  tv_result_2 <- calculate_tv_distance_knownfn(distribution_2, x2, n_test)
  
  # Perform Kolmogorov-Smirnov Test
  ks_test_result_1 <- perform_ks_test_knownfn(distribution_1, x1, n_test)
  ks_test_result_2 <- perform_ks_test_knownfn(distribution_2, x2, n_test)
  
  # Perform Cramer-von Mises Test
  cvm_test_result_1 <- perform_cvm_test_knownfn(distribution_1, x1, n_test)
  cvm_test_result_2 <- perform_cvm_test_knownfn(distribution_2, x2, n_test)
  
  # Assess correlation changes
  corr_change_test <- evaluate_performance_corr(x1, x2, target_corr_kendall)
  
  # Return all results
  results <- list(
    inv_cdf_eval_1 = inv_cdf_eval_1,
    inv_cdf_eval_2 = inv_cdf_eval_2,
    tv_result_1 = tv_result_1,
    tv_result_2 = tv_result_2,
    ks_test_result_1 = ks_test_result_1,
    ks_test_result_2 = ks_test_result_2,
    cvm_test_result_1 = cvm_test_result_1,
    cvm_test_result_2 = cvm_test_result_2,
    corr_change_test = corr_change_test
  )
  
  return(results)
}


set.seed(123)
distribution_normal <- function(n) rnorm(n, mean = 0, sd = 1)
distribution_exponential <- function(n) rexp(n, rate = 1)
distribution_uniform <- function(n) runif(n, min = 0, max = 1)
distribution_log_normal <- function(n) rlnorm(n, meanlog = 0, sdlog = 1)
distribution_beta <- function(n) rbeta(n, shape1 = 2, shape2 = 5)
distribution_t <- function(n) rt(n, df = 2)
distribution_mixture1 <- function(n) {
  prob <- runif(n)
  rnorm(n, mean = ifelse(prob > 0.5, 3, -2), sd = 1)
}
distribution_mixture2 <- function(n) {
  prob <- runif(n)
  ifelse(prob > 0.5, rnorm(n, mean = 0, sd = 1), rlnorm(n, meanlog = 0, sdlog = 0.5))
}

# PARAM SPACE
target_corr_kendall <- seq(-1, 1, 0.1)
distributions <- list(
  normal = distribution_normal,
  exponential = distribution_exponential,
  uniform = distribution_uniform,
  log_normal = distribution_log_normal,
  beta = distribution_beta,
  student_t = distribution_t,
  mixture1 = distribution_mixture1,
  mixture2 = distribution_mixture2
)
inv_cdf_types <- c("quantile", "linear", "poly")
copula_types <- c("gaussian", "t")
n_range <- c(100, 1000, 10000)

results <- NULL
iter = 0
for (n_val in n_range){
  for (tau in target_corr_kendall) {
    for (dist_name_1 in names(distributions)) {
      for (dist_name_2 in names(distributions)) {
        for (inv_cdf_type in inv_cdf_types) {
          for (copula_type in copula_types) {
            
            result <- evaluate_copulas_and_inverse_cdfs_parametric(
              target_corr_kendall = tau,
              distribution_1 = distributions[[dist_name_1]],
              distribution_2 = distributions[[dist_name_2]],
              copula_type = copula_type,
              inv_cdf_type = inv_cdf_type,
              n = n_val
            )
            print(iter)
            iter = iter + 1
            # Store results in a data frame
            temp_result <- data.frame(
              n = n_val,
              target_corr_kendall = tau,
              distribution_1 = dist_name_1,
              distribution_2 = dist_name_2,
              inv_cdf_type = inv_cdf_type,
              copula_type = copula_type,
              mse_inv_cdf_1 = result$inv_cdf_eval_1$MSE,
              mae_inv_cdf_1 = result$inv_cdf_eval_1$MAE,
              mse_inv_cdf_2 = result$inv_cdf_eval_2$MSE,
              mae_inv_cdf_2 = result$inv_cdf_eval_2$MAE,
              tv_val_1 = result$tv_result_1$tv_distance,
              tv_val_2 = result$tv_result_2$tv_distance,
              ks_stat_1 = result$ks_test_result_1$ks_test_result$statistic,
              ks_stat_2 = result$ks_test_result_2$ks_test_result$statistic,
              cvm_stat_1 = result$cvm_test_result_1$CramerVonMises$statistic,
              cvm_stat_2 = result$cvm_test_result_2$CramerVonMises$statistic,
              err_kendall = result$corr_change_test$err_kendall
            )
            results <- bind_rows(results, temp_result)
          }
        }
      }
    }
  }
}



# Save the results
write.csv(results, "results/final_eval_known_fn.csv")





################################ AGAINST COPULA ################################ 

# MSE
results %>% group_by(copula_type, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_mse_inv_cdf)) +
  facet_grid(copula_type~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# MAE
results %>% group_by(copula_type, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_mae_inv_cdf)) +
  facet_grid(copula_type~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# TV DIST
results %>% group_by(copula_type, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_tv)) +
  facet_grid(copula_type~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# KS STAT
results %>% group_by(copula_type, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_ks_stat)) +
  facet_grid(copula_type~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# CVM STAT
results %>% group_by(copula_type, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_cvm_stat)) +
  facet_grid(copula_type~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")



################################## AGAINST N ################################### 

# MSE
results %>% group_by(n, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_mse_inv_cdf)) +
  facet_grid(n~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# MAE
results %>% group_by(n, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_mae_inv_cdf)) +
  facet_grid(n~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# TV DIST
results %>% group_by(n, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_tv)) +
  facet_grid(n~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# KS STAT
results %>% group_by(n, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_ks_stat)) +
  facet_grid(n~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# CVM STAT
results %>% group_by(n, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_tv = mean(c(tv_val_1, tv_val_2)),
            mean_ks_stat = mean(c(ks_stat_1, ks_stat_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_cvm_stat)) +
  facet_grid(n~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")













################################################################################
############################# ACTUAL DATA BOOTSTRAP ############################
################################################################################







evaluate_copulas_and_inverse_cdfs_bootstrap <- function(target_corr_kendall,
                                                        sampled_x1, sampled_x2, 
                                                        act_x1, act_x2, 
                                                        copula_type, inv_cdf_type) {
  
  n = length(sampled_x1)
  
  # Generate copula samples
  if (copula_type == "gaussian") {
    target_rho = sin(target_corr_kendall*pi/2) # Greiner's equality 
    res <- gauss_copula_2(n, target_rho) 
  } else if (copula_type == "t") {
    target_rho = sin(target_corr_kendall*pi/2) # Greiner's equality 
    res <- t_copula_2(n, target_rho)
  } else if (copula_type == "clayton") {
    res <- clayton_copula_2(n, target_corr_kendall)
  } else if (copula_type == "gumbel") {
    res <- gumbel_copula_2(n, target_corr_kendall)
  } else if (copula_type == "amh") {
    res <- amh_copula_2(n, target_corr_kendall)
  } else {
    stop("Unsupported copula type")
  }
  
  if (all(is.na(res))){
    results <- list(
      tv_result_1 = list(NA),
      tv_result_2 = list(NA),
      ks_test_result_1 = list(NA),
      ks_test_result_2 = list(NA),
      cvm_test_result_1 = list(NA),
      cvm_test_result_2 = list(NA),
      corr_change_test = list(NA)
    )
    return(results)
  }
  
  # Apply chosen inverse CDF
  inv_cdf_d1 <- switch(inv_cdf_type,
                       "quantile" = genCDFInv_quantile(sampled_x1),
                       "linear" = genCDFInv_linear(sampled_x1),
                       "akima" = genCDFInv_akima(sampled_x1),
                       "poly" = genCDFInv_poly(sampled_x1))
  
  inv_cdf_d2 <- switch(inv_cdf_type,
                       "quantile" = genCDFInv_quantile(sampled_x2),
                       "linear" = genCDFInv_linear(sampled_x2),
                       "akima" = genCDFInv_akima(sampled_x2),
                       "poly" = genCDFInv_poly(sampled_x2))
  
  F1Inv <- Vectorize(inv_cdf_d1)
  F2Inv <- Vectorize(inv_cdf_d2)
  
  x1 <- F1Inv(res[,1])
  x2 <- F2Inv(res[,2])
  
  # Calculate TV Distance
  tv_result_1 <- calculate_tv_distance_empirical(act_x1, x1)
  tv_result_2 <- calculate_tv_distance_empirical(act_x2, x2)
  
  # Perform Kolmogorov-Smirnov Test
  ks_test_result_1 <- perform_ks_test_empirical(act_x1, x1)
  ks_test_result_2 <- perform_ks_test_empirical(act_x2, x2)
  
  # Perform Cramer-von Mises Test
  cvm_test_result_1 <- perform_cvm_test_empirical(act_x1, x1)
  cvm_test_result_2 <- perform_cvm_test_empirical(act_x2, x2)
  
  # Assess correlation changes
  corr_change_test <- evaluate_performance_corr(x1, x2, target_corr_kendall)
  
  # Return all results
  results <- list(
    tv_result_1 = tv_result_1,
    tv_result_2 = tv_result_2,
    ks_test_result_1 = ks_test_result_1,
    ks_test_result_2 = ks_test_result_2,
    cvm_test_result_1 = cvm_test_result_1,
    cvm_test_result_2 = cvm_test_result_2,
    corr_change_test = corr_change_test
  )
  
  return(results)
}


# PARAM SPACE
target_corr_kendall <- seq(-1, 1, 0.1)
datasets <- list()
inv_cdf_types <- c("quantile", "linear", "poly")
copula_types <- c("gaussian", "t")
bootstap_size = 50

results <- NULL
iter = 0
for (tau in target_corr_kendall) {
  for (inv_cdf_type in inv_cdf_types) {
    for (copula_type in copula_types) {
      for (dataset in names(datasets)) {
        for (boot_iter in 1:bootstap_size){
          df = datasets[[dataset]]
          result <- evaluate_copulas_and_inverse_cdfs_bootstrap(
            target_corr_kendall = tau,
            sampled_x1 = sample(df$x1, replace = T),
            sampled_x2 = sample(df$x2, replace = T), 
            act_x1 = df$x1,
            act_x2 = df$x2,
            copula_type = copula_type,
            inv_cdf_type = inv_cdf_type
          )
          print(iter)
          iter = iter + 1
          # Store results in a data frame
          temp_result <- data.frame(
            boot_iter = boot_iter,
            target_corr_kendall = tau,
            dataset = dataset,
            inv_cdf_type = inv_cdf_type,
            copula_type = copula_type,
            tv_val_1 = result$tv_result_1$tv_distance,
            tv_val_2 = result$tv_result_2$tv_distance,
            ks_stat_1 = result$ks_test_result_1$ks_test_result$statistic,
            ks_stat_2 = result$ks_test_result_2$ks_test_result$statistic,
            cvm_stat_1 = result$cvm_test_result_1$CramerVonMises$statistic,
            cvm_stat_2 = result$cvm_test_result_2$CramerVonMises$statistic,
            err_kendall = result$corr_change_test$err_kendall
          )
          results <- bind_rows(results, temp_result)
        }
      }
    }
  }
}





# Save the results
write.csv(results, "results/final_eval_bootstap.csv")



