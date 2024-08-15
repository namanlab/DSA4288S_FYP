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


################################################################################
################################ COPULAS FUNC ##################################
################################################################################


# Function to generate i.i.d. samples from a Gaussian copula
generate_gaussian_copula_samples <- function(n, d, rho_matrix) {
  # Step 1: Generate multivariate normal samples
  mean_vector <- rep(0, d)  # Mean vector for multivariate normal
  mvn_samples <- mvrnorm(n = n, mu = mean_vector, Sigma = rho_matrix)
  
  # Step 2: Transform samples to uniform using the CDF of the standard normal
  uniform_samples <- pnorm(mvn_samples)
  
  return(uniform_samples)
}
# Generate samples from Gaussian copula and plot pairwise scatterplot matrix
gauss_copula_2 <- function(n, p) {
  rho_matrix <- matrix(c(1, p, p, 1), nrow = 2, ncol = 2)
  smpls <- generate_gaussian_copula_samples(n, 2, rho_matrix)
  return(smpls)
}

# Function to generate i.i.d. samples from a t-copula
generate_t_copula_samples <- function(n, d, rho_matrix, df) {
  # Step 1: Generate multivariate t samples
  mvn_samples <- rmvt(n = n, sigma = rho_matrix, df = df)
  
  # Step 2: Transform samples to uniform using the CDF of the t-distribution
  uniform_samples <- pt(mvn_samples, df = df)
  
  return(uniform_samples)
}
# Generate samples from t-copula and plot pairwise scatterplot matrix
t_copula_2 <- function(n, p) {
  rho_matrix <- matrix(c(1, p, p, 1), nrow = 2, ncol = 2)
  smpls <- generate_t_copula_samples(n, 2, rho_matrix, 5)
  return(smpls)
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
evaluate_inv_cdf <- function(est_inv_cdf, N_test, true_gen_fn) {
  x_test <- true_gen_fn(N_test)
  x_vals <- seq(0, 1, length.out = N_test)
  samples <- quantile(x_test, x_vals)
  estimated <- est_inv_cdf(x_vals)
  mse <- mean((samples - estimated)^2)
  mae <- mean(abs(samples - estimated))
  list(MSE = mse, MAE = mae)
}

# 2. Function to Perform Kolmogorov-Smirnov Test
perform_ks_test <- function(original_data, generated_data) {
  ks_test_result <- ks.test(original_data, generated_data)
  list(ks_test_result = ks_test_result)
}

# 3. Function to Perform Cramer-von Mises Test
perform_goodness_of_fit_tests <- function(generated_data, N_test, true_gen_fn) {
  x_test <- true_gen_fn(N_test)
  Fn = ecdf(x_test)
  cvm_test <- cvm.test(generated_data, Fn)
  list(CramerVonMises = cvm_test)
}

# 4. Function to Evaluate Correlation Performance
evaluate_performance_corr <- function(x1, x2, target_corr_kendall) {
  kendall_tau <- cor(x1, x2, method = "kendall")
  err_kendall <- kendall_tau - target_corr_kendall
  return(list(err_kendall = err_kendall))
}



################################################################################
################################ FINAL EVAL FN #################################
################################################################################



evaluate_copulas_and_inverse_cdfs <- function(target_corr_kendall,
                                              distribution_1, distribution_2, copula_type, 
                                              inv_cdf_type, n = 1000, n_test = 10000) {
  
  # Step 1: Generate original data
  d1 <- distribution_1(n)
  d2 <- distribution_2(n)
  
  # Step 2: Generate copula samples
  if (copula_type == "gaussian") {
    target_rho = sin(target_corr_kendall*pi/2) # Greiner's equality 
    res <- gauss_copula_2(n, target_rho) 
  } else if (copula_type == "t") {
    target_rho = sin(target_corr_kendall*pi/2)# Replace with some form ???
    res <- t_copula_2(n, target_rho)
  } else {
    stop("Unsupported copula type")
  }
  
  # Step 3: Apply chosen inverse CDF
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
  
  # Step 4: Evaluate performance of inverse CDFs
  inv_cdf_eval_1 <- evaluate_inv_cdf(inv_cdf_d1, n_test, distribution_1)
  inv_cdf_eval_2 <- evaluate_inv_cdf(inv_cdf_d2, n_test, distribution_2)
  
  # Step 5: Perform Kolmogorov-Smirnov Test
  ks_test_result_1 <- perform_ks_test(d1, x1)
  ks_test_result_2 <- perform_ks_test(d2, x2)
  
  # Step 6: Perform Cramer-von Mises Test
  goodness_of_fit_tests_1 <- perform_goodness_of_fit_tests(x1, n_test, distribution_1)
  goodness_of_fit_tests_2 <- perform_goodness_of_fit_tests(x2, n_test, distribution_2)
  
  # Step 7: Assess correlation changes
  corr_change_test <- evaluate_performance_corr(x1, x2, target_corr_kendall)
  
  # Step 8: Return all results
  results <- list(
    inv_cdf_eval_1 = inv_cdf_eval_1,
    inv_cdf_eval_2 = inv_cdf_eval_2,
    ks_test_result_1 = ks_test_result_1,
    ks_test_result_2 = ks_test_result_2,
    goodness_of_fit_tests_1 = goodness_of_fit_tests_1,
    goodness_of_fit_tests_2 = goodness_of_fit_tests_2,
    corr_change_test = corr_change_test
  )
  
  return(results)
}


################################################################################
############################## FINAL EVAL CASES ################################
################################################################################

set.seed(123)

distribution_normal <- function(n) rnorm(n, mean = 0, sd = 1)
distribution_exponential <- function(n) rexp(n, rate = 1)
distribution_uniform <- function(n) runif(n, min = 0, max = 1)
distribution_log_normal <- function(n) rlnorm(n, meanlog = 0, sdlog = 1)
distribution_beta <- function(n) rbeta(n, shape1 = 2, shape2 = 5)
distribution_mixture1 <- function(n) {
  prob <- runif(n)
  rnorm(n, mean = ifelse(prob > 0.5, 3, -2), sd = 1)
}
distribution_mixture2 <- function(n) {
  prob <- runif(n)
  ifelse(prob > 0.5, rnorm(n, mean = 0, sd = 1), rlnorm(n, meanlog = 0, sdlog = 0.5))
}

target_corr_kendall <- seq(-1, 1, 0.1)
distributions <- list(
  normal = distribution_normal,
  exponential = distribution_exponential,
  uniform = distribution_uniform,
  log_normal = distribution_log_normal,
  beta = distribution_beta,
  mixture1 = distribution_mixture1,
  mixture2 = distribution_mixture2
)
inv_cdf_types <- c("quantile", "linear", "poly")
copula_types <- c("gaussian", "t")


library(dplyr)

results <- NULL

iter = 0
for (tau in target_corr_kendall) {
  for (dist_name_1 in names(distributions)) {
    for (dist_name_2 in names(distributions)) {
      for (inv_cdf_type in inv_cdf_types) {
        for (copula_type in copula_types) {
          
          result <- evaluate_copulas_and_inverse_cdfs(
            target_corr_kendall = tau,
            distribution_1 = distributions[[dist_name_1]],
            distribution_2 = distributions[[dist_name_2]],
            copula_type = copula_type,
            inv_cdf_type = inv_cdf_type
          )
          
          print(iter)
          iter = iter + 1
          
          # Store results in a data frame
          temp_result <- data.frame(
            target_corr_kendall = tau,
            distribution_1 = dist_name_1,
            distribution_2 = dist_name_2,
            inv_cdf_type = inv_cdf_type,
            copula_type = copula_type,
            mse_inv_cdf_1 = result$inv_cdf_eval_1$MSE,
            mae_inv_cdf_1 = result$inv_cdf_eval_1$MAE,
            mse_inv_cdf_2 = result$inv_cdf_eval_2$MSE,
            mae_inv_cdf_2 = result$inv_cdf_eval_2$MAE,
            ks_stat_1 = result$ks_test_result_1$ks_test_result$statistic,
            ks_pval_1 = result$ks_test_result_1$ks_test_result$p.value,
            ks_stat_2 = result$ks_test_result_2$ks_test_result$statistic,
            ks_pval_2 = result$ks_test_result_2$ks_test_result$p.value,
            cvm_stat_1 = result$goodness_of_fit_tests_1$CramerVonMises$statistic,
            cvm_stat_2 = result$goodness_of_fit_tests_2$CramerVonMises$statistic,
            err_kendall = result$corr_change_test$err_kendall
          )
          
          results <- bind_rows(results, temp_result)
        }
      }
    }
  }
}


# View the results
print(results)
  
results %>% dplyr::select(inv_cdf_type, distribution_1, copula_type, mse_inv_cdf_1, mae_inv_cdf_1) %>%
  pivot_longer(cols = c(mse_inv_cdf_1, mae_inv_cdf_1), names_to = "err_type", values_to = "err_val") %>%
  ggplot() +
  geom_boxplot(aes(x = inv_cdf_type, y = err_val, fill = copula_type)) +
  scale_y_log10() +
  facet_wrap(~err_type) +
  theme_bw()

options(scipen = 999)
# MSE
results %>% group_by(copula_type, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_ks_pval = mean(c(ks_pval_1, ks_pval_2)),
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
            mean_ks_pval = mean(c(ks_pval_1, ks_pval_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_mae_inv_cdf)) +
  facet_grid(copula_type~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# KS P
results %>% group_by(copula_type, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_ks_pval = mean(c(ks_pval_1, ks_pval_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_ks_pval)) +
  facet_grid(copula_type~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

# CVM
results %>% group_by(copula_type, inv_cdf_type, distribution_1, distribution_2) %>%
  summarise(mean_mse_inv_cdf = mean(c(mse_inv_cdf_1, mse_inv_cdf_2)),
            mean_mae_inv_cdf = mean(c(mae_inv_cdf_1, mae_inv_cdf_2)),
            mean_ks_pval = mean(c(ks_pval_1, ks_pval_2)),
            mean_cvm_stat = mean(c(cvm_stat_1, cvm_stat_2))) %>%
  ggplot() +
  geom_tile(aes(x = distribution_1, y = distribution_2, fill = mean_cvm_stat)) +
  facet_grid(copula_type~inv_cdf_type) +
  guides(fill = guide_legend(override.aes = list(size = 10)))  +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(trans = "log10", low = "white", high = "darkred", mid = "red") +
  labs(fill = "Err", x = "", y = "")

