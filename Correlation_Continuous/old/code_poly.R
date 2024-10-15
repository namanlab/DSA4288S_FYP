# Load necessary library
library(MCMCpack)
library(GGally)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(akima)

# Function to generate i.i.d. samples from a Gaussian copula
generate_gaussian_copula_samples <- function(n, d, rho_matrix) {
  # Step 1: Generate multivariate normal samples
  mean_vector <- rep(0, d)  # Mean vector for multivariate normal
  mvn_samples <- mvrnorm(n = n, mu = mean_vector, Sigma = rho_matrix)
  print(dim(mvn_samples))
  # Step 2: Transform samples to uniform using the CDF of the standard normal
  uniform_samples <- pnorm(mvn_samples)
  return(uniform_samples)
}
gauss_copula_2 <- function(n, p){
  rho_matrix <- matrix(c(
    1, p,
    p, 1
  ), nrow = 2, ncol = 2)
  smpls = generate_gaussian_copula_samples(n, 2, rho_matrix)
  print(ggpairs(as.data.frame(smpls)))
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
t_copula_2 <- function(n, p){
  rho_matrix <- matrix(c(
    1, p,
    p, 1
  ), nrow = 2, ncol = 2)
  smpls = generate_t_copula_samples(n, 2, rho_matrix, 5)
  print(ggpairs(as.data.frame(smpls)))
  return(smpls)
}

################################################################################
################################ INV FUNCTION ##################################
################################################################################

# Function to generate inverse CDF using polynomial fitting and polyroot
genCDFInv <- function(data, degree = 10) {
  x <- sort(unique(data))
  y <- ecdf(data)(x)
  poly_fit <- lm(x ~ poly(y, degree = degree, raw = TRUE))
  cdf_poly <- function(y) predict(poly_fit, newdata = data.frame(y = y))
  return(cdf_poly)
}


################################################################################
############################### SIMPLE EXAMPLE #################################
################################################################################

set.seed(42)
N <- 300
d1 <- rnorm(N, mean = 2, sd = 3)
d2 <- rexp(N, rate = 0.4)

# Plotting the distributions
par(mfrow = c(1, 2))  # Set up a 1x2 plotting area
hist(d1, breaks = 100, probability = TRUE, 
     main = "Complex Mixture Distribution 1", 
     xlab = "Value", col = "skyblue", border = "white")
lines(density(d1), col = "darkblue", lwd = 2)
hist(d2, breaks = 100, probability = TRUE, 
     main = "Complex Mixture Distribution 2", 
     xlab = "Value", col = "lightgreen", border = "white")
lines(density(d2), col = "darkgreen", lwd = 2)

inv_cdf_d1 <- genCDFInv(d1)
inv_cdf_d2 <- genCDFInv(d2)
F1Inv <- Vectorize(inv_cdf_d1)
F2Inv <- Vectorize(inv_cdf_d2)


################################# PLOT CHECKS ##################################

set.seed(42)
N_test <- 10000
d1_test <- rnorm(N_test, mean = 2, sd = 3)
d2_test <- rexp(N_test, rate = 0.4)

x_vals <- seq(0, 1, length.out = N_test)
ecdf_vals_d1 <- quantile(d1_test, x_vals)
cdf_poly_vals_d1 <- F1Inv(x_vals)
df1 = data.frame(x = rep(x_vals, 2),
                 y = c(ecdf_vals_d1, cdf_poly_vals_d1),
                 type = rep(c("Empirical ICDF", "Polynomial ICDF"), each = length(x_vals)), var = "1") 
x_vals <- seq(0, 1, length.out = N_test)
ecdf_vals_d2 <- quantile(d2_test, x_vals)
cdf_poly_vals_d2 <- F2Inv(x_vals)
df2 = data.frame(x = rep(x_vals, 2),
                 y = c(ecdf_vals_d2, cdf_poly_vals_d2),
                 type = rep(c("Empirical ICDF", "Polynomial ICDF"), each = length(x_vals)), var = "2") 
df <- rbind(df1, df2)
rmse_d1 <- sqrt(mean((ecdf_vals_d1 - cdf_poly_vals_d1)^2))
max_error_d1 <- max(abs(ecdf_vals_d1 - cdf_poly_vals_d1))
rmse_d2 <- sqrt(mean((ecdf_vals_d2 - cdf_poly_vals_d2)^2))
max_error_d2 <- max(abs(ecdf_vals_d2 - cdf_poly_vals_d2))
facet_labels <- c(paste0("Var 1\nRMSE: ", round(rmse_d1, 4), "\nMax Error: ", round(max_error_d1, 4)),
                  paste0("Var 2\nRMSE: ", round(rmse_d2, 4), "\nMax Error: ", round(max_error_d2, 4)))
names(facet_labels) <- c("1", "2")
ggplot(df, aes(x = x, y = y, color = type)) +
  geom_line(size = 1) +
  labs(title = "Empirical CDF vs Polynomial-fitted CDF",
       x = "x", y = "CDF", color = "Type") + 
  theme_minimal() + 
  facet_wrap(~var, scales = "free", labeller = as_labeller(facet_labels)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 60))

################################# COPULA GEN. ##################################

n <- 1000 
res = gauss_copula_2(n, -0.8)
x1 = F1Inv(res[,1])
x2 = F2Inv(res[,2])
ggpairs(data.frame(x1 = x1, x2 = x2))

n <- 1000 
res = t_copula_2(n, 0.8)
x1 = F1Inv(res[,1])
x2 = F2Inv(res[,2])
ggpairs(data.frame(x1 = x1, x2 = x2))


################################################################################
################################ MORE COMPLEX ##################################
################################################################################

set.seed(42)
N <- 300
# First complex mixture distribution
d1 <- c(rnorm(N * 0.4, mean = 2, sd = 3), 
        rexp(N * 0.3, rate = 0.4), 
        runif(N * 0.3, min = -10, max = 10))

# Second complex mixture distribution
d2 <- c(rnorm(N * 0.3, mean = -5, sd = 2), 
        rnorm(N * 0.3, mean = 5, sd = 2), 
        rchisq(N * 0.4, df = 4))

# Plotting the distributions
par(mfrow = c(1, 2))  # Set up a 1x2 plotting area
hist(d1, breaks = 100, probability = TRUE, 
     main = "Complex Mixture Distribution 1", 
     xlab = "Value", col = "skyblue", border = "white")
lines(density(d1), col = "darkblue", lwd = 2)
hist(d2, breaks = 100, probability = TRUE, 
     main = "Complex Mixture Distribution 2", 
     xlab = "Value", col = "lightgreen", border = "white")
lines(density(d2), col = "darkgreen", lwd = 2)


inv_cdf_d1 <- genCDFInv(d1)
inv_cdf_d2 <- genCDFInv(d2)
F1Inv <- Vectorize(inv_cdf_d1)
F2Inv <- Vectorize(inv_cdf_d2)

################################# PLOT CHECKS ##################################

set.seed(42)
N_test <- 10000
d1_test <- c(rnorm(N_test * 0.4, mean = 2, sd = 3), 
             rexp(N_test * 0.3, rate = 0.4), 
             runif(N_test * 0.3, min = -10, max = 10))
d2_test <- c(rnorm(N_test * 0.3, mean = -5, sd = 2), 
             rnorm(N_test * 0.3, mean = 5, sd = 2), 
             rchisq(N_test * 0.4, df = 4))

x_vals <- seq(0, 1, length.out = N_test)
ecdf_vals_d1 <- quantile(d1_test, x_vals)
cdf_poly_vals_d1 <- F1Inv(x_vals)
df1 = data.frame(x = rep(x_vals, 2),
                 y = c(ecdf_vals_d1, cdf_poly_vals_d1),
                 type = rep(c("Empirical ICDF", "Polynomial ICDF"), each = length(x_vals)), var = "1") 
x_vals <- seq(0, 1, length.out = N_test)
ecdf_vals_d2 <- quantile(d2_test, x_vals)
cdf_poly_vals_d2 <- F2Inv(x_vals)
df2 = data.frame(x = rep(x_vals, 2),
                 y = c(ecdf_vals_d2, cdf_poly_vals_d2),
                 type = rep(c("Empirical ICDF", "Polynomial ICDF"), each = length(x_vals)), var = "2") 
df <- rbind(df1, df2)
rmse_d1 <- sqrt(mean((ecdf_vals_d1 - cdf_poly_vals_d1)^2))
max_error_d1 <- max(abs(ecdf_vals_d1 - cdf_poly_vals_d1))
rmse_d2 <- sqrt(mean((ecdf_vals_d2 - cdf_poly_vals_d2)^2))
max_error_d2 <- max(abs(ecdf_vals_d2 - cdf_poly_vals_d2))
facet_labels <- c(paste0("Var 1\nRMSE: ", round(rmse_d1, 4), "\nMax Error: ", round(max_error_d1, 4)),
                  paste0("Var 2\nRMSE: ", round(rmse_d2, 4), "\nMax Error: ", round(max_error_d2, 4)))
names(facet_labels) <- c("1", "2")
ggplot(df, aes(x = x, y = y, color = type)) +
  geom_line(size = 1) +
  labs(title = "Empirical CDF vs Polynomial-fitted CDF",
       x = "x", y = "CDF", color = "Type") + 
  theme_minimal() + 
  facet_wrap(~var, scales = "free", labeller = as_labeller(facet_labels)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 60))



################################# COPULA GEN. ##################################


n <- 10000 
res = gauss_copula_2(n, -0.8)
x1 = F1Inv(res[,1])
x2 = F2Inv(res[,2])
ggpairs(data.frame(x1 = x1, x2 = x2))


n <- 10000 
res = t_copula_2(n, 0.8)
x1 = F1Inv(res[,1])
x2 = F2Inv(res[,2])
ggpairs(data.frame(x1 = x1, x2 = x2))



# TRY: 
# - USE SPEARMAN RANK INSTEDA OF PEARSON FOR. ALL, AND CONVERT BY FORMULAS
# - Use gaussian quadrature  or chebyshev polynomials for interpolation
# what is akima interpolation
# find any other suitable fucntions or inverse cdfs ...
# - STEFFEN SPLINES: https://math.stackexchange.com/questions/4162828/interpolation-with-exact-inverse


