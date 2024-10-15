# Load necessary library
library(MCMCpack)
library(GGally)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(akima)
library(rstan)

# Define the log of the joint PDF function
log_joint_pdf <- function(data) {
  # Extract x and y from data
  x <- data[1]
  y <- data[2]
  
  # Compute F1, F2, f1, f2
  F1_vals <- F1(x)
  F2_vals <- F2(y)
  f1_vals <- f1(x)
  f2_vals <- f2(y)
  
  # Define the joint PDF based on the provided distribution
  log_joint_pdf <- log(f1_vals) + log(f2_vals) + log((1 + a * (1 - 2 * F1_vals) * (1 - 2 * F2_vals)))
  
  # Return log of the joint PDF
  return(log_joint_pdf)
}

log_posterior <- function(params, data) {
  x <- params[1]
  y <- params[2]
  return(log_joint_pdf(c(x, y)))
}

mcmc_sample <- function(n, d1, d2){
  # Initial values for x and y
  x_init <- mean(d1)
  y_init <- mean(d2)
  
  # Combine x and y into a vector for initial values
  init_params <- c(x_init, y_init)
  
  # Run MCMC sampling using MCMCmetrop1R
  mcmc_result <- MCMCmetrop1R(
    log_posterior,
    theta.init = init_params,
    data = NULL,
    mcmc = n,
    burnin = 1000,
    thin = 1,
    tune = 0.5,
    verbose = TRUE
  )
  return(mcmc_result)
}

mcmc_sample_my <- function(n, d1, d2){
  # Initial values for x and y
  x_init <- mean(d1)
  y_init <- mean(d2)
  
  # Combine x and y into a vector for initial values
  init_params <- c(x_init, y_init)
  
  # Initialize the chain
  n_params <- 2
  samples <- matrix(0, nrow=n, ncol=n_params)
  samples[1, ] <- init_params
  
  # Metropolis-Hastings algorithm
  for (t in 2:n) {
    
    # Propose new state from the proposal distribution
    proposed_values <- rnorm(n_params, mean=samples[t-1, ], sd=10)
    
    # Calculate log-posterior for current and proposed states
    current_log_posterior <- log_posterior(samples[t-1, ])
    proposed_log_posterior <- log_posterior(proposed_values)
    
    # Calculate acceptance ratio
    acceptance_ratio <- exp(proposed_log_posterior - current_log_posterior)
    
    # Accept or reject the proposed state
    if (runif(1) <= acceptance_ratio) {
      samples[t, ] <- proposed_values
    } else {
      samples[t, ] <- samples[t-1, ]
    }
  }
  
  # Return the samples
  return(samples)
}


# Define the log-posterior function in a Stan model
stan_model_code <- "
functions {
  real log_joint_pdf(vector params) {
    real x = params[1];
    real y = params[2];
    return log(f1(x)) + log(f2(y)) + log(1 + a * (1 - 2 * F1(x)) * (1 - 2 * F2(y)));
  }
}
data {
  // No data required for this example
}
parameters {
  vector[2] params;
}
model {
  target += log_joint_pdf(params);  // Use log_joint_pdf to calculate log-posterior
}
"
stan_model <- stan_model(model_code = stan_model_code)



################################################################################
################################ INV FUNCTION ##################################
################################################################################

genCDFInv <- function(X) {
  ecdf1 <- ecdf(X)
  U <- sort(unique(ecdf1(X)))
  Finv <- sort(X)
  ak2 <- approxfun(U, Finv, method = "linear", rule = 2)
  return(ak2)
}

genCDF <- function(X) {
  ecdf1 <- ecdf(X)
  U <- sort(unique(ecdf1(X)))
  Finv <- sort(X)
  ak2 <- approxfun(Finv, U, method = "linear", rule = 2)
  return(ak2)
}

# Define the numerical derivative function
numerical_derivative <- function(F, x, h = 1e-5) {
  (F(x + h) - F(x - h)) / (2 * h)
}

# Function to create an empirical PDF from a dataset
get_empirical_pdf <- function(data, adjust = 1) {
  # Estimate the density using the density function
  density_est <- density(data, adjust = adjust)
  
  # Create an interpolation function based on the density estimate
  empirical_pdf <- approxfun(density_est$x, density_est$y, method = "linear", rule = 2)
  
  # Return the empirical PDF function
  return(empirical_pdf)
}

################################################################################
############################### SIMPLE EXAMPLE #################################
################################################################################

set.seed(42)
N <- 1000
d1 <- rnorm(N, mean = 2, sd = 3)
d2 <- rexp(N, rate = 0.4)

inv_cdf_d1 <- genCDFInv(d1)
inv_cdf_d2 <- genCDFInv(d2)
F1Inv <- Vectorize(inv_cdf_d1)
F2Inv <- Vectorize(inv_cdf_d2)

cdf_d1 <- genCDF(d1)
cdf_d2 <- genCDF(d2)
F1 <- Vectorize(cdf_d1)
F2 <- Vectorize(cdf_d2)

# Define f1 and f2 using the numerical derivative
# f1 <- function(x) numerical_derivative(F1, x)
# f2 <- function(y) numerical_derivative(F2, y)

# Define f1 and f2 using the empircal denisty
pdf_d1 <- get_empirical_pdf(d1)
pdf_d2 <- get_empirical_pdf(d2)
f1 <- Vectorize(pdf_d1)
f2 <- Vectorize(pdf_d2)

# Generate a sequence of x and y values to evaluate f1 and f2
x_vals <- seq(min(d1), max(d1), length.out = 100)
y_vals <- seq(min(d2), max(d2), length.out = 100)

# Compute f1 and f2 for these values
f1_vals <- sapply(x_vals, f1)
f2_vals <- sapply(y_vals, f2)

# Plot the empirical density of d1 and overlay f1
ggplot(data = data.frame(x = d1), aes(x = x)) +
  geom_density(aes(y = ..density..), color = "blue", fill = "blue", alpha = 0.3) +
  geom_line(data = data.frame(x = x_vals, f1 = f1_vals), aes(x = x, y = f1), color = "red", size = 1) +
  ggtitle("Empirical Density of d1 with f1(x) Overlay") +
  xlab("x") +
  ylab("Density / f1(x)")

# Plot the empirical density of d2 and overlay f2
ggplot(data = data.frame(y = d2), aes(x = y)) +
  geom_density(aes(y = ..density..), color = "blue", fill = "blue", alpha = 0.3) +
  geom_line(data = data.frame(y = y_vals, f2 = f2_vals), aes(x = y, y = f2), color = "red", size = 1) +
  ggtitle("Empirical Density of d2 with f2(y) Overlay") +
  xlab("y") +
  ylab("Density / f2(y)")


################################# COPULA GEN. ##################################

n <- 10000 
a <- 0.2
res = mcmc_sample(n, d1, d2)
x1 = res[,1]
x2 = res[,2]
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



