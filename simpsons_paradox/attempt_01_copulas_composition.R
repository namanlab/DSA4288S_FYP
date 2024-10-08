# Load necessary libraries
library(gridExtra)
library(ggExtra)
library(clue)
library(MCMCpack)
library(stringr)
library(GGally)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(akima)
library(actuar)
library(copula)
library(goftest)
options(scipen = 999)

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

genCDFInv_linear <- function(X) {
  ecdf1 <- ecdf(X)
  U <- sort((ecdf1(X)))
  Finv <- sort(X)
  ak2 <- approxfun(U, Finv, method = "linear", rule = 2)
  return(ak2)
}

# Step 1: Generate Data
set.seed(123)
n <- 1000

# Generate data for two categories of z
z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 1)
y <- -2*x + rnorm(n, 5, sd = 1)

F1Inv <- Vectorize(genCDFInv_linear(x))
F2Inv <- Vectorize(genCDFInv_linear(y))


t1 = -0.8
t2 = 0.8
t3 = -0.8
p1 = sin(t1*pi/2) 
p2 = sin(t2*pi/2) 
p3 = sin(t3*pi/2) 
p = c(p1, p2, p3)


################################################################################
################################################################################
################################################################################


## Approach 1: Baseline: fx = p1fx + p2fx + ... pnfx
## pi is pro of zi in z

res1 = gauss_copula_2(sum(z == "A"), p1)
x1 = F1Inv(res1[,1])
y1 = F2Inv(res1[,2])
z1 = rep("A", sum(z == "A"))
df_res1 = data.frame(x = x1, y = y1, z = z1)
res2 = gauss_copula_2(sum(z == "B"), p2)
x2 = F1Inv(res2[,1])
y2 = F2Inv(res2[,2])
z2 = rep("B", sum(z == "B"))
df_res2 = data.frame(x = x2, y = y2, z = z2)
res3 = gauss_copula_2(sum(z == "C"), p3)
x3 = F1Inv(res3[,1])
y3 = F2Inv(res3[,2])
z3 = rep("C", sum(z == "C"))
df_res3 = data.frame(x = x3, y = y3, z = z3)
df_res <- rbind(df_res1, df_res2, df_res3)

pl1 = data.frame(x = x, y = y, z= z) %>% ggplot() +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Original Distribution")
pl2 = df_res %>% ggplot() +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Modified Distribution")
pl1 = ggMarginal(pl1, type="density")
pl2 = ggMarginal(pl2, type="density")
grid.arrange(pl1, pl2, ncol = 2)

for (i in unique(df_res$z)){
  cur_df = df_res %>% filter(z == i)
  print(cor(cur_df$x, cur_df$y, method = "kendall"))
}
print(cor(df_res$x, df_res$y, method = "kendall"))
  

################################################################################
################################################################################
################################################################################


## Approach 2: fx = f1 + f1 + ... fn
## fi is between region iwth area = pi

# Split the distribution of x into quantile slices based on categories in z
z_levels <- sort(unique(z))
z <- sort(z)
x_quantiles <- quantile(x, probs = cumsum(table(z)/length(z)))

df_res <- data.frame()

for (i in seq_along(z_levels)) {
  z_cat <- z_levels[i]
  prob_val = as.vector(table(z)[i]/length(z))
  
  # Generate samples using the Gaussian copula
  p_corr <- p[i]
  n_samples <- sum(z == z_cat)
  copula_samples <- gauss_copula_2(n_samples, p_corr)
  
  # Filter x_samples based on x_range
  filtered_idx <- if (i == 1) {
    which(x <= x_quantiles[i])
  } else {
    which(x > x_quantiles[i - 1] & x <= x_quantiles[i])
  }
  x_samples <- x[filtered_idx]
  
  F1Inv <- Vectorize(genCDFInv_linear(x_samples))
  F2Inv <- Vectorize(genCDFInv_linear(y))
  
  # Transform samples using the inverse CDFs
  x_samples <- F1Inv(copula_samples[, 1])
  y_samples <- F2Inv(copula_samples[, 2])
  
  # Append results
  df_res <- rbind(df_res, data.frame(x = x_samples, y = y_samples, z = z_cat))
}


pl1 = data.frame(x = x, y = y, z= z) %>% ggplot() +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Original Distribution")
pl2 = df_res %>% ggplot() +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Modified Distribution")
pl1 = ggMarginal(pl1, type="density")
pl2 = ggMarginal(pl2, type="density")
grid.arrange(pl1, pl2, ncol = 2)

for (i in unique(df_res$z)){
  cur_df = df_res %>% filter(z == i)
  print(cor(cur_df$x, cur_df$y, method = "kendall"))
}
print(cor(df_res$x, df_res$y, method = "kendall"))


################################################################################
################################################################################
################################################################################


## Approach 3: fx = f1 + f1 + ... fn, fy = f1 + ... + fn

# Split the distribution of x into quantile slices based on categories in z
z_levels <- sort(unique(z))
z <- sort(z)
x_quantiles <- quantile(x, probs = cumsum(table(z)/length(z)))
y_quantiles <- quantile(y, probs = cumsum(table(z)/length(z)))

df_res <- data.frame()

for (i in seq_along(z_levels)) {
  z_cat <- z_levels[i]
  prob_val = as.vector(table(z)[i]/length(z))
  
  # Generate samples using the Gaussian copula
  p_corr <- p[i]
  n_samples <- sum(z == z_cat)
  copula_samples <- gauss_copula_2(n_samples, p_corr)
  
  # Filter x_samples based on x_range
  filtered_idx_x <- if (i == 1) {
    which(x <= x_quantiles[i])
  } else {
    which(x > x_quantiles[i - 1] & x <= x_quantiles[i])
  }
  filtered_idx_y <- if (i == 1) {
    which(y <= y_quantiles[i])
  } else {
    which(y > y_quantiles[i - 1] & y <= y_quantiles[i])
  }
  x_samples <- x[filtered_idx_x]
  y_samples <- y[filtered_idx_y]
  
  F1Inv <- Vectorize(genCDFInv_linear(x_samples))
  F2Inv <- Vectorize(genCDFInv_linear(y_samples))
  
  # Transform samples using the inverse CDFs
  x_samples <- F1Inv(copula_samples[, 1])
  y_samples <- F2Inv(copula_samples[, 2])
  
  # Append results
  df_res <- rbind(df_res, data.frame(x = x_samples, y = y_samples, z = z_cat))
}

pl1 = data.frame(x = x, y = y, z= z) %>% ggplot() +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Original Distribution")
pl2 = df_res %>% ggplot() +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Modified Distribution")
pl1 = ggMarginal(pl1, type="density")
pl2 = ggMarginal(pl2, type="density")
grid.arrange(pl1, pl2, ncol = 2)

for (i in unique(df_res$z)){
  cur_df = df_res %>% filter(z == i)
  print(cor(cur_df$x, cur_df$y, method = "kendall"))
}
print(cor(df_res$x, df_res$y, method = "kendall"))




################################################################################
################################################################################
################################################################################


## Approach 4: 3.5 + allow sum variable boundaries

set.seed(123)
n <- 1000

# Generate data for two categories of z
z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 1)
y <- -2*x + rnorm(n, 5, sd = 1)

t = c(-0.8, -0.8, -0.8) # c(0.8, -0.8, 0.6)
p = sin(t*pi/2) 

# Split the distribution of x into quantile slices based on categories in z
z_levels <- sort(unique(z))
z <- sort(z)
x_quantiles <- quantile(x, probs = cumsum(table(z)/length(z)))
y_quantiles <- quantile(y, probs = cumsum(table(z)/length(z)))
eps_x <- 1
eps_y <- 1

df_res <- data.frame()

for (i in seq_along(z_levels)) {
  z_cat <- z_levels[i]
  prob_val = as.vector(table(z)[i]/length(z))
  
  # Generate samples using the Gaussian copula
  p_corr <- p[i]
  n_samples <- sum(z == z_cat)
  copula_samples <- gauss_copula_2(n_samples, p_corr)
  
  # Filter x_samples based on x_range
  filtered_idx_x <- if (i == 1) {
    which(x <= x_quantiles[i] + eps_x)
  } else {
    which(x > x_quantiles[i - 1] - eps_x & x <= x_quantiles[i] + eps_x)
  }
  filtered_idx_y <- if (i == 1) {
    which(y <= y_quantiles[i] + eps_y)
  } else {
    which(y > y_quantiles[i - 1] - eps_y  & y <= y_quantiles[i] + eps_y)
  }
  x_samples <- x[filtered_idx_x]
  y_samples <- y[filtered_idx_y]
  
  F1Inv <- Vectorize(genCDFInv_linear(x_samples))
  F2Inv <- Vectorize(genCDFInv_linear(y_samples))
  
  # Transform samples using the inverse CDFs
  x_samples <- F1Inv(copula_samples[, 1])
  y_samples <- F2Inv(copula_samples[, 2])
  
  # Append results
  df_res <- rbind(df_res, data.frame(x = x_samples, y = y_samples, z = z_cat))
}

p = df_res %>% ggplot() +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() +
  theme(legend.position = "bottom")
ggMarginal(p, type="density")
for (i in unique(df_res$z)){
  cur_df = df_res %>% filter(z == i)
  print(cor(cur_df$x, cur_df$y, method = "kendall"))
}
print(cor(df_res$x, df_res$y, method = "kendall"))



## Approach 5: 4 but preserve densities

set.seed(123)
n <- 1000

# Generate data for two categories of z
z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 1)
y <- -2*x + rnorm(n, 5, sd = 1)

t = c(-0.8, -0.8, -0.8) # c(0.8, -0.8, 0.6)
p = sin(t*pi/2) 

# Split the distribution of x into quantile slices based on categories in z
z_levels <- sort(unique(z))
z <- sort(z)
x_quantiles <- quantile(x, probs = cumsum(table(z)/length(z)))
y_quantiles <- quantile(y, probs = cumsum(table(z)/length(z)))
eps_x <- 2
eps_y <- 2

df_res <- data.frame()

for (i in seq_along(z_levels)) {
  z_cat <- z_levels[i]
  prob_val = as.vector(table(z)[i]/length(z))
  
  # Generate samples using the Gaussian copula
  p_corr <- p[i]
  n_samples <- sum(z == z_cat)
  copula_samples <- gauss_copula_2(n_samples, p_corr)
  
  # Filter x_samples based on x_range
  filtered_idx_x <- if (i == 1) {
    which(x <= x_quantiles[i] + eps_x)
  } else {
    which(x > x_quantiles[i - 1] - eps_x & x <= x_quantiles[i] + eps_x)
  }
  filtered_idx_y <- if (i == 1) {
    which(y <= y_quantiles[i] + eps_y)
  } else {
    which(y > y_quantiles[i - 1] - eps_y  & y <= y_quantiles[i] + eps_y)
  }
  x_samples <- x[filtered_idx_x]
  y_samples <- y[filtered_idx_y]
  
  F1Inv <- Vectorize(genCDFInv_linear(x_samples))
  F2Inv <- Vectorize(genCDFInv_linear(y_samples))
  
  # Transform samples using the inverse CDFs
  x_samples <- F1Inv(copula_samples[, 1])
  y_samples <- F2Inv(copula_samples[, 2])
  
  # Append results
  df_res <- rbind(df_res, data.frame(x = x_samples, y = y_samples, z = z_cat))
}

p = df_res %>% ggplot() +
  geom_point(aes(x = x, y = y, color = z)) +
  theme_bw() +
  theme(legend.position = "bottom")
ggMarginal(p, type="density")
for (i in unique(df_res$z)){
  cur_df = df_res %>% filter(z == i)
  print(cor(cur_df$x, cur_df$y, method = "kendall"))
}
print(cor(df_res$x, df_res$y, method = "kendall"))

