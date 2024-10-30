devtools::load_all(".")

# Libraries needed
if (!require(mvtnorm)) install.packages("mvtnorm")
library(mvtnorm)  # For multivariate normal distribution
if (!require(ggplot2)) install.packages("ggplot2")
library(ggplot2) # For box plot distribution
if (!require(dplyr)) install.packages("dplyr")
library(dplyr)
if (!require(tidyr)) install.packages("tidyr")
library(tidyr)
if (!require(pryr)) install.packages("pryr")
library(pryr) # For memory usage

set.seed(123)

## Data generation process
# Simulation parameters
n <- 600              # Total population size per simulation
K <- 50               # Number of simulations
beta <- c(4, -2, -1)
covmat <- matrix(0.5, ncol = 2, nrow = 2)
diag(covmat) <- c(0.7, 1.3)
g1 <- 20              # Number of domains for d1
g2 <- 30              # Number of domains for d2

# Preallocate storage for simulation results
population_data <- vector("list", K)
sampled_data_uninformative <- vector("list", K)
sampled_data_informative <- vector("list", K)

# Perform simulations
for (k in 1:K) {
  X1           <- rnorm(n, mean = 2)
  X2           <- rexp(n, rate = 1)

  group1 <- rep(1:g1, length.out = n)
  group2 <- rep(1:g2, length.out = n)
  re1    <- rnorm(g1, sd = 2)
  re2    <- rmvnorm(g2, sigma = covmat)

  modX  <- model.matrix( ~ X1 + X2)
  modZ1 <- model.matrix( ~ -1 + as.factor(group1))
  modZ2A <- model.matrix( ~ -1 + as.factor(group2))
  modZ2B <- model.matrix( ~ -1 + X1:as.factor(group2))

  eta    <- modX %*% beta + modZ1 %*% re1 +
    modZ2A %*% as.vector(re2[,1]) +
    modZ2B %*% as.vector(re2[,2])

  lin <- eta + rnorm(n, sd = 2.3)
  prob <- 1 / (1 + exp(-eta))
  bin <- rbinom(n, size = 1, prob = prob)

  dfsamp <- data.frame(group1 = group1, group2 = group2, X1 = X1, X2 = X2, lin = lin, bin = bin, eta = eta)
  population_data[[k]] <- dfsamp
}

# Predefined lists to store results
results_list_uninf <- vector("list", K)
time_stats_uninf <- vector("list", K)
memory_stats_uninf <- vector("list", K)
cpu_usage_stats_uninf <- vector("list", K)
results_list_inf <- vector("list", K)
time_stats_inf <- vector("list", K)
memory_stats_inf <- vector("list", K)
cpu_usage_stats_inf <- vector("list", K)

# Step 2: Sampling
for (k in 1:K) {
  print(k)
  data_k <- population_data[[k]]

  # Non-informative sampling
  sampled_data_uninformative[[k]] <- do.call("rbind", lapply(split(data_k, data_k$group1), function(df) {
    probs <- df$X2 / sum(df$X2)

    sampled_indices <- sample(nrow(df), 5, replace = FALSE, prob = probs)
    sampled_df <- df[sampled_indices, ]

    sampled_df$inclusion_prob <- probs[sampled_indices]
    sampled_df$weights <- 1 / sampled_df$inclusion_prob

    return(sampled_df)
  }))

  # Using system.time to measure execution time
  result <- system.time({
    # Start CPU time tracking
    start_time <- proc.time()
    # Check memory usage before running the code
    initial_mem <- mem_used()

    # Function call here
    out <- wglmm( bin ~ X1 + X2 + (1|group1) + (1+X1|group2),
                  data = sampled_data_uninformative[[k]], trace = TRUE, iter1 = 250,
                  iter2 = 1001, MI = 2000, family = binomial(), e_step_method = "e_step_gibbs_sampling")

    # End CPU time tracking
    end_time <- proc.time()
    # Calculate CPU time used
    cpu_usage <- end_time - start_time
    # Check memory usage after running the code
    final_mem <- mem_used()
    # Calculate the difference in memory usage
    mem_change <- final_mem - initial_mem
    mem_change_bytes <- as.numeric(mem_change)
    mem_change_kb <- mem_change_bytes / 1024
  })

  results_list_uninf[[k]] <- out
  time_stats_uninf[[k]] <- result
  memory_stats_uninf[[k]] <- mem_change_kb
  cpu_usage_stats_uninf[[k]] <- cpu_usage

  # Informative sampling
  sampled_data_informative[[k]] <- do.call("rbind", lapply(split(data_k, data_k$group1), function(df) {
    residuals <- abs(df$lin - df$eta)

    quantiles <- quantile(residuals, probs = c(0.25, 0.5, 0.75))

    probs <- ifelse(residuals <= quantiles[1], 0.1,
                    ifelse(residuals <= quantiles[2], 0.2, 0.4))

    sampled_indices <- sample(nrow(df), 5, replace = FALSE, prob = probs)
    sampled_df <- df[sampled_indices, ]

    sampled_df$inclusion_prob <- probs[sampled_indices] * (5/30)
    sampled_df$weights <- 1 / sampled_df$inclusion_prob

    return(sampled_df)
  }))

  # Using system.time to measure execution time
  result <- system.time({
    # Start CPU time tracking
    start_time <- proc.time()
    # Check memory usage before running the code
    initial_mem <- mem_used()

    out <- wglmm( bin ~ X1 + X2 + (1|group1) + (1+X1|group2),
                  data = sampled_data_informative[[k]], trace = TRUE, iter1 = 250,
                  iter2 = 1001, MI = 2000, family = binomial(), e_step_method = "e_step_gibbs_sampling")

    # End CPU time tracking
    end_time <- proc.time()
    # Calculate CPU time used
    cpu_usage <- end_time - start_time
    # Check memory usage after running the code
    final_mem <- mem_used()
    # Calculate the difference in memory usage
    mem_change <- final_mem - initial_mem
    mem_change_bytes <- as.numeric(mem_change)
    mem_change_kb <- mem_change_bytes / 1024
  })
  results_list_inf[[k]] <- out
  time_stats_inf[[k]] <- result
  memory_stats_inf[[k]] <- mem_change_kb
  cpu_usage_stats_inf[[k]] <- cpu_usage
}

# # IMP LIN
# out <- wglmm( lin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#               data = subset(sampled_data_uninformative[[k]], select = -c(inclusion_prob, weights)), trace = TRUE, iter1 = 250,
#               iter2 = 1001, MI = 2000, family = gaussian(), e_step_method = "e_step_importance_sampling", weights = sampled_data_uninformative[[k]]$weights)
# out <- wglmm( lin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#               data = subset(sampled_data_informative[[k]], select = -c(inclusion_prob, weights)), trace = TRUE, iter1 = 250,
#               iter2 = 1001, MI = 2000, family = gaussian(), e_step_method = "e_step_importance_sampling", weights = sampled_data_informative[[k]]$weights)
# # IMP BIN
# out <- wglmm( bin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#               data = subset(sampled_data_uninformative[[k]], select = -c(inclusion_prob, weights)), trace = TRUE, iter1 = 250,
#               iter2 = 1001, MI = 2000, family = binomial(), e_step_method = "e_step_importance_sampling", weights = sampled_data_uninformative[[k]]$weights)
# out <- wglmm( bin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#               data = subset(sampled_data_informative[[k]], select = -c(inclusion_prob, weights)), trace = TRUE, iter1 = 250,
#               iter2 = 1001, MI = 2000, family = binomial(), e_step_method = "e_step_importance_sampling", weights = sampled_data_informative[[k]]$weights)
# # GIBBS LIN
# out <- wglmm( lin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#               data = sampled_data_uninformative[[k]], trace = TRUE, iter1 = 250,
#               iter2 = 1001, MI = 2000, family = gaussian(), e_step_method = "e_step_gibbs_sampling")
# out <- wglmm( lin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#               data = sampled_data_informative[[k]], trace = TRUE, iter1 = 250,
#               iter2 = 1001, MI = 2000, family = gaussian(), e_step_method = "e_step_gibbs_sampling")
# # GIBBS BIN
# out <- wglmm( bin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#               data = sampled_data_uninformative[[k]], trace = TRUE, iter1 = 250,
#               iter2 = 1001, MI = 2000, family = binomial(), e_step_method = "e_step_gibbs_sampling")
# out <- wglmm( bin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#               data = sampled_data_informative[[k]], trace = TRUE, iter1 = 250,
#               iter2 = 1001, MI = 2000, family = binomial(), e_step_method = "e_step_gibbs_sampling")

### For plotting the box-plots (IMPORTANCE) ###
## Non-informative ##
# Extracting coefficients and creating a tidy data frame
df_coefs <- do.call(rbind, lapply(seq_along(results_list_uninf), function(k) {
  tibble(
    iteration = k,
    coef1 = results_list_uninf[[k]]$coef[1],
    coef2 = results_list_uninf[[k]]$coef[2],
    coef3 = results_list_uninf[[k]]$coef[3]
  )
}))

# Reshape data from wide to long format
df_long <- pivot_longer(df_coefs, cols = starts_with("coef"), names_to = "Coefficient", values_to = "Value")

df_long$Coefficient <- factor(df_long$Coefficient,
                              levels = c("coef1", "coef2", "coef3"),
                              labels = c(expression(beta[0]), expression(beta[1]), expression(beta[2])))

# Create the boxplot with specified aesthetic changes
p <- ggplot(df_long, aes(x = Coefficient, y = Value)) +
  geom_boxplot(outlier.colour = "gray", outlier.shape = 1, outlier.size = 2) +
  annotate("segment", x = 0.6, xend = 1.4, y = 4, yend = 4, color = "red", size = 1) +  # Line for coef1
  annotate("segment", x = 1.6, xend = 2.4, y = -2, yend = -2, color = "red", size = 1) +  # Line for coef2
  annotate("segment", x = 2.6, xend = 3.4, y = -1, yend = -1, color = "red", size = 1) +  # Line for coef3
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = NULL) +  # Remove y-axis label but keep the ticks
  scale_x_discrete(labels = c(expression(beta[0]), expression(beta[1]), expression(beta[2])))

print(p)
ggsave("lin_imp_non_informative_boxplot_coeff.png", plot = p, width = 10, height = 6, dpi = 300)

# Extracting variance and creating a tidy data frame
df_vars <- do.call(rbind, lapply(seq_along(results_list_uninf), function(k) {
  tibble(
    iteration = k,
    var1 = results_list_uninf[[k]]$VarCov$group1[1,1],
    var2 = results_list_uninf[[k]]$VarCov$group2[1,1],
    var3 = results_list_uninf[[k]]$VarCov$group2[1,2],
    var4 = results_list_uninf[[k]]$VarCov$group2[2,2]
  )
}))

# Reshape data from wide to long format
df_long <- pivot_longer(df_vars, cols = starts_with("var"), names_to = "Variance", values_to = "Value")

df_long$Variance <- factor(df_long$Variance,
                              levels = c("var1", "var2", "var3", "var4"),
                              labels = c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4])))

# Create the boxplot with specified aesthetic changes
p <- ggplot(df_long, aes(x = Variance, y = Value)) +
  geom_boxplot(outlier.colour = "gray", outlier.shape = 1, outlier.size = 2) +
  annotate("segment", x = 0.6, xend = 1.4, y = 4, yend = 4, color = "red", size = 1) +  # Line for var1
  annotate("segment", x = 1.6, xend = 2.4, y = 0.7, yend = 0.7, color = "red", size = 1) +  # Line for var2
  annotate("segment", x = 2.6, xend = 3.4, y = 0.5, yend = 0.5, color = "red", size = 1) +  # Line for var3
  annotate("segment", x = 3.6, xend = 4.4, y = 1.3, yend = 1.3, color = "red", size = 1) +  # Line for var4
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = NULL) +  # Remove y-axis label but keep the ticks
  scale_x_discrete(labels = c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4])))

print(p)
ggsave("lin_imp_non_informative_boxplot_var.png", plot = p, width = 10, height = 6, dpi = 300)

## Informative ##
# Extracting coefficients and creating a tidy data frame
df_coefs <- do.call(rbind, lapply(seq_along(results_list_inf), function(k) {
  tibble(
    iteration = k,
    coef1 = results_list_inf[[k]]$coef[1],
    coef2 = results_list_inf[[k]]$coef[2],
    coef3 = results_list_inf[[k]]$coef[3]
  )
}))

# Reshape data from wide to long format
df_long <- pivot_longer(df_coefs, cols = starts_with("coef"), names_to = "Coefficient", values_to = "Value")

df_long$Coefficient <- factor(df_long$Coefficient,
                              levels = c("coef1", "coef2", "coef3"),
                              labels = c(expression(beta[0]), expression(beta[1]), expression(beta[2])))

# Create the boxplot with specified aesthetic changes
p <- ggplot(df_long, aes(x = Coefficient, y = Value)) +
  geom_boxplot(outlier.colour = "gray", outlier.shape = 1, outlier.size = 2) +
  annotate("segment", x = 0.6, xend = 1.4, y = 4, yend = 4, color = "red", size = 1) +  # Line for coef1
  annotate("segment", x = 1.6, xend = 2.4, y = -2, yend = -2, color = "red", size = 1) +  # Line for coef2
  annotate("segment", x = 2.6, xend = 3.4, y = -1, yend = -1, color = "red", size = 1) +  # Line for coef3
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = NULL) +  # Remove y-axis label but keep the ticks
  scale_x_discrete(labels = c(expression(beta[0]), expression(beta[1]), expression(beta[2])))

print(p)
ggsave("lin_imp_informative_boxplot_coeff.png", plot = p, width = 10, height = 6, dpi = 300)

# Extracting variance and creating a tidy data frame
df_vars <- do.call(rbind, lapply(seq_along(results_list_inf), function(k) {
  tibble(
    iteration = k,
    var1 = results_list_inf[[k]]$VarCov$group1[1,1],
    var2 = results_list_inf[[k]]$VarCov$group2[1,1],
    var3 = results_list_inf[[k]]$VarCov$group2[1,2],
    var4 = results_list_inf[[k]]$VarCov$group2[2,2]
  )
}))

# Reshape data from wide to long format
df_long <- pivot_longer(df_vars, cols = starts_with("var"), names_to = "Variance", values_to = "Value")

df_long$Variance <- factor(df_long$Variance,
                           levels = c("var1", "var2", "var3", "var4"),
                           labels = c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4])))

# Create the boxplot with specified aesthetic changes
p <- ggplot(df_long, aes(x = Variance, y = Value)) +
  geom_boxplot(outlier.colour = "gray", outlier.shape = 1, outlier.size = 2) +
  annotate("segment", x = 0.6, xend = 1.4, y = 4, yend = 4, color = "red", size = 1) +  # Line for var1
  annotate("segment", x = 1.6, xend = 2.4, y = 0.7, yend = 0.7, color = "red", size = 1) +  # Line for var2
  annotate("segment", x = 2.6, xend = 3.4, y = 0.5, yend = 0.5, color = "red", size = 1) +  # Line for var3
  annotate("segment", x = 3.6, xend = 4.4, y = 1.3, yend = 1.3, color = "red", size = 1) +  # Line for var4
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = NULL) +  # Remove y-axis label but keep the ticks
  scale_x_discrete(labels = c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4])))

print(p)
ggsave("lin_imp_informative_boxplot_var.png", plot = p, width = 10, height = 6, dpi = 300)

### For printing means ###
mean_time_stat_uninf <- mean(sapply(time_stats_uninf, function(x) x[[3]]), na.rm = TRUE)
print(mean_time_stat_uninf)
mean_memory_stat_uninf <- mean(sapply(memory_stats_uninf, function(x) x[[1]]), na.rm = TRUE)
print(mean_memory_stat_uninf)
mean_cpu_usage_stat_uninf <- mean(sapply(cpu_usage_stats_uninf, function(x) x[[3]]), na.rm = TRUE)
print(mean_cpu_usage_stat_uninf)

mean_time_stat_inf <- mean(sapply(time_stats_inf, function(x) x[[3]]), na.rm = TRUE)
print(mean_time_stat_inf)
mean_memory_stat_inf <- mean(sapply(memory_stats_inf, function(x) x[[1]]), na.rm = TRUE)
print(mean_memory_stat_inf)
mean_cpu_usage_stat_inf <- mean(sapply(cpu_usage_stats_inf, function(x) x[[3]]), na.rm = TRUE)
print(mean_cpu_usage_stat_inf)


### For plotting the box-plots (GIBBS) ###
## Non-informative ##
# Extracting coefficients and creating a tidy data frame
df_coefs <- do.call(rbind, lapply(seq_along(results_list_uninf), function(k) {
  tibble(
    iteration = k,
    coef1 = results_list_uninf[[k]]$solutions[1,1],
    coef2 = results_list_uninf[[k]]$solutions[2,1],
    coef3 = results_list_uninf[[k]]$solutions[3,1]
  )
}))

# Reshape data from wide to long format
df_long <- pivot_longer(df_coefs, cols = starts_with("coef"), names_to = "Coefficient", values_to = "Value")

df_long$Coefficient <- factor(df_long$Coefficient,
                              levels = c("coef1", "coef2", "coef3"),
                              labels = c(expression(beta[0]), expression(beta[1]), expression(beta[2])))

# Create the boxplot with specified aesthetic changes
p <- ggplot(df_long, aes(x = Coefficient, y = Value)) +
  geom_boxplot(outlier.colour = "gray", outlier.shape = 1, outlier.size = 2) +
  annotate("segment", x = 0.6, xend = 1.4, y = 4, yend = 4, color = "red", size = 1) +  # Line for coef1
  annotate("segment", x = 1.6, xend = 2.4, y = -2, yend = -2, color = "red", size = 1) +  # Line for coef2
  annotate("segment", x = 2.6, xend = 3.4, y = -1, yend = -1, color = "red", size = 1) +  # Line for coef3
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = NULL) +  # Remove y-axis label but keep the ticks
  scale_x_discrete(labels = c(expression(beta[0]), expression(beta[1]), expression(beta[2])))

print(p)
ggsave("bin_gibbs_non_informative_boxplot_coeff.png", plot = p, width = 10, height = 6, dpi = 300)

# Extracting variance and creating a tidy data frame
df_vars <- do.call(rbind, lapply(seq_along(results_list_uninf), function(k) {
  tibble(
    iteration = k,
    var1 = results_list_uninf[[k]]$Gcovariances[1,1],
    var2 = results_list_uninf[[k]]$Gcovariances[2,1],
    var3 = results_list_uninf[[k]]$Gcovariances[3,1],
    var4 = results_list_uninf[[k]]$Gcovariances[5,1]
  )
}))

# Reshape data from wide to long format
df_long <- pivot_longer(df_vars, cols = starts_with("var"), names_to = "Variance", values_to = "Value")

df_long$Variance <- factor(df_long$Variance,
                           levels = c("var1", "var2", "var3", "var4"),
                           labels = c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4])))

# Create the boxplot with specified aesthetic changes
p <- ggplot(df_long, aes(x = Variance, y = Value)) +
  geom_boxplot(outlier.colour = "gray", outlier.shape = 1, outlier.size = 2) +
  annotate("segment", x = 0.6, xend = 1.4, y = 4, yend = 4, color = "red", size = 1) +  # Line for var1
  annotate("segment", x = 1.6, xend = 2.4, y = 0.7, yend = 0.7, color = "red", size = 1) +  # Line for var2
  annotate("segment", x = 2.6, xend = 3.4, y = 0.5, yend = 0.5, color = "red", size = 1) +  # Line for var3
  annotate("segment", x = 3.6, xend = 4.4, y = 1.3, yend = 1.3, color = "red", size = 1) +  # Line for var4
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = NULL) +  # Remove y-axis label but keep the ticks
  scale_x_discrete(labels = c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4])))

print(p)
ggsave("bin_gibbs_non_informative_boxplot_var.png", plot = p, width = 10, height = 6, dpi = 300)

## Informative ##
# Extracting coefficients and creating a tidy data frame
df_coefs <- do.call(rbind, lapply(seq_along(results_list_inf), function(k) {
  tibble(
    iteration = k,
    coef1 = results_list_inf[[k]]$solutions[1,1],
    coef2 = results_list_inf[[k]]$solutions[2,1],
    coef3 = results_list_inf[[k]]$solutions[3,1]
  )
}))

# Reshape data from wide to long format
df_long <- pivot_longer(df_coefs, cols = starts_with("coef"), names_to = "Coefficient", values_to = "Value")

df_long$Coefficient <- factor(df_long$Coefficient,
                              levels = c("coef1", "coef2", "coef3"),
                              labels = c(expression(beta[0]), expression(beta[1]), expression(beta[2])))

# Create the boxplot with specified aesthetic changes
p <- ggplot(df_long, aes(x = Coefficient, y = Value)) +
  geom_boxplot(outlier.colour = "gray", outlier.shape = 1, outlier.size = 2) +
  annotate("segment", x = 0.6, xend = 1.4, y = 4, yend = 4, color = "red", size = 1) +  # Line for coef1
  annotate("segment", x = 1.6, xend = 2.4, y = -2, yend = -2, color = "red", size = 1) +  # Line for coef2
  annotate("segment", x = 2.6, xend = 3.4, y = -1, yend = -1, color = "red", size = 1) +  # Line for coef3
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = NULL) +  # Remove y-axis label but keep the ticks
  scale_x_discrete(labels = c(expression(beta[0]), expression(beta[1]), expression(beta[2])))

print(p)
ggsave("bin_gibbs_informative_boxplot_coeff.png", plot = p, width = 10, height = 6, dpi = 300)

# Extracting variance and creating a tidy data frame
df_vars <- do.call(rbind, lapply(seq_along(results_list_inf), function(k) {
  tibble(
    iteration = k,
    var1 = results_list_inf[[k]]$Gcovariances[1,1],
    var2 = results_list_inf[[k]]$Gcovariances[2,1],
    var3 = results_list_inf[[k]]$Gcovariances[3,1],
    var4 = results_list_inf[[k]]$Gcovariances[5,1]
  )
}))

# Reshape data from wide to long format
df_long <- pivot_longer(df_vars, cols = starts_with("var"), names_to = "Variance", values_to = "Value")

df_long$Variance <- factor(df_long$Variance,
                           levels = c("var1", "var2", "var3", "var4"),
                           labels = c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4])))

# Create the boxplot with specified aesthetic changes
p <- ggplot(df_long, aes(x = Variance, y = Value)) +
  geom_boxplot(outlier.colour = "gray", outlier.shape = 1, outlier.size = 2) +
  annotate("segment", x = 0.6, xend = 1.4, y = 4, yend = 4, color = "red", size = 1) +  # Line for var1
  annotate("segment", x = 1.6, xend = 2.4, y = 0.7, yend = 0.7, color = "red", size = 1) +  # Line for var2
  annotate("segment", x = 2.6, xend = 3.4, y = 0.5, yend = 0.5, color = "red", size = 1) +  # Line for var3
  annotate("segment", x = 3.6, xend = 4.4, y = 1.3, yend = 1.3, color = "red", size = 1) +  # Line for var4
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = NULL) +  # Remove y-axis label but keep the ticks
  scale_x_discrete(labels = c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4])))

print(p)
ggsave("bin_gibbs_informative_boxplot_var.png", plot = p, width = 10, height = 6, dpi = 300)
