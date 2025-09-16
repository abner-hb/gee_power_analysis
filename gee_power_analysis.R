# Setup ####
library(mvtnorm) # to simulate from multivariate normal
library(geepack) # to fit GEE regression
library(broom) # For easy extraction of confidence intervals
source("aux_functions.R") # load functions to simulate data
set.seed(2025) # to increase reproducibility
# Example of a single simulation for a hypothesis test ####
# Simulate covariates and outcomes from model
data <- simulate_full_data(
  llfs_obs = 2500,
  hrs_obs = 200,
  beta_coefs = c(0, 1, 0.5, 2), # intercept, age_stand, llfs, age_stand:llfs
  my_sigma = 1.3, # standard deviation of residuals
  cor_coef = 0.8 # within-cluster correlation coefficient
)
# Fit the model
model <- geeglm(
  y_sim ~ age_stand + llfs + age_stand:llfs,
  data = data,
  id = cluster_id,
  family = gaussian,
  corstr = "exchangeable"
)
# Check model results to verify the fitting works
summary(model)
# Set alpha for confidence level
alpha <- 0.05
# Extract lower limit of confidence interval
results_mod <- tidy(model, conf.int = TRUE, conf.level = 1 - alpha)
lower_lim <- as.double(
  results_mod[results_mod$term == "age_stand:llfs", "conf.low"]
)
# Is the lower limit greater than zero?
test_result <- lower_lim > 0
test_result
# Run multiple simulations to estimate power ####
# The function below runs multiple simulations in the background and returns
# the estimated power
my_power <- estimate_power(
  llfs_obs = 2500,
  hrs_obs = 200,
  beta_coefs = c(0, 1, 0.5, 0.2), # intercept, age_stand, llfs, age_stand:llfs
  my_sigma = 1.3, # standard deviation of residuals
  cor_coef = 0.1, # within-cluster correlation coefficient
  alpha = 0.05,
  num_sims = 100 # For quick results. Increase later for more accuracy
)
my_power
