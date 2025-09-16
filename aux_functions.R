simulate_covariates <- function(llfs_obs, hrs_obs, num_clusters = 2500 / 9.4) {
  # Calculate total number of observations across both study types
  total_obs <- llfs_obs + hrs_obs
  # Generate standardized age values from normal distribution for all
  # observations
  age_stand_vec <- rnorm(total_obs)
  # Create binary indicator: 1 for LLFS participants, 0 for HRS participants
  llfs_vec <- c(rep(1L, llfs_obs), rep(0L, hrs_obs))
  # Generate cluster IDs for LLFS participants (sampling clusters with
  # replacement)
  cluster_id_llfs <- paste(
    "llfs",
    sample(1:num_clusters, size = llfs_obs, replace = TRUE),
    sep = "_"
  )
  # Generate sequential cluster IDs for HRS participants
  cluster_id_hrs <- paste(
    "hrs",
    1:hrs_obs,
    sep = "_"
  )
  # Combine LLFS and HRS cluster IDs into single vector
  cluster_id_vec <- c(cluster_id_llfs, cluster_id_hrs)
  # Create data frame with all simulated covariates
  covariates_df <- data.frame(
    age_stand = age_stand_vec,
    llfs = llfs_vec,
    cluster_id = as.factor(cluster_id_vec)
  )
  # Sort by cluster ID
  covariates_df <- covariates_df[order(covariates_df$cluster_id), ]
  return(covariates_df)
}

compute_lin_pred <- function(data, beta_coefs) {
  # Create design matrix
  # Assuming intercept + age_stand + llfs
  X <- model.matrix(~ age_stand + llfs + age_stand:llfs, data = data)
  # Error handling: Check dimension compatibility
  if (ncol(X) != length(beta_coefs)) {
    stop(paste0(
      "Dimension mismatch: Design matrix X has ",
      ncol(X),
      " columns but beta_coefs has ",
      length(beta_coefs),
      " elements. These must be equal for matrix multiplication."
    ))
  }
  # Additional check: ensure beta_coefs is numeric
  if (!is.numeric(beta_coefs)) {
    stop("beta_coefs must be a numeric vector")
  }
  # Calculate linear predictor (mean structure)
  linear_pred <- X %*% beta_coefs
  return(linear_pred)
}

simulate_cor_errs <- function(stand_dev, cor_coef, dim_size) {
  # Error handling: Check that correlation coefficient is valid
  if (!is.numeric(cor_coef) || length(cor_coef) != 1) {
    stop("cor_coef must be a single numeric value")
  }
  if (cor_coef < -1 || cor_coef > 1) {
    stop(paste0(
      "cor_coef must be in the interval [-1, 1]. ",
      "Provided value: ",
      cor_coef
    ))
  }
  # Exchangeable correlation matrix
  correl_matrix <- matrix(cor_coef, dim_size, dim_size)
  diag(correl_matrix) <- 1
  # Covariance matrix
  covar_matrix <- stand_dev^2 * correl_matrix
  # Simulate correlated errors
  correlated_errors <- rmvnorm(
    n = 1,
    mean = rep(0, dim_size),
    sigma = covar_matrix
  )
  return(as.vector(correlated_errors))
}

simulate_gee_outcome <- function(data, beta_coefs, my_sigma, cor_coef) {
  # Calculate linear predictor (mean structure)
  linear_pred_all <- compute_lin_pred(data = data, beta_coefs = beta_coefs)
  # Initialize the outcome vector
  y_sim <- numeric(nrow(data))
  # Get cluster information
  cluster_ids <- as.numeric(data$cluster_id)
  unique_clusters <- unique(cluster_ids)
  # Simulate errors for each cluster
  for (i in seq_along(unique_clusters)) {
    cluster_idx <- which(cluster_ids == unique_clusters[i])
    cluster_size <- length(cluster_idx)
    # Simulate errors using an exchangeable correlation matrix
    cluster_errors <- simulate_cor_errs(
      stand_dev = my_sigma,
      cor_coef = cor_coef,
      dim_size = cluster_size
    )
    # Add cluster means with their errors
    cluster_means <- linear_pred_all[cluster_idx]
    y_sim[cluster_idx] <- cluster_means + cluster_errors
  }
  return(y_sim)
}

simulate_full_data <- function(
  llfs_obs,
  hrs_obs,
  beta_coefs,
  my_sigma,
  cor_coef
) {
  data <- simulate_covariates(llfs_obs = llfs_obs, hrs_obs = hrs_obs)
  data$y_sim <- simulate_gee_outcome(
    data = data,
    beta_coefs = beta_coefs, # intercept, age_stand, llfs, age_stand:llfs
    my_sigma = my_sigma, # standard deviation of residuals
    cor_coef = cor_coef # within-cluster correlation coefficient
  )
  return(data)
}

simulate_hyp_test <- function(
  llfs_obs,
  hrs_obs,
  beta_coefs,
  my_sigma,
  cor_coef,
  alpha = 0.05
) {
  # Simulate covariates and outcomes from model
  data <- simulate_full_data(
    llfs_obs = llfs_obs,
    hrs_obs = hrs_obs,
    beta_coefs = beta_coefs, # intercept, age_stand, llfs, age_stand:llfs
    my_sigma = my_sigma, # standard deviation of residuals
    cor_coef = cor_coef # within-cluster correlation coefficient
  )
  # Fit the model
  model <- geeglm(
    y_sim ~ age_stand + llfs + age_stand:llfs,
    data = data,
    id = cluster_id,
    family = gaussian,
    corstr = "exchangeable"
  )
  # Extract lower limit of confidence interval
  results_mod <- tidy(model, conf.int = TRUE, conf.level = 1 - alpha)
  lower_lim <- as.double(
    results_mod[results_mod$term == "age_stand:llfs", "conf.low"]
  )
  # Return whether the lower limit is above zero
  return(lower_lim > 0)
}

estimate_power <- function(
  llfs_obs,
  hrs_obs,
  beta_coefs,
  my_sigma,
  cor_coef,
  alpha = 0.05,
  num_sims = 1000
) {
  tests_results <- replicate(
    n = num_sims,
    expr = simulate_hyp_test(
      llfs_obs = llfs_obs,
      hrs_obs = hrs_obs,
      beta_coefs = beta_coefs, # intercept, age_stand, llfs, age_stand:llfs
      my_sigma = my_sigma, # standard deviation of residuals
      cor_coef = cor_coef # within-cluster correlation coefficient
    ),
    simplify = "vector"
  )
  power <- mean(tests_results)
  return(power)
}
