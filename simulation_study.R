library(tidyverse)
library(refund)
library(mgcv)

# Set up -----------------------------------------------------------------------

# Root folder of project
setwd("/users/9/pai00032/step_up_sim")

dat <- readRDS("data/dat.RDS")
source("code/fosr_fpcareg.R")

steps_data <- dat$steps_w[, 1:24]
all_na_int <- rowSums(!is.na(steps_data)) <= 3
steps_f <- steps_data[!all_na_int, ] # Remove rows with few values
fpca_res_int <- fpca.face(unclass(steps_f), knots = 15, pve = .95)
steps_data <- dat$steps_w[, 25:36]
all_na_post <- rowSums(!is.na(steps_data)) <= 3

# Helper function
unAsIs <- function(X) {
  if("AsIs" %in% class(X)) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}

full_data_res_int <- fosr_coef_function(dat[!all_na_int, ], "int")

# Functions --------------------------------------------------------------------

# IMPT: Simulation is only for intervention period

# Data that's not all NA from the full data
dat_not_na <- dat[!params[["int"]][["all_na"]], ]

## Get the fixed effects from the real data
arm_coef_df <- full_data_res_int[[1]]
# Matrix with cols collab, comp, supp
arm_coef_mat <- arm_coef_df %>%
  pivot_wider(names_from = arm, values_from = value) %>%
  select(-week) %>%
  as.matrix()
intercept <- full_data_res_int[[2]]
baseline_steps <- full_data_res_int[[3]]
covar_coef <- full_data_res_int[[8]]
coef_mat <- cbind(intercept, arm_coef_mat, baseline_steps) # 24 x 5

## Get error variance
error_var <- full_data_res_int[[5]]$sig2

## Get eigenfunctions/values
ef <- fpca_res_int$efunctions # 24 x 2
eval <- fpca_res_int$evalues

generate_data <- function(N) {
  # Get samples from the data
  sample <- dat_not_na %>% sample_n(N, replace = T)
  
  # Assemble fixed effects
  X_mat <- model.matrix(~ 1 + arm + baseline_steps + start_day + 
                          ESE_sticking_to_it + age + gender
                          , data = sample) # 602 x 5
  fixed <- X_mat %*% t(cbind(coef_mat, matrix(rep(covar_coef, each = 24), nrow = 24)))
  
  ## Assemble random effects
  xi <- sapply(eval, function(ev) {
    rnorm(nrow(sample), mean = 0, sd = sqrt(ev))
  }) # 602 x 2
  random <- xi %*% t(ef)
  
  # Errors
  error <- rnorm(nrow(sample) * nrow(ef), mean = 0, sd = sqrt(error_var))
  
  # Assemble
  W <- fixed + random + error
  
  # Sample missingness patterns
  miss_patterns <- dat_not_na %>% 
    sample_n(N, replace = T) %>%
    select(steps_w_int) %>% 
    unAsIs() %>% 
    is.na()
  W <- ifelse(miss_patterns, NA, W)
  
  # Put in form of original df so that we can use same FoSR function
  sample$steps_w[, 1:24] <- W
  sample
}

ise <- function(est, truth) {
  # IMPT: It's set up for intervention period
  sfsmisc::integrate.xy(x = 1:24, fx = (est - truth)^2)
}

# Run simulation ---------------------------------------------------------------

intermed_res_folder <- "simulation_results/take_cov_2"
intermed_res_folder_fpca <- "simulation_results/take_cov_2_fpca"

N_options <- c(600, 900, 1200, 1500, 1800, 2100, 2400)
job_number <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # Task number in job

for (N in N_options) {
  # Set seed so that every simulation is on a different seed
  set.seed(job_number * 50 + N / 100)
  
  # Generate data
  sample <- generate_data(N)
  
  #Fit FoSR model
  res <- fosr_coef_function(sample, "int", npc = 2)
  arm_coef_mat_hat <- res[[1]] %>%
    pivot_wider(names_from = arm, values_from = value) %>%
    select(-week) %>%
    as.matrix()
  intercept_hat <- res[[2]]
  baseline_steps_hat <- res[[3]]
  coef_mat_hat <- cbind("week" = 1:24, "intercept" = intercept_hat, arm_coef_mat_hat,
                        "baseline_steps" = baseline_steps_hat,
                        "N" = N,
                        "job_number" = job_number,
                        "seed" = job_number * 50 + N / 100)
  covar_coef_hat <- c(N = N,
                      job_number = job_number,
                      seed = job_number * 50 + N / 100,
                      res[[8]])

  saveRDS(coef_mat_hat, paste0(intermed_res_folder, "/", N, "-", job_number, ".RDS"))
  saveRDS(covar_coef_hat, paste0(intermed_res_folder,"/", N, "-", job_number, "-covar.RDS"))
  
  # Fit FPCA + regression model
  res <- fpca_coef_function(sample, "int", npc = 2)
  arm_coef_mat_hat <- res[[1]] %>%
    pivot_wider(names_from = arm, values_from = value) %>%
    select(-week) %>%
    as.matrix()
  intercept_hat <- res[[2]]
  baseline_steps_hat <- res[[3]]
  coef_mat_hat <- cbind("week" = 1:24, "intercept" = intercept_hat, arm_coef_mat_hat, 
                        "baseline_steps" = baseline_steps_hat,
                        "N" = N,
                        "job_number" = job_number,
                        "seed" = job_number * 50 + N / 100)
  covar_coef_hat <- cbind(N = N,
                          job_number = job_number,
                          seed = job_number * 50 + N / 100,
                          res[[8]])
  
  saveRDS(coef_mat_hat, paste0(intermed_res_folder_fpca, "/", N, "-", job_number, ".RDS"))
  saveRDS(covar_coef_hat, paste0(intermed_res_folder_fpca,"/", N, "-", job_number, "-covar.RDS"))
}


# Helpful code to analyze results -----------------------------------------------

# Calculate ISE of fixed effects
fixed_ise <- sapply(1:(S * length(N_options) * 2), function(s) {
  s_res <- sim_res[((s - 1) * 24 + 1):(s * 24), ]
  c(sapply(1:ncol(coef_mat), function(k) {
    ise(s_res[, k + 1], coef_mat[, k])
  }), s_res[1, "job_number"], s_res[1, "N"], s_res[1, "method"])
}) %>% t()

