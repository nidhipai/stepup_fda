# Parameters -------------------------------------------------------------------

# Parameters for FPCA + reg and FoSR for intervention and follow-up periods
params <- list(
  int = list(start = 1, end = 24, fpca_k = 15, fosr_k = 20, all_na = all_na_int),
  post = list(start = 25, end = 36, fpca_k = 8, fosr_k = 7, all_na = all_na_post)
)

# Fit FoSR ---------------------------------------------------------------------
# Written as a function for both periods and subsampled/simulated data

fosr_coef_function <- function(dat_sample, period, npc = 2) {
  start <- params[[period]][["start"]]
  end <- params[[period]][["end"]]
  
  ef_names <- paste0("ef_", 1:npc)
  
  # Fit FPCA model
  steps_data <- dat_sample$steps_w[, start:end]
  fpca_res <- fpca.face(unclass(steps_data), knots = params[[period]][["fpca_k"]], npc = npc)
  fpca_ef <- data.frame(fpca_res$efunctions)
  colnames(fpca_ef) <- ef_names
  fpca_ef$week <- start:end
  covariates <- c("start_day", "ESE_sticking_to_it", "age", "gender")
  covariates_formula <- paste0(covariates, collapse = "+")
  fosr_df <- dat_sample %>%
    mutate(re_id = factor(1:n())) %>%
    select(id, re_id, arm, baseline_steps, covariates) %>% # Necessary to unnest
    cbind(unAsIs(dat_sample$steps_w[, start:end])) %>%
    pivot_longer(cols = -c(id, re_id, arm, baseline_steps, covariates), names_to = "week", values_to = "log_steps") %>%
    mutate(week = as.integer(week)) %>%
    mutate(week = week + start - 1) %>% # Colnames are indexes
    left_join(fpca_ef[, c("week", ef_names)])
  
  # Fit FoSR model
  ef_formula <- paste0("s(re_id, bs = \"re\", by = ", ef_names, ")", collapse = " + ")
  fosr_formula <- paste0("log_steps ~ arm + s(week, by = arm, k = params[[period]][[\"fosr_k\"]], bs = \"ps\") + s(week, by = baseline_steps, bs = \"ps\") + ",
                         covariates_formula, "+", ef_formula)
  fosr <- bam(as.formula(fosr_formula), method = "fREML", data = fosr_df)
  summary_fosr <- summary(fosr)
  
  # Extract intercept
  predict_res <- fosr_df %>%
    mutate(arm = "Control", baseline_steps = 1) %>%
    slice(1:(end - start + 1)) %>% # First set of start to end
    predict(fosr, newdata = ., type = "terms")
  p_intercept <- summary_fosr$p.coeff["(Intercept)"]
  intercept <- p_intercept + predict_res[, "s(week):armControl"]
  # Extract coef for baseline steps
  baseline_coef <- predict_res[, "s(week):baseline_steps"]
  
  # Extract coefficients for arm
  control_coef <- predict_res[, "s(week):armControl"]
  arm_coef_df <- data.frame()
  p_coef <- list()
  for (curr_arm in levels(dat$arm)) {
    if (curr_arm == "Control") next
    predict_res <- fosr_df %>%
      mutate(arm = curr_arm) %>%
      slice(1:(end - start + 1)) %>%
      predict(fosr, newdata = ., type = "terms", se.fit = F)
    parametric_coef <- summary_fosr$p.coeff[paste0("arm", curr_arm)]
    arm_df <- data.frame(value = predict_res[, paste0("s(week):arm", curr_arm)] +
                           parametric_coef - control_coef) %>%
      mutate(arm = curr_arm, week = start:end)
    arm_coef_df <- rbind(arm_coef_df, arm_df)
    p_coef[[curr_arm]] <- parametric_coef
  }
  covariates_dummy <- c("start_day", "ESE_sticking_to_it", "age", "genderMale")
  covar_coef <- summary_fosr$p.coeff[covariates_dummy]

  return(list(arm_coef_df, intercept, baseline_coef, fosr_df, fosr, p_intercept, p_coef, covar_coef))
}


# FPCA + regression as coefficient function --------------------------------------------------

fpca_coef_function <- function(dat_sample, period, npc = 2) {
  start <- params[[period]][["start"]]
  end <- params[[period]][["end"]]
  
  # Fit FPCA model
  steps_data <- dat_sample$steps_w[, start:end]
  fpca_res <- fpca.face(unclass(steps_data), knots = params[[period]][["fpca_k"]], npc = npc)
  
  num_ef <- fpca_res$npc
  scores <- fpca_res$scores
  efuncs <- fpca_res$efunctions
  
  # Construct beta matrix
  beta <- matrix(nrow = num_ef, ncol = 9) # Intercept + 3 arms + baseline_steps
  for (ef in 1:num_ef) { # Get regression coefficients for covariates and efs
    dat_sample$scores <- scores[, ef]
    lm_res <- lm(scores ~ arm + baseline_steps + start_day + ESE_sticking_to_it + age + gender, data = dat_sample) 
    beta[ef, ] <- coef(lm_res)
  }
  beta_pr <- efuncs %*% beta
  arm_beta <- data.frame(beta_pr)
  colnames(arm_beta) <- c("Intercept", "Collaborative", "Competitive", "Supportive",
                          "Baseline", "start_day", "ESE_sticking_to_it", "age", "gender")
  intercept <- arm_beta[, "Intercept"] + fpca_res$mu 
  baseline_steps <- arm_beta[, "Baseline"]
  covar_coef <- arm_beta[, 6:9]
  arm_beta <- arm_beta[, 2:4]
  
  # Long format
  arm_beta <- arm_beta %>%
    mutate(week = start:end) %>%
    pivot_longer(cols = -week, names_to = "arm", values_to = "value") %>%
    select(value, arm, week) # Reorder columns
  
  list(arm_beta, intercept, baseline_steps, NA, NA, NA, NA, covar_coef)
}


# Bootstrap --------------------------------------------------------------------

# Generic function for bootstrapping coefficients
bootstrap <- function(B, period, method, npc = 2) {
  coef_function <- ifelse(method == "fosr", fosr_coef_function, fpca_coef_function)
  dat_not_na <- dat[!params[[period]][["all_na"]], ]
  start <- params[[period]][["start"]]
  end <- params[[period]][["end"]]
  
  period_length <- end - start + 1
  arm_bs_res <- data.frame(matrix(nrow = B * 3 * period_length, ncol = 4))
  intercept_bs_res <- data.frame(matrix(nrow = B * period_length, ncol = 3))
  baseline_coef_bs_res <- data.frame(matrix(nrow = B * period_length, ncol = 3))
  if (method == "fosr") {
    covar_coef_res <- data.frame(matrix(nrow = B, ncol = 5)) # 4 covar
  } else {
    covar_coef_res <- data.frame(matrix(nrow = B * period_length, ncol = 6))
  }
  l_ <- l <- 1 # index counter
  for (b in 1:B) {
    print(b)
    # Bootstrap sample
    dat_sample <- sample_n(dat_not_na, size = nrow(dat_not_na), replace = T)
    b_res <- do.call(coef_function, list(dat_sample, period, npc))
    arm_bs_res[l:(l + 3 * period_length - 1), ] <- b_res[[1]] %>% mutate(b = b)
    intercept_bs_res[l_:(l_ + period_length - 1), ] <- data.frame(week = start:end, 
                                                                  value = b_res[[2]], b = b)
    baseline_coef_bs_res[l_:(l_ + period_length - 1), ] <- data.frame(week = start:end, 
                                                                      value = b_res[[3]], b = b)
    if (method == "fosr") {
      covar_coef_res[b, ] <- c(b, b_res[[8]])
    } else {
      covar_coef_res[l_:(l_ + period_length - 1), ] <- cbind(week = start:end, b_res[[8]], b = b)
    }
    l <- l + 3 * period_length
    l_ <- l_ + period_length
  }
  colnames(arm_bs_res) <- c("value", "arm", "week", "b")
  colnames(intercept_bs_res) <- c("week", "value", "b")
  colnames(baseline_coef_bs_res) <- c("week", "value", "b")
  if (method == "fosr") {
    colnames(covar_coef_res) <- c("b", names(b_res[[8]]))
  } else {
    colnames(covar_coef_res) <- c("week", colnames(b_res[[8]]), "b")
  }
  list(arm_bs_res, intercept_bs_res, baseline_coef_bs_res, covar_coef_res)
}


# Plots and CMA calculations -------------------------------------------------------------

# Plots relating to coefficients for arm
bootstrap_plots_arm <- function(period, arm_bs_res, full_arm_coef, method, plot_suffix = "") {
  start <- params[[period]][["start"]]
  end <- params[[period]][["end"]]
  
  # Coefficients and SE from BS ------------------------------------------------
  se_res <- arm_bs_res %>%
    group_by(arm, week) %>%
    summarize(se = sd(value)) %>%
    inner_join(full_arm_coef, by = c("arm", "week")) %>%
    rename(full_value = value) # This is value computed with all data
  
  # Unadjusted Wald at each point
  unadjusted <- se_res %>% mutate(q = 1.96)
  # CMA adjusted intervals
  cma <- arm_bs_res %>% # First get q
    left_join(se_res, by = c("arm", "week")) %>%
    group_by(arm, b) %>%
    summarize(db = max(abs(value - full_value) / se)) %>%
    group_by(arm) %>%
    summarize(q = quantile(db, .95)) %>% 
    right_join(se_res, by = "arm")
  
  # Plot of unadjusted and adjusted together
  se_df <- cma %>% 
    mutate(lower = full_value - q * se, upper = full_value + q * se, method = "CMA") %>%
    bind_rows(unadjusted %>% 
                mutate(lower = full_value - q * se, 
                       upper = full_value + q * se, method = "Unadjusted"))
  fosr_coef_ci_plot <- se_df %>%
    ggplot(aes(x = week, y = full_value)) +
    facet_wrap(~ arm) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 1) + 
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    theme_bw() +
    labs(x = "Week", y = "Coefficient") +
    guides(fill = guide_legend(title = "Method")) +
    scale_fill_manual(values = colors_light_dark) +
    theme(strip.background = element_rect(fill = "transparent", 
                                          color = "transparent")) +
    theme(legend.position = "bottom")
  
  # Construct contrasts b/w arms -------------------------------------------------
  # Get contrast for each BS iteration in long format
  bs_res_contrast <- arm_bs_res %>%
    group_by(week, b) %>%
    summarize(Competitive_Collaborative = value[arm == "Competitive"] - value[arm == "Collaborative"],
              Competitive_Supportive = value[arm == "Competitive"] - value[arm == "Supportive"],
              Collaborative_Supportive = value[arm == "Collaborative"] - value[arm == "Supportive"]) %>%
    pivot_longer(cols = -c(week, b), names_to = "pair", values_to = "contrast") 
  
  # Switch full data results into contrast format
  full_data_contrast <- full_arm_coef %>%
    group_by(week) %>%
    summarize(Competitive_Collaborative = value[arm == "Competitive"] - value[arm == "Collaborative"],
              Competitive_Supportive = value[arm == "Competitive"] - value[arm == "Supportive"],
              Collaborative_Supportive = value[arm == "Collaborative"] - value[arm == "Supportive"]) %>%
    pivot_longer(cols = -c(week), names_to = "pair", values_to = "contrast")
  
  # BS SE for each contrast
  se_res_contrast <- bs_res_contrast %>%
    group_by(pair, week) %>%
    summarize(se = sd(contrast)) %>%
    inner_join(full_data_contrast, by = c("pair", "week")) %>%
    rename(full_contrast = contrast) %>% # This is value computed with all data
    ungroup()
  
  # Unadjusted
  unadjusted_contrast <- se_res_contrast %>% mutate(q = 1.96)
  # CMA
  cma_contrast <- bs_res_contrast %>%
    left_join(se_res_contrast, by = c("pair", "week")) %>%
    group_by(pair, b) %>%
    summarize(db = max(abs(contrast - full_contrast) / se)) %>%
    group_by(pair) %>%
    summarize(q = quantile(db, .95)) %>% 
    right_join(se_res_contrast, by = "pair")
  
  # Combine CMA and unadjusted CIs
  se_df_contrast <- cma_contrast %>% 
    mutate(lower = full_contrast - q * se, upper = full_contrast + q * se, method = "CMA") %>%
    bind_rows(unadjusted_contrast %>% 
                mutate(lower = full_contrast - q * se, 
                       upper = full_contrast + q * se, method = "Unadjusted"))
  # Aesthetic purposes
  pair_labeller <- function(variable, value) {
    # Get the fake IDs for display
    gsub("_", " - ", value)
  } 
  # Plot contrasts and CI
  contrast_plot <- se_df_contrast %>%
    mutate(pair = forcats::fct_relevel(pair, c(#"Collaborative_Control", 
                                               #"Competitive_Control", "Supportive_Control", 
                                               "Competitive_Collaborative", "Competitive_Supportive", 
                                               "Collaborative_Supportive"))) %>%
    ggplot(aes(x = week, y = full_contrast)) +
    facet_wrap(~ pair, labeller = pair_labeller) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 1) + 
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#636363") + # TODO might be better behind
    theme_bw() +
    labs(x = "Week", y = "Coefficient") +
    scale_fill_manual(values = colors_light_dark) +
    theme(strip.background = element_rect(fill = "transparent", 
                                          color = "transparent")) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title = "Method"))
  
  # Coef p-values ------------------------------------------------------------
  
  # Find p-values
  coef_db <- arm_bs_res %>%
    left_join(se_res, by = c("arm", "week")) %>%
    group_by(arm, b) %>%
    summarize(db = max(abs(value - full_value) / se))
  
  alpha_options <- seq(from = 0, to = 1, length.out = 1000)
  # Takes time
  for (r in 1:nrow(se_res)) { # Iterate over arm-weeks
    arm <- unlist(se_res[r, "arm"])
    se <- se_res[r, "se"]
    full_value <- se_res[r, "full_value"]
    # Get list of all db for pair-week to look for quantiles
    db_list <- unlist(coef_db[coef_db$arm == arm, "db"])
    # Find all alphas such that 0 isn't contained
    valid_alphas <- c()
    for (alpha in alpha_options) {
      q <- quantile(db_list, 1 - alpha)
      lower <- full_value - q * se
      upper <- full_value + q * se
      if (lower > 0 | upper < 0) { # Doesn't contain 0
        valid_alphas <- c(valid_alphas, alpha)
      }
    }
    
    if (length(valid_alphas) == 0) {
      se_res[r, "p_value"] <- 1
    } else {
      se_res[r, "p_value"] <- min(valid_alphas)
    }
  }
  
  # Global p-values = min(pointwise)
  global_p_values <- se_res %>%
    group_by(arm) %>%
    summarize(global_p = min(p_value))
  
  pvalue_df <- se_res %>%
    ungroup()

  # Plot of p-values, just so we can see them
  pvalue_plot <- pvalue_df %>%
    ggplot(aes(x = week, y = p_value)) +
    geom_point() +
    facet_wrap(~ arm) +
    geom_hline(yintercept =  .05, linetype = "dashed", color = "#636363") +
    theme_bw() +
    labs(x = "Week", y = "p-value") +
    theme(strip.background = element_rect(fill = "transparent", 
                                          color = "transparent"))
  
  return(list(se_df, se_df_contrast, pvalue_df))
}


# Plots relating to intercept and baseline steps
bootstrap_plots_other <- function(period, bs_res, full_data_res, method, plot_suffix= "") {
  intercept_bs_res <- bs_res[[2]]
  baseline_coef_bs_res <- bs_res[[3]]
  start <- params[[period]][["start"]]
  end <- params[[period]][["end"]]
  
  full_intercept <- data.frame(value = full_data_res[[2]], week = start:end)
  full_baseline_coef <- data.frame(value = full_data_res[[3]], week = start:end)
  if (method == "fosr") {
    full_p_intercept <- full_data_res[[6]]
  }
  
  # Intercept --------------------------------------------------------------------
  se_res_intercept <- intercept_bs_res %>%
    group_by(week) %>%
    summarize(se = sd(value)) %>%
    inner_join(full_intercept, by = c("week")) %>%
    rename(full_value = value) # This is value computed with all data
  
  # Unadjusted Wald at each point
  unadjusted_intercept <- se_res_intercept %>% mutate(q = 1.96)
  # CMA adjusted intervals (I think)
  q_cma_intercept <- intercept_bs_res %>% # First get q
    left_join(se_res_intercept, by = c("week")) %>%
    group_by(b) %>%
    summarize(db = max(abs(value - full_value) / se)) %>%
    select(db) %>% unlist() %>%
    quantile(.95)
  cma_intercept <- se_res_intercept %>% mutate(q = q_cma_intercept)
  
  # Plot of unadjusted and adjusted together
  se_df_intercept <- cma_intercept %>% 
    mutate(lower = full_value - q * se, upper = full_value + q * se, method = "CMA") %>%
    bind_rows(unadjusted_intercept %>% 
                mutate(lower = full_value - q * se, 
                       upper = full_value + q * se, method = "Unadjusted"))
  intercept_plot <- se_df_intercept %>%
    ggplot(aes(x = week, y = full_value)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 1) + 
    geom_line() +
    theme_bw() +
    labs(x = "Week", y = "Coefficient", title = "Intercept") +
    guides(fill = guide_legend(title = "Method")) +
    scale_fill_manual(values = colors_light_dark) +
    theme(legend.position = "bottom")
  if (method == "fosr") {
    intercept_plot <- intercept_plot + 
      geom_hline(yintercept = full_p_intercept, linetype = "dashed", color = "#636363")
  }
  
  # Coefficient for baseline step count ------------------------------------------
  
  se_res_baseline <- baseline_coef_bs_res %>%
    group_by(week) %>%
    summarize(se = sd(value)) %>%
    inner_join(full_baseline_coef, by = c("week")) %>%
    rename(full_value = value) # This is value computed with all data
  
  # Unadjusted Wald at each point
  unadjusted_baseline <- se_res_baseline %>% mutate(q = 1.96)
  # CMA adjusted intervals (I think)
  q_cma_baseline <- baseline_coef_bs_res %>% # First get q
    left_join(se_res_baseline, by = c("week")) %>%
    group_by(b) %>%
    summarize(db = max(abs(value - full_value) / se)) %>%
    select(db) %>% unlist() %>%
    quantile(.95)
  cma_baseline <- se_res_baseline %>% mutate(q = q_cma_baseline)
  
  # Plot of unadjusted and adjusted together
  se_df_baseline <- cma_baseline %>% 
    mutate(lower = full_value - q * se, upper = full_value + q * se, method = "CMA") %>%
    bind_rows(unadjusted_baseline %>% 
                mutate(lower = full_value - q * se, 
                       upper = full_value + q * se, method = "Unadjusted"))
  baseline_plot <- se_df_baseline %>%
    ggplot(aes(x = week, y = full_value)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 1) + 
    geom_line() +
    theme_bw() +
    labs(x = "Week", y = "Coefficient", title = "Baseline steps beta(s)") +
    guides(fill = guide_legend(title = "Method")) +
    scale_fill_manual(values = colors_light_dark) +
    theme(legend.position = "bottom")
  
  return(list(se_df_intercept, se_df_baseline))
}

















































