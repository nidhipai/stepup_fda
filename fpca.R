# FPCA

# Intervention period ----------------------------------------------------------
start <- 1
end <- 24
steps_data <- dat$steps_w[, start:end]

all_na_int <- rowSums(!is.na(steps_data)) <= 3
steps_f <- steps_data[!all_na_int, ] # Remove rows with few values
fpca_res_int <- fpca.face(unclass(steps_f), knots = 15, pve = .95)

# Post-intervention period -----------------------------------------------------
start <- 25
end <- 36
steps_data <- dat$steps_w[, start:end]

all_na_post <- rowSums(!is.na(steps_data)) <= 3
steps_f <- steps_data[!all_na_post, ] # Remove rows with few values
fpca_res_post <- fpca.face(unclass(steps_f), knots = 8, pve = .95)

# Regression of scores on covariates ---------------------------------------------------------

fpca_coefs <- function(fpca_res, all_na) {
  num_ef <- fpca_res$npc
  scores <- fpca_res$scores
  fpca_dat <- dat[!all_na, ]
  
  covariates <- c("age", "gender", "baseline_steps", "wearabledev", "psqi_quality", 
                  "bf_extroversion", "bf_Agreeableness", "bf_Conscientiousness",
                  "bf_Neuroticism", "bf_openness", "grit_score", "ESE_sticking_to_it",
                  "dospert_health_safety", "dospert_social", "mosss_overall", "start_day")
  covariates_formula <- paste0(covariates, collapse = "+")
  
  # Compare models with arm * covariates, with arm, and without arm
  for (ef in 1:num_ef) {
    print("-----------------------------------------------")
    print(ef)
    fpca_dat$scores <- scores[, ef]
    no_arm_model <- lm(as.formula(paste0("scores ~ ", covariates_formula)), 
                       data = fpca_dat) 
    arm_model <- lm(as.formula(paste0("scores ~ ", covariates_formula, " + arm")),
                    data = fpca_dat) 
    full_model <- lm(as.formula(paste0("scores ~ arm * (", covariates_formula, ")")),
                    data = fpca_dat)
    print(anova(no_arm_model, arm_model, full_model))
  }
  
  # Regression of scores on covariates + arm
  coefs_all <- data.frame(matrix(nrow = 0, ncol = 8)) 
  for (ef in 1:num_ef) {
    fpca_dat$scores <- scores[, ef]
    lm_res <- lm(as.formula(paste0("scores ~ ", covariates_formula, " + arm")),
                 data = fpca_dat)
    lm_coefs <- data.frame(summary(lm_res)$coefficients)
    colnames(lm_coefs) <- c("estimate", "std_error", "t_value", "p_value")
    lm_coefs$term <- rownames(lm_coefs)
    rownames(lm_coefs) <- 1:nrow(lm_coefs)
    lm_coefs$eigenfunction <- paste0("Eigenfunction ", ef)
    lm_coefs$p_value_adj <- p.adjust(lm_coefs$p_value, method = "holm")
    coefs_all <- rbind(coefs_all, lm_coefs)
  }
  coefs_all <- coefs_all[, c("eigenfunction", "term", "estimate", "std_error", 
                             "t_value", "p_value", "p_value_adj")] # Reorder cols
  coefs_all <- coefs_all[coefs_all$term != "(Intercept)", ]
  coefs_all
}

# Intervention period
int_coefs <- fpca_coefs(fpca_res_int, all_na_int)
View(int_coefs)

# Follow-up period
post_coefs <- fpca_coefs(fpca_res_post, all_na_post)
View(post_coefs)

# Test the contrast of each arm in intervention period
fpca_dat <- dat[!all_na_int, ]
fpca_dat$scores <- fpca_res_int$scores[, 1]
lm_res <- lm(as.formula(paste0("scores ~ ", covariates_formula, " + arm")),
             data = fpca_dat)
emmeans_res <- emmeans::emmeans(lm_res, c("arm"))
contrasts_list <- list(col_comp = c(0, 1, -1, 0), col_supp = c(0, 1, 0, -1), comp_supp = c(0, 0, 1, -1))
emmeans::contrast(emmeans_res, contrasts_list, adjust = "scheffe")




















