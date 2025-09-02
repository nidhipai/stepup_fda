library(geepack)

gee_df <- dat %>%
  select(id, arm, age, gender, baseline_steps, start_day, start_date, cohort_id) %>%
  mutate(cohort_id = as.numeric(cohort_id) + 700,
         cohort_id_2 = ifelse(is.na(cohort_id), id, cohort_id)) %>%
  cbind(unAsIs(dat$steps_w_int)) %>%
  pivot_longer(cols = -c(id, arm, age, gender, baseline_steps, start_day, 
                         start_date, cohort_id, cohort_id_2), names_to = "week", values_to = "log_steps") %>%
  mutate(week = as.numeric(week) - 1,
         date = week * 7 + start_date,
         month = factor(month(date)))


# Fit a GEE model with exchangeable covariance and the same covariates & outcome 
# as the FoSR model
# Use factor effects coding since our model includes interaction terms,
# ensuring that tests correspond to main effects.

gee_model <- geeglm(log_steps ~ week * arm + age + gender + start_day +
                      baseline_steps,
                    corstr = "exchangeable", id = id, data = gee_df)
summary(gee_model)$coefficients

# Plot coefficients over week and SE bars --------------------------------------

# Calculate coef by week and arm + variance
M <- data.frame(Collaborative = c(rep(0, 24), rep(1, 24), rep(0, 48)),
                Competitive = c(rep(0, 48), rep(1, 24), rep(0, 24)),
                Supportive = c(rep(0, 72), rep(1, 24)),
                week = rep(0:23, 4)) %>%
  mutate(Collaborative_w = Collaborative * week,
         Competitive_w = Competitive * week,
         Supportive_w = Supportive * week) %>% as.matrix()
coef_names_vec <- c("armCollaborative", "armCompetitive", "armSupportive", 
                    "week", "week:armCollaborative", "week:armCompetitive",
                    "week:armSupportive")
coef <- gee_model$coef[coef_names_vec]
pred <- M %*% coef
vcov_gee <- vcov(gee_model)[coef_names_vec, coef_names_vec]
covar_pred <- diag(M %*% vcov_gee %*% t(M))

pred_df <- data.frame(pred, var = covar_pred) %>%
  mutate(arm = rep(c("Control", "Collaborative", "Competitive", "Supportive"),
                   each = 24), week = rep(1:24, 4)) %>%
  mutate(lower = pred - 1.96 * sqrt(var),
         upper = pred + 1.96 * sqrt(var)) %>%
  select(arm, week, pred, lower, upper)
pred_df %>%
  ggplot(aes(x = week, y = pred, color = arm)) +
  geom_hline(yintercept = 0, color = "#636363", linetype = "longdash") +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "Week", y = "Coefficient (arm + week + interaction)") +
  scale_color_manual(values = colors, guide = guide_legend(title = "Arm")) 
# With standard error (unadjusted), faceted
gee_coef_plot <- pred_df %>%
  ggplot(aes(x = week, y = pred)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1, 
              fill = colors_light_dark[2]) +
  #geom_point() +
  geom_hline(yintercept = 0, color = "#636363", linetype = "longdash") +
  geom_line() +
  facet_wrap(~arm) +
  labs(x = "Week (intervention period)", y = "Arm effect (arm + week + interaction)") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent", color = "transparent")) 