# Predict post-intervention steps from intervention steps + arm + baseline_steps

all_na_int <- rowSums(!is.na(dat$steps_w[, 1:24])) <= 3
all_na_post <- rowSums(!is.na(dat$steps_w[, 25:36])) <= 3
all_na <- all_na_int | all_na_post

# FoFR main function -----------------------------------------------------------

fofr_function <- function(sample_ids) {
  # dat for sample (can contain repeat rows)
  imp_dat <- data.frame(id = sample_ids) %>%
    left_join(dat, by = "id")
  
  # Impute data with FPCA -------------------------------------------------------
  pop_avg <- colMeans(imp_dat$steps_w, na.rm = T)
  
  # Intervention period
  start <- 1
  end <- 24
  steps_data <- imp_dat$steps_w[, start:end]
  fpca_res <- fpca.face(unclass(steps_data), knots = 15, pve = .95)
  int_imp <- fpca_res$scores %*% t(fpca_res$efunctions) + pop_avg[start:end]
  # If missing, use FPCA imputation
  imp_dat$steps_w[, start:end] <- ifelse(is.na(imp_dat$steps_w[, start:end]),
                                         int_imp, imp_dat$steps_w[, start:end])
  
  # Follow-up period - not used for imputation
  start <- 25
  end <- 36
  steps_data <- imp_dat$steps_w[, start:end]
  fpca_res_post_bs <- fpca.face(unclass(steps_data), knots = 8, pve = .95)
  post_imp <- fpca_res_post_bs$scores %*% t(fpca_res_post_bs$efunctions) + pop_avg[start:end]
  imp_dat$steps_w[, start:end] <- ifelse(is.na(imp_dat$steps_w[, start:end]),
                                         post_imp, imp_dat$steps_w[, start:end])
  
  # Fix after imputation
  imp_dat$steps_w_int <- I(imp_dat$steps_w[, 1:24])
  imp_dat$steps_w_post <- I(imp_dat$steps_w[, 25:36])
  
  # Prepare data matrices for mgcv directly --------------------------------------
  
  # Construct matrices and df
  n <- nrow(imp_dat)
  S <- 25:36 # Outcome domain
  U <- 1:24 # Predictor domain
  W <- imp_dat$steps_w_post # Outcome
  X <- imp_dat$steps_w_int # Predictor
  
  # Predictor side
  Xv <- kronecker(X, matrix(1, nrow = length(S), ncol = 1)) # Predictor matrix
  q <- (U[length(U)] - U[1]) / length(U) / 3 * 
    c(1, rep(c(4, 2), length = length(U) - 2), 1) # Quadrature weights
  L <- kronecker(matrix(1, nrow = n * length(S), ncol = 1), t(q)) # n |S| x |U|
  XvL <- Xv * L # "Predictors", X * q, n |S| x |U|
  
  # Outcome side
  wv <- as.vector(t(W)) # n * |S|
  
  # Domains
  Uv <- kronecker(matrix(1, nrow = n * length(S), ncol = 1), t(U)) # n |S| x |U|
  sv <- kronecker(matrix(1, nrow = n, ncol = 1), S) # n |S| x 1
  Sv <- kronecker(sv, matrix(1, nrow = 1, ncol = length(U))) # n |S| x |U|
  
  # Combine into df
  fofr_df <- data.frame(wv, I(XvL), sv, I(Sv), I(Uv))
  
  # Add other covariates
  # Refit model with arm and baseline steps
  #fofr_df$baseline_steps <- rep(imp_dat$baseline_steps, each = 12)
  #fofr_mgcv_df$arm_long <- rep(imp_dat$arm, each = 12)
  fofr_df$start_day <- rep(imp_dat$start_day, each = 12)
  fofr_df$baseline_steps <- rep(imp_dat$ESE_sticking_to_it, each = 12)
  fofr_df$age <- rep(imp_dat$age, each = 12)
  fofr_df$gender <- rep(as.numeric(imp_dat$gender == "Male"), each = 12)
  
  # Helper function that returns mask of dim XvL
  # Probably unnecessary
  match_XvL <- function(long_dat) {
    long_dat %>%
      rep(each = 12) %>% # 12 = |S|
      matrix(ncol = 1) %>%
      kronecker(matrix(1, nrow = 1, ncol = 24)) # 24 = |U|
  }
  
  fofr_df$collaborative <- XvL * match_XvL(imp_dat$arm == "Collaborative")
  fofr_df$competitive <- XvL * match_XvL(imp_dat$arm == "Competitive")
  fofr_df$supportive <- XvL * match_XvL(imp_dat$arm == "Supportive")
  fofr_df$control <- XvL * match_XvL(imp_dat$arm == "Control")
  #fofr_mgcv_df$baseline_steps <- XvL * match_XvL(imp_dat$baseline_steps)
  
  # Necessary for residual modeling
  # Can't be the original id because then R.E. are shared
  fofr_df$re_id <- factor(rep(1:length(sample_ids), each = 12))
  
  # Add post-intervention period eigenfunctions
  ef_df <- data.frame(fpca_res_post_bs$efunctions)
  colnames(ef_df) <- paste0("ef_", 1:ncol(ef_df))
  ef_df$week <- 25:36
  fofr_df <- fofr_df %>%
    left_join(ef_df, by = join_by(Sv == week))
  
  # TODO can change
  #sv <- as.vector(sv)
  
  # Fit model ----------------------------------------------------------------
  fofr_sample <- bam(wv ~ s(sv, bs = "ps", k = 8) +
                       # te(Uv, Sv, by = XvL, bs = "ps", m = list(c(2, 1), c(2, 1)), #sp = c(1.3, 3.9),
                       #    k = c(5, 5)) +
                       te(Uv, Sv, by = collaborative, bs = "ps", #m = list(c(2, 1), c(2, 1)), #sp = c(1.3, 3.9),
                          k = c(5, 5)) +
                       te(Uv, Sv, by = competitive, bs = "ps", #m = list(c(2, 1), c(2, 1)), #sp = c(1.2, 3.9),
                          k = c(5, 5)) +
                       te(Uv, Sv, by = supportive, bs = "ps", #m = list(c(2, 1), c(2, 1)), #sp = c(1.1, 3.9),
                          k = c(5, 5)) +
                       te(Uv, Sv, by = control, bs = "ps", #m = list(c(2, 1), c(2, 1)), #sp = c(1.2, 3.9),
                          k = c(5, 5)) +
                       #s(sv, by = baseline_steps, bs = "ps", k = 10) +
                       s(re_id, by = ef_1, bs = "re") + s(re_id, by = ef_2, bs = "re") +
                       start_day + age + gender + baseline_steps,
                     method = "fREML", #chunk.size = 10000, 
                     data = fofr_df)
  summary_fofr <- summary(fofr_sample)
  
  # Extract coefficients-------------------------------------------------------
  uv_pred <- seq(1, 24, length = length(U))
  sv_pred <- seq(25, 36, length = length(S))
  df_pred <- expand.grid(Uv = uv_pred, Sv = sv_pred) %>%
    mutate(sv = rep(S, length(U)), XvL = 1, control = 1, collaborative = 1, supportive = 1, 
           competitive = 1, steps_long = 1, ef_1 = 1, ef_2 = 1, 
           re_id = fofr_df$re_id[1], start_day = 1, gender = 1,  age = 1, baseline_steps = 1)
  pred_res <- predict(fofr_sample, type = "terms", newdata = df_pred, se = F)
  
  # Extract intercept
  fixed_int <- summary_fofr$p.coeff[[1]]
  intercept <- data.frame(s = 25:36, 
                          value = fixed_int + pred_res[1:12, "s(sv)"])
  # Extract arm
  arm_coef <- data.frame(collaborative = pred_res[, "te(Uv,Sv):collaborative"],
                         competitive = pred_res[, "te(Uv,Sv):competitive"],
                         supportive = pred_res[, "te(Uv,Sv):supportive"],
                         control = pred_res[, "te(Uv,Sv):control"])
  arm_coef <- cbind(arm_coef, expand.grid(u = uv_pred, s = sv_pred))
  
  list(intercept = intercept, arm_coef = arm_coef, model = fofr_sample)
}

# Bootstrap --------------------------------------------------------------------

B <- 300
bs_arm_coef <- data.frame()
bs_intercept <- data.frame()
bs_baseline_coef <- data.frame()

set.seed(103856)
b <- 1
while (b <= B) {
  print(b)
  tryCatch({
    # Sample ids instead of reconstructing dfs
    bs_ids <- sample(dat$id[!all_na], size = sum(!all_na), replace = T)
    # Run model and coefficients
    fofr_res <- fofr_function(bs_ids)
    # Append results
    bs_arm_coef <- rbind(bs_arm_coef, fofr_res$arm_coef %>% mutate(b = b))
    bs_intercept <- rbind(bs_intercept, fofr_res$intercept %>% mutate(b = b))
    bs_baseline_coef <- rbind(bs_baseline_coef, fofr_res$baseline_coef %>% mutate(b = b))
    b <- b + 1
  }, error = function(e) {
    print(e)
  })
}

# Full data results & plots ----------------------------------------------------

full_fofr_res <- fofr_function(dat$id[!all_na])
full_intercept <- full_fofr_res$intercept
full_arm_coef <- full_fofr_res$arm_coef
full_arm_coef_long <- full_arm_coef %>%
  pivot_longer(cols = c("collaborative", "competitive", "supportive", "control"),
               names_to = "arm", values_to = "value") %>%
  mutate(arm = fct_recode(arm, 
                          "Collaborative" = "collaborative",
                          "Competitive" = "competitive", "Supportive" = "supportive", 
                          "Control" = "control")) %>%
  mutate(arm = fct_relevel(arm,
                           "Collaborative", "Competitive", "Supportive", "Control"))

# Intercept
full_intercept %>%
  ggplot(aes(x = s, y = value)) +
  geom_line() +
  theme_bw() +
  labs(x = "Week", y = "Intercept")

# Define a custom scale so it stays consistent
custom_gradient <- scale_fill_gradientn(limits = c(-.5, .5),
                                        colors = c("#3E502B","#8cae68","#c6d7b4","#ffffff","#b2d0d3","#64a0a6", "#34575B"),
                                        values = scales::rescale(c(-.5, -.18, -.05, 0, .05, .18, .5)))

# Smaller range (good for coefficient without CI)
custom_gradient_lim <- scale_fill_gradientn(limits = c(-.31, .31),
                                            colors = c("#3E502B","#8cae68","#c6d7b4","#ffffff","#b2d0d3","#64a0a6", "#34575B"),
                                            values = scales::rescale(c(-.31, -.18, -.08, 0, .08, .18, .31)))

# Bootstrap CI and plots -------------------------------------------------------

# Put arm dfs in long form
bs_arm_coef_long <- bs_arm_coef %>%
  pivot_longer(cols = c("collaborative", "competitive", "supportive", "control"),
               names_to = "arm", values_to = "value") %>%
  mutate(arm = fct_recode(arm, 
                          "Collaborative" = "collaborative",
                          "Competitive" = "competitive", "Supportive" = "supportive", 
                          "Control" = "control")) %>%
  mutate(arm = fct_relevel(arm,
                           "Collaborative", "Competitive", "Supportive", "Control"))

se_res <- bs_arm_coef_long %>%
  group_by(arm, s, u) %>%
  summarize(se = sd(value)) %>%
  inner_join(full_arm_coef_long, by = c("arm", "s", "u")) %>%
  rename(full_value = value) # This is value computed with all data

# Unadjusted Wald at each point
unadjusted <- se_res %>% mutate(q = 1.96)
# CMA adjusted intervals
cma <- bs_arm_coef_long %>% # First get q
  left_join(se_res, by = c("arm", "s", "u")) %>%
  group_by(arm, b) %>%
  # Max over s and u instead of just s
  summarize(db = max(abs(value - full_value) / se)) %>%
  group_by(arm) %>%
  summarize(q = quantile(db, .95)) %>% 
  right_join(se_res, by = "arm")

se_df <- cma %>% 
  mutate(method = "CMA") %>%
  bind_rows(unadjusted %>% mutate(method = "Unadjusted")) %>%
  mutate(lower = full_value - q * se, upper = full_value + q * se)

## Plot lower and upper bounds of CI surface-----------------------------------

combined_fofr_ci <- se_df %>%
  select(-c(q, se)) %>%
  pivot_longer(cols = c(full_value, lower, upper), names_to = "type", values_to = "value") %>%
  mutate(type = recode(type, "full_value" = "Point estimate", 
                       "lower" = "Lower", "upper" = "Upper")) %>%
  ggplot(aes(x = u, y = s, fill = value)) +
  facet_grid(arm ~ type) +
  geom_tile() +
  coord_equal(expand = F) + # Scale x and y same, prevent border
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Week (intervention period)", y = "Week (follow-up period)",
       title = "FoFR coefficients", fill = "Value") +
  theme(strip.background = element_rect(fill = "transparent", color = "transparent")) +
  custom_gradient

# Mask out cells that aren't significant
fofr_p_plot <- se_df %>%
  filter(method == "CMA") %>%
  mutate(full_value = ifelse((lower < 0) & (upper > 0), 0, full_value)) %>%
  ggplot(aes(x = u, y = s, fill = full_value)) +
  facet_wrap(~ arm, nrow = 1) +
  geom_tile() +
  coord_equal(expand = F) + # Scale x and y same, prevent border
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Week (intervention period)", y = "Week (follow-up period)",
       fill = "Coefficient") +
  theme(strip.background = element_rect(fill = "transparent", color = "transparent")) +
  custom_gradient_lim

# Mask out plot + full coef plot
combined_p_full_fofr_plot <- se_df %>%
  filter(method == "CMA") %>%
  mutate(full_value = ifelse((lower < 0) & (upper > 0), 0, full_value)) %>% 
  mutate(value = full_value) %>%
  select(arm, s, u, value) %>%
  mutate(row = "Significant") %>%
  bind_rows(full_arm_coef_long %>% mutate(row = "All")) %>%
  ggplot(aes(x = u, y = s, fill = value)) +
  facet_grid(row ~ arm) +
  geom_tile() +
  coord_equal(expand = F) + # Scale x and y same, prevent border
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Week (intervention period)", y = "Week (follow-up period)",
       fill = "Coefficient") +
  theme(strip.background = element_rect(fill = "transparent", color = "transparent")) +
  theme(legend.key.width = unit(2, "cm")) +
  custom_gradient_lim


# Cross-section plot -----------------------------------------------------------


target_s <- 26

cross_sec_plot <- se_df %>%
  filter(s == target_s) %>%
  ggplot(aes(x = u, y = full_value)) +
  facet_wrap(~ arm, nrow = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 1) + 
  geom_line() +
  theme_bw() +
  labs(x = "Week", y = "Coefficient") +
  guides(fill = guide_legend(title = "Method")) +
  scale_fill_manual(values = colors_light_dark) +
  #theme(legend.position = "bottom") +
  theme(panel.spacing.x = unit(.75, "lines")) +
  theme(strip.background = element_rect(fill = "transparent", color = "transparent"))

# Plot BS replications individually
bs_cross_sec_plot <- bs_arm_coef_long %>%
  mutate(type = "Bootstrap") %>%
  bind_rows(full_arm_coef_long %>% mutate(type = "Full data", b = 0)) %>%
  filter(abs(s - target_s) < .1) %>%
  filter(b %in% 0:20) %>%
  mutate(b = factor(b, levels = c(1:20, 0))) %>%
  ggplot(aes(x = u, y = value, group = b, color = type, size = type)) +
  scale_size_manual(values = c(.7, .7)) +
  scale_color_manual(values = c("#B2B2B2", "black")) +
  geom_line(alpha = .8) +
  facet_wrap(~ arm, nrow = 1) +
  labs(x = "Week", y = "Coefficient", color = "Data sample", size = "Data sample") +
  theme_bw() +
  #guides(size = element_blank()) +
  theme(#legend.position = "bottom",
    panel.spacing.x = unit(.75, "lines")) +
  theme(strip.background = element_rect(fill = "transparent", color = "transparent"))

# Patchwork together the two
combined_bs_cross_sec_plot <- (cross_sec_plot +
                                 theme(axis.title.x = element_blank(),
                                       axis.text.x = element_blank(),
                                       axis.ticks.x = element_blank())) /
  (bs_cross_sec_plot + theme(strip.text.x = element_blank())) 
