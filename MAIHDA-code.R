
# ==============================================================================
# Title: Multilevel Analysis of Individual Heterogeneity and Discriminatory 
#        Accuracy (MAIHDA) 
# Description: Bayesian ordinal multilevel modeling, LOO cross-validation, 
#              interaction effects, and trajectory analysis for perceived health.
# Author: Belachew et al.
# Date: 2026-02-06
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
# Install missing packages automatically

libs <- c("haven", "tidyverse", "brms", "tidybayes", "bayesplot", "patchwork", 
          "openxlsx", "matrixStats", "ggrepel", "scales", "RColorBrewer")
new_libs <- libs[!(libs %in% installed.packages()[,"Package"])]
if(length(new_libs)) install.packages(new_libs)

invisible(lapply(libs, library, character.only = TRUE))

# Create structured output directories 
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/diagnostics", showWarnings = FALSE)

# Load data 
cat("\nLoading data...\n")
df <- readRDS("MAIHDA_data_final.rds") 

# 2. DATA PREPARATION ----------------------------------------------------------
cat("Preparing variables and intersectional strata...\n")
df <- df %>%
  mutate(
    Study_age = factor(Study_age, levels = c("perceived_health16", "perceived_health33"), labels = c(1, 2)),
    across(c(gender, family_type, parental_educ, income_tertiles, Close_friends), as.numeric),
    # Constructing unique intersectional stratum IDs
    stratum = gender*10000 + family_type*1000 + parental_educ*100 + income_tertiles*10 + Close_friends,
    across(c(gender, family_type, parental_educ, income_tertiles, Close_friends, stratum), as.factor)
  )

# 3. MODEL IMPLEMENTATION & CONVERGENCE DIAGNOSTICS ----------------------------
# Setting up shared brms parameters, including saving all parameters for robust LOO
brms_settings <- list(
  family = cumulative("logit"), warmup = 1000, iter = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15), seed = 2024,
  save_pars = save_pars(all = TRUE) # Crucial for moment_match LOO comparisons
)

cat("\nFitting Model 1: Null Model...\n")
null_model <- do.call(brm, c(list(
  formula = Perceived_health ~ Study_age + (1|stratum/project_ID), data = df), 
  brms_settings))

trace_null <- mcmc_plot(null_model, type = "trace", pars = c("^b_", "^sd_")) +
  theme_minimal() + labs(title = "Trace Plots: Null Model Convergence")

cat("\nFitting Model 2: Main Intersectional Model...\n")
main_formula <- Perceived_health ~ (gender + family_type + parental_educ + 
                                  income_tertiles + Close_friends) * Study_age + (1|stratum/project_ID)

main_model <- do.call(brm, c(list(formula = main_formula, data = df), brms_settings))

trace_main <- mcmc_plot(main_model, type = "trace", pars = c("^b_", "^sd_")) +
  theme_minimal() + labs(title = "Trace Plots: Main Model Convergence")

# 4. MODEL COMPARISON (LOO-CV) -------------------------------------------------
cat("\nConducting Model Comparison (LOO-CV)...\n")
loo_null <- loo(null_model, moment_match = TRUE)
loo_main <- loo(main_model, moment_match = TRUE)

model_comparison <- loo_compare(loo_null, loo_main)
print(model_comparison)

# Save LOO results to a text file for record
sink("outputs/diagnostics/Model_Comparison_Results.txt")
print("LOO Cross-Validation Comparison:")
print(model_comparison)
sink()

# 5. POSTERIOR PREDICTIONS & TIDYING -------------------------------------------
cat("\nGenerating Posterior Predictions...\n")
posterior_probs <- posterior_epred(main_model, newdata = df, re_formula = NULL, summary = FALSE)
good_very_probs <- posterior_probs[, , 2] + posterior_probs[, , 3]

df_preds <- df %>%
  ungroup() %>%
  mutate(
    pred_well = colMedians(good_very_probs),
    pred_well_lower = colQuantiles(good_very_probs, probs = 0.025),
    pred_well_upper = colQuantiles(good_very_probs, probs = 0.975),
    Age = factor(Study_age, levels = c(1, 2), labels = c("Adolescence (16 yrs)", "Early adulthood (33 yrs)")),
    Sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    Income = factor(income_tertiles, levels = c(1, 2, 3), labels = c("Low", "Middle", "High")),
    Education = factor(parental_educ, levels = c(1, 2, 3), labels = c("Basic", "Secondary", "Tertiary")),
    Friends = factor(Close_friends, levels = c(1, 2, 3), labels = c("No", "1-2", ">2")),
    Strata = factor(stratum)
  )

# 6. INTERACTION PLOT (Fig 1) --------------------------------------------------
cat("\nPlotting Figure 1 (Interaction Effects)...\n")
conditions <- data.frame(
  income_tertiles = unique(df$income_tertiles),
  gender = factor("1", levels = levels(df$gender)),  
  family_type = factor("1", levels = levels(df$family_type)),  
  parental_educ = factor("1", levels = levels(df$parental_educ)),  
  Close_friends = factor("1", levels = levels(df$Close_friends))  
)

health_effects <- conditional_effects(
  main_model, effects = "Study_age", conditions = conditions,
  categorical = TRUE, prob = 0.95
)

ef_df <- health_effects$Study_age
reshaped_df <- data.frame()
unique_combos <- unique(ef_df[, c("Study_age", "income_tertiles")])

for (i in 1:nrow(unique_combos)) {
  age <- unique_combos$Study_age[i]
  income <- unique_combos$income_tertiles[i]
  subset_df <- ef_df[ef_df$Study_age == age & ef_df$income_tertiles == income, ]
  
  new_row <- data.frame(
    Study_age = age, income_tertiles = income,
    prob_poor_moderate = subset_df$estimate__[subset_df$cats__ == "Poor/Moderate"],
    prob_good = subset_df$estimate__[subset_df$cats__ == "Good"],
    prob_very_good = subset_df$estimate__[subset_df$cats__ == "Very good"],
    lower_poor_moderate = subset_df$lower__[subset_df$cats__ == "Poor/Moderate"],
    lower_good = subset_df$lower__[subset_df$cats__ == "Good"],
    lower_very_good = subset_df$lower__[subset_df$cats__ == "Very good"],
    upper_poor_moderate = subset_df$upper__[subset_df$cats__ == "Poor/Moderate"],
    upper_good = subset_df$upper__[subset_df$cats__ == "Good"],
    upper_very_good = subset_df$upper__[subset_df$cats__ == "Very good"]
  )
  
  new_row$prob_2_or_above <- new_row$prob_good + new_row$prob_very_good
  new_row$lower_2_or_above <- 1 - new_row$upper_poor_moderate
  new_row$upper_2_or_above <- 1 - new_row$lower_poor_moderate
  reshaped_df <- rbind(reshaped_df, new_row)
}

reshaped_df$x_pos <- as.numeric(as.character(reshaped_df$Study_age))
reshaped_df$x_pos <- reshaped_df$x_pos + as.numeric(reshaped_df$income_tertiles) * 0.15 - 0.3

p1 <- ggplot(reshaped_df, aes(x = x_pos, y = prob_2_or_above, color = income_tertiles, shape = income_tertiles)) +
  geom_point(size = 3) + geom_errorbar(aes(ymin = lower_2_or_above, ymax = upper_2_or_above), width = 0.05) +
  labs(y = "Probability of good or very good perceived health", x = "Study age") +
  scale_x_continuous(breaks = c(1, 2), labels = c("Adolescence (16 yrs)", "Early adulthood (33 yrs)")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  scale_color_brewer(palette = "Set2", name = "Income Level", labels = c("1" = "Low", "2" = "Medium", "3" = "High")) +
  scale_shape_manual(name = "Income Level", values = c(16, 17, 15), labels = c("1" = "Low", "2" = "Medium", "3" = "High")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major = element_line(color = "grey70", linewidth = 0.3), panel.grid.minor = element_line(color = "grey85", linewidth = 0.1),
        panel.border = element_rect(color = "grey50", fill = NA), plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"), plot.margin = margin(t = 10, r = 10, b = 10, l = 45))

# 7. TREND FIGURES (Fig 2 & 3) ------------------------------------------
cat("Plotting Figures 2 & 3 (Trends over study age)...\n")
plot_data_long <- df_preds %>% group_by(Age, Sex, Income, Strata) %>% summarize(mean_prob = mean(pred_well), .groups = "drop")

p2 <- ggplot(plot_data_long, aes(x = Age, y = mean_prob, group = Strata, color = Income)) +
  geom_line(linewidth = 0.8, alpha = 0.6) + scale_color_manual(values = c("Low"="#FF9999", "Middle"="#66B2FF", "High"="#66CC99")) +
  labs(x = "Study age", y = "Prob. of Good/Very Good Health") + theme_minimal(base_size = 10) +
  theme(axis.line = element_line(color = "gray70"), panel.border = element_rect(color = "grey85", fill = NA))

p3 <- ggplot(plot_data_long, aes(x = Age, y = mean_prob, group = Strata, color = Sex)) +
  geom_line(linewidth = 0.8, alpha = 0.6) + scale_color_manual(values = c("Male"="#66B2FF", "Female"="#FF9999")) +
  labs(x = "Study age", y = "Prob. of Good/Very Good Health", color = "Sex") + theme_minimal(base_size = 10) +
  theme(axis.line = element_line(color = "gray70"), panel.border = element_rect(color = "grey50", fill = NA))

# 8. BOXPLOT FIGURES (Fig S1 & S2) -----------------------------------------------
cat("Plotting Figures 4 & 5 (Boxplots)...\n")
educ_colors <- c("Basic"="#FF9999", "Secondary"="#66B2FF", "Tertiary"="#66CC99")
p4 <- ggplot(df_preds, aes(x = Education, y = pred_well, fill = Education)) +
  facet_wrap(~ Age) + geom_boxplot(alpha = 0.7, width = 0.6, size = 0.3, outlier.shape = NA) +
  geom_jitter(color = "darkgray", width = 0.2, alpha = 0.3, size = 0.8) + scale_fill_manual(values = educ_colors) +
  labs(x = "Parental Education Level", y = "Predicted Probability") + theme_minimal() + 
  theme(legend.position = "none", axis.line = element_line(color = "gray70"), panel.border = element_rect(color = "gray70", fill = NA))

friends_colors <- c("No" = "#FF9999", "1-2" = "#66B2FF", ">2" = "#66CC99")
p5 <- ggplot(df_preds, aes(x = Friends, y = pred_well, fill = Friends)) +
  facet_wrap(~ Age) + geom_boxplot(alpha = 0.7, width = 0.6, size = 0.3, outlier.shape = NA) +
  geom_jitter(color = "darkgray", width = 0.2, alpha = 0.3, size = 0.8) + scale_fill_manual(values = friends_colors) +
  labs(x = "Number of Close Friends", y = "Predicted Probability") + theme_minimal() +
  theme(legend.position = "none", axis.line = element_line(color = "gray70"), panel.border = element_rect(color = "gray70", fill = NA))

# 9. STRATUM RESIDUALS (CATERPILLAR PLOTS) (Figure S3) -------------------------------------
cat("Plotting Residuals...\n")
get_residuals <- function(model) {
  ranef(model)$stratum %>% 
    as_tibble(rownames = "stratum") %>%
    rename(estimate = Estimate.Intercept, lower = Q2.5.Intercept, upper = Q97.5.Intercept) %>%
    mutate(rank = rank(estimate))
}

re_null <- get_residuals(null_model)
re_main <- get_residuals(main_model)
re_main$rank <- re_null$rank[match(re_main$stratum, re_null$stratum)]

create_clean_re_plot <- function(data, y_label) {
  y_lims <- c(min(re_null$lower, re_main$lower) - 0.2, max(re_null$upper, re_main$upper) + 0.2)
  ggplot(data, aes(x = rank, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, color = "gray30", alpha = 0.4) +
    geom_point(color = "black", size = 0.7) + coord_cartesian(ylim = y_lims) + theme_bw() +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          axis.line = element_line(color = "black", linewidth = 0.4), 
          axis.title = element_text(face = "plain"), aspect.ratio = 1.2) + labs(x = "Stratum Rank", y = y_label)
}

combined_residuals <- (create_clean_re_plot(re_null, "Residuals (Null Model)") + 
                       create_clean_re_plot(re_main, "Residuals (Main Model)")) + 
                       plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "plain", size = 14))

# 10. TRAJECTORY ANALYSIS & PLOTTING (Figure 4) -------------------------------------------
cat("\nRunning Trajectory Analysis...\n")
stratum_sizes <- df %>% group_by(stratum) %>% summarize(size = n(), .groups = "drop")
large_strata <- stratum_sizes %>% filter(size >= 5) %>% pull(stratum)

stratum_predictions <- df_preds %>%
  filter(stratum %in% large_strata) %>% group_by(stratum, Study_age) %>%
  summarize(pred_mean = mean(pred_well, na.rm = TRUE), pred_lower = mean(pred_well_lower, na.rm = TRUE),
            pred_upper = mean(pred_well_upper, na.rm = TRUE), n_obs = n(), .groups = "drop") %>%
  left_join(df %>% dplyr::select(stratum, gender, family_type, parental_educ, income_tertiles, Close_friends) %>% distinct(), by = "stratum") %>%
  mutate(Gender = factor(gender, levels = c(1, 2), labels = c("Boys", "Girls")),
         Study_age_label = factor(Study_age, levels = c(1, 2), labels = c("Adolescence (16 yrs)", "Early adulthood (33 yrs)")),
         Income = factor(income_tertiles, levels = c(1, 2, 3), labels = c("Low", "Middle", "High")))

trajectory_data <- stratum_predictions %>%
  pivot_wider(id_cols = c(stratum, gender, Gender, family_type, parental_educ, Income, Close_friends, income_tertiles),
              names_from = Study_age, values_from = c(pred_mean, pred_lower, pred_upper)) %>%
  mutate(difference = pred_mean_2 - pred_mean_1)

increasing_strata <- trajectory_data %>% filter(difference > 0) %>% arrange(desc(difference)) %>% head(5) 
declining_strata <- trajectory_data %>% filter(difference < 0) %>% arrange(difference) %>% head(5) 

increasing_plot_data <- stratum_predictions %>% filter(stratum %in% increasing_strata$stratum) %>%
  left_join(increasing_strata %>% dplyr::select(stratum, difference), by = "stratum") %>% mutate(trajectory_type = "Increasing")

declining_plot_data <- stratum_predictions %>% filter(stratum %in% declining_strata$stratum) %>%
  left_join(declining_strata %>% dplyr::select(stratum, difference), by = "stratum") %>% mutate(trajectory_type = "Decreasing")

combined_plot_data <- bind_rows(increasing_plot_data, declining_plot_data)
top_increasing_stratum <- increasing_strata %>% head(1) %>% pull(stratum)
top_declining_stratum <- declining_strata %>% head(1) %>% pull(stratum)

increasing_label_data <- increasing_plot_data %>% filter(Study_age == "1") %>% arrange(desc(difference)) %>%
  mutate(label = ifelse(stratum == top_increasing_stratum, paste0("*** ", stratum), as.character(stratum)))

declining_label_data <- declining_plot_data %>% filter(Study_age == "1") %>% arrange(difference) %>%
  mutate(label = ifelse(stratum == top_declining_stratum, paste0("*** ", stratum), as.character(stratum)))

combined_label_data <- bind_rows(increasing_label_data, declining_label_data)

combined_trajectory_plot <- ggplot() +
  geom_line(data = combined_plot_data, aes(x = Study_age_label, y = pred_mean, group = stratum, color = Income), linewidth = 1.0) +
  geom_text_repel(data = combined_label_data, aes(x = Study_age_label, y = pred_mean, label = label, color = Income),
                  direction = "y", nudge_x = -0.2, hjust = 1, segment.size = 0.2, segment.color = "gray50", min.segment.length = 0, 
                  size = 3.2, fontface = "bold", force = 2, box.padding = 0.4, point.padding = 0.2, show.legend = FALSE) +
  facet_wrap(~ trajectory_type, scales = "free_y", ncol = 2,
             labeller = labeller(trajectory_type = c("Increasing" = "Top 5 Increasing", "Decreasing" = "Top 5 Decreasing"))) +
  scale_color_manual(values = c("Low" = "#FF9999", "Middle" = "#66B2FF", "High" = "#66CC99"), name = "Income Level") +
  scale_y_continuous(breaks = breaks_width(0.02), labels = label_number(accuracy = 0.01)) +
  scale_x_discrete(limits = c("Adolescence (16 yrs)", "Early adulthood (33 yrs)"),
                   labels = c("Adolescence (16 yrs)" = "Adolescence", "Early adulthood (33 yrs)" = "Early adulthood")) +
  labs(caption = "*** indicates the stratum with the largest increase/decrease", x = "Study age", y = "Probability of good/very good health") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
        plot.caption = element_text(size = 10, hjust = 0.5, margin = margin(t = 10)), legend.position = "bottom",
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray90"), panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
        axis.text.y = element_text(size = 9), strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightyellow", color = "gray70"), plot.margin = margin(10, 60, 10, 10))

# 11. EXPORT RANKING TO EXCEL (Table S2)--------------------------------------------------
cat("\nExporting Datasets and Ranking...\n")
raw_data_export <- df_preds %>%
  group_by(stratum) %>%
  summarize(
    p_mean = mean(pred_well), p_lower = mean(pred_well_lower), p_upper = mean(pred_well_upper),
    n_total = n(), gender = first(gender), family = first(family_type), 
    educ = first(parental_educ), inc = first(income_tertiles), fr = first(Close_friends), .groups = "drop"
  ) %>% filter(n_total >= 5) %>%
  mutate(`Prob (95% CI)` = sprintf("%.2f (%.2f–%.2f)", p_mean, p_lower, p_upper), Rank = rank(desc(p_mean), ties.method = "first")) %>% arrange(Rank)

write.xlsx(raw_data_export, "outputs/Appendix_Full_Ranking.xlsx", overwrite = TRUE)

# 12. SAVE ALL FIGURES & DIAGNOSTICS -------------------------------------------
cat("\nSaving all figures to directory...\n")
ggsave("outputs/figures/Figure_1_Interaction.tiff", p1, width = 6, height = 4, dpi = 600, compression = "lzw")
ggsave("outputs/figures/Figure_2_Income.tiff", p2, width = 6, height = 4, dpi = 600, compression = "lzw")
ggsave("outputs/figures/Figure_3_Sex.tiff", p3, width = 6, height = 4, dpi = 600, compression = "lzw")
ggsave("outputs/figures/Figure_5_Education.tiff", p5, width = 8, height = 5, dpi = 600, compression = "lzw")
ggsave("outputs/figures/Figure_6_Friends.tiff", p6, width = 8, height = 5, dpi = 600, compression = "lzw")
ggsave("outputs/figures/Figure_Residuals.tiff", combined_residuals, width = 10, height = 7, dpi = 600, compression = "lzw")
ggsave("outputs/figures/Figure_Trajectory.tiff", combined_trajectory_plot, width = 6, height = 4, dpi = 600, compression = "lzw")

ggsave("outputs/diagnostics/Appendix_Trace_Null.tiff", trace_null, width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave("outputs/diagnostics/Appendix_Trace_Main.tiff", trace_main, width = 8, height = 6, dpi = 300, compression = "lzw")

