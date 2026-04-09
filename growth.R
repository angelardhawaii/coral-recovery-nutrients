library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcompView)
library(gridExtra)

# ===============================
# Growth Analysis
# ===============================

# -------------------------------
# 1. Load and prepare data
# -------------------------------
data_growth <- read.csv("data/input/surv_bw.csv") %>%
  filter(!calc_recovery %in% "D", !trt2 %in% "dead", surv_days >= 10) %>%
  mutate(
    trt2        = factor(trt2, levels = c("CTL", "NH4", "NO3")),
    Acclimation = as.numeric(calc_acclimation),
    Heating     = as.numeric(calc_heating),
    Recovery    = as.numeric(calc_recovery)
  )

# -------------------------------
# 2. Descriptive statistics
# -------------------------------
desc_stats_g <- data_growth %>%
  pivot_longer(
    cols = c("Acclimation", "Heating", "Recovery"),
    names_to = "Phase",
    values_to = "Growth"
  ) %>%
  group_by(trt2, Phase) %>%
  summarise(
    N      = n(),
    Mean   = round(mean(Growth, na.rm = TRUE), 3),
    SD     = round(sd(Growth, na.rm = TRUE), 3),
    SE     = round(sd(Growth, na.rm = TRUE) / sqrt(n()), 3),
    Median = round(median(Growth, na.rm = TRUE), 3),
    Min    = round(min(Growth, na.rm = TRUE), 3),
    Max    = round(max(Growth, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(desc_stats_g)

# -------------------------------
# 3. Phase-wise ANOVA + Tukey HSD
# -------------------------------
do_anova <- function(phase) {
  model <- aov(reformulate("trt2", response = phase), data = data_growth)
  tk <- TukeyHSD(model)
  cld <- multcompLetters4(model, tk)
  
  cld_df <- data_growth %>%
    group_by(trt2) %>%
    summarise(
      mean  = mean(.data[[phase]], na.rm = TRUE),
      quant = quantile(.data[[phase]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean))
  
  cld_df$cld <- as.data.frame.list(cld$trt2)$Letters
  
  list(
    model = model,
    tukey = tk,
    letters = cld_df
  )
}

res_a <- do_anova("Acclimation")
res_h <- do_anova("Heating")
res_r <- do_anova("Recovery")

cat("\n==============================\n")
cat("Acclimation ANOVA\n")
cat("==============================\n")
print(summary(res_a$model))
print(res_a$tukey)
print(res_a$letters)

cat("\n==============================\n")
cat("Heating ANOVA\n")
cat("==============================\n")
print(summary(res_h$model))
print(res_h$tukey)
print(res_h$letters)

cat("\n==============================\n")
cat("Recovery ANOVA\n")
cat("==============================\n")
print(summary(res_r$model))
print(res_r$tukey)
print(res_r$letters)

# -------------------------------
# 4. Optional phase comparison LMM
#    (kept from your original code)
# -------------------------------
data_long <- data_growth %>%
  pivot_longer(
    cols = c("calc_acclimation", "calc_heating", "calc_recovery"),
    names_to = "phase",
    values_to = "measurement",
    names_prefix = "calc_"
  ) %>%
  mutate(
    phase = factor(phase, levels = c("acclimation", "heating", "recovery")),
    measurement = as.numeric(measurement),
    nubbin_num = as.factor(nubbin_num)
  )

model_phase <- lmer(measurement ~ phase + (1 | nubbin_num), data = data_long)

cat("\n==============================\n")
cat("Phase LMM\n")
cat("==============================\n")
print(summary(model_phase))
print(emmeans(model_phase, pairwise ~ phase, adjust = "tukey"))

# -------------------------------
# 5. Plot settings
# -------------------------------
trt2_color <- c(
  CTL = "dodgerblue3",
  NH4 = "firebrick3",
  NO3 = "goldenrod2"
)

global_y_min <- min(
  data_growth$Acclimation,
  data_growth$Heating,
  data_growth$Recovery,
  na.rm = TRUE
)

global_y_max <- max(
  data_growth$Acclimation,
  data_growth$Heating,
  data_growth$Recovery,
  na.rm = TRUE
)

# -------------------------------
# 6. Phase-wise boxplot function
# -------------------------------
plot_phase <- function(df, phase, title, ylab, letters_df) {
  ggplot(df, aes(x = trt2, y = .data[[phase]], fill = trt2)) +
    geom_boxplot(alpha = 0.75) +
    scale_fill_manual(values = trt2_color) +
    labs(title = title, x = "", y = ylab) +
    theme_bw(base_size = 16) +
    theme(
      legend.position  = "none",
      panel.border     = element_blank(),
      panel.background = element_rect(fill = "grey100", color = NA),
      axis.title       = element_text(size = 16),
      axis.text        = element_text(size = 14),
      plot.title       = element_text(size = 16, face = "bold")
    ) +
    ylim(global_y_min, global_y_max) +
    geom_text(
      data = letters_df,
      aes(label = cld, x = trt2, y = global_y_max - 0.5),
      inherit.aes = FALSE,
      vjust = -2.5,
      size = 5
    )
}

# -------------------------------
# 7. Three-panel phase figure
# -------------------------------
p_acclimation <- plot_phase(
  data_growth,
  "Acclimation",
  "Acclimation",
  expression(Growth ~ (mm %.% d^-1)),
  res_a$letters
)

p_heating <- plot_phase(
  data_growth,
  "Heating",
  "Heating",
  "",
  res_h$letters
)

p_recovery_phase <- plot_phase(
  data_growth,
  "Recovery",
  "Recovery",
  "",
  res_r$letters
)

combined_plot <- grid.arrange(
  p_acclimation,
  p_heating,
  p_recovery_phase,
  ncol = 3
)

print(combined_plot)

# -------------------------------
# 8. Recovery-focused figure
# -------------------------------
mean_acclimation <- mean(data_growth$Acclimation, na.rm = TRUE)

p_recovery <- ggplot(data_growth, aes(x = trt2, y = Recovery, color = trt2)) +
  geom_boxplot(fill = NA, linewidth = 0.65, outlier.shape = NA) +
  geom_jitter(
    aes(fill = trt2),
    width = 0.15,
    shape = 21,
    size = 2.5,
    alpha = 0.6,
    stroke = 0.3
  ) +
  geom_hline(
    yintercept = mean_acclimation,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  ) +
  scale_color_manual(values = trt2_color) +
  scale_fill_manual(values = trt2_color) +
  scale_x_discrete(
    labels = c(
      "CTL" = "Control",
      "NH4" = "NH4+",
      "NO3" = "NO3-"
    )
  ) +
  labs(
    x = "",
    y = expression(Growth ~ (mm %.% d^-1))
  ) +
  geom_text(
    data = res_r$letters,
    aes(label = cld, x = trt2, y = global_y_max - 0.9),
    inherit.aes = FALSE,
    vjust = -1,
    size = 6
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position  = "none",
    panel.border     = element_blank(),
    panel.background = element_rect(fill = "grey100", color = NA),
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 14),
    axis.ticks.x     = element_blank(),
    axis.title.y     = element_text(size = 16),
    axis.text.y      = element_text(size = 14)
  ) +
  ylim(global_y_min, global_y_max)

print(p_recovery)

# -------------------------------
# 9. Optional export
# -------------------------------
# write.csv(desc_stats_g, "outputs/growth_descriptive_stats.csv", row.names = FALSE)
# ggsave("growth_phases.png", combined_plot, width = 8, height = 5, dpi = 300, bg = "white")
# ggsave("p_recovery.png", p_recovery, width = 8, height = 7, dpi = 300, bg = "white")