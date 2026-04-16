# =============================================================================
# mcap_populate_temp.R
# Calculates daily mean temperatures and midday temperature (MMD) throughout
# the experiment (acclimation to recovery), runs summary statistics, and
# produces a daily temperature line plot with phase background shading.
#
# Author : Ji Hoon Han
# Updated: 2026-04-12
#
# Input  (data/input/):
#   mcap_temp.csv    -- hourly temperature logger data
#                       DateTime: M/D/YYYY HH:MM
#                       E4, E1, E2: water bath (bin) sensors
#                       J1-J7: jar sensors
#                       bin_temp, Jar_temp: pre-calculated hourly averages
#                       Phase: acclimation / heating / recovery
#
# Outputs:
#   data/transformed/mcap_temp_daily.csv
#   data/transformed/mcap_temp_MMD.csv
#   figs/fig_temperature.png
# =============================================================================

library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Load and prepare
# -----------------------------------------------------------------------------
mcap_temp <- read.csv("data/input/mcap_temp.csv", sep = ",",
                      stringsAsFactors = FALSE) %>%
  mutate(
    DateTime = as.POSIXct(DateTime, format = "%m/%d/%Y %H:%M"),
    Date     = as.Date(DateTime),
    Phase    = tolower(trimws(Phase))
  ) %>%
  mutate(across(c(E4, E1, E2, J1:J7, bin_temp, Jar_temp),
                ~ suppressWarnings(as.numeric(na_if(as.character(.), "NA")))))

message("Loaded ", nrow(mcap_temp), " rows | ",
        "Date range: ", min(mcap_temp$Date), " to ", max(mcap_temp$Date))
message("Phases: ", paste(sort(unique(mcap_temp$Phase)), collapse = ", "))

# -----------------------------------------------------------------------------
# 2. Daily means
# -----------------------------------------------------------------------------
daily_summary <- mcap_temp %>%
  group_by(Date, Phase) %>%
  summarise(
    across(c(E4, E1, E2, J1:J7), ~ mean(.x, na.rm = TRUE)),
    bin_temp = mean(bin_temp, na.rm = TRUE),
    Jar_temp = mean(Jar_temp, na.rm = TRUE),
    .groups  = "drop"
  )

write.csv(daily_summary, "data/transformed/mcap_temp_daily.csv",
          row.names = FALSE)
message("Done — mcap_temp_daily.csv written: ", nrow(daily_summary), " days")

# -----------------------------------------------------------------------------
# 3. Midday temperature (MMD): mean of bin_temp 10am-2pm
# -----------------------------------------------------------------------------
mmd_summary <- mcap_temp %>%
  filter(format(DateTime, "%H") %in% c("10", "11", "12", "13", "14")) %>%
  group_by(Date, Phase) %>%
  summarise(MMD = mean(bin_temp, na.rm = TRUE), .groups = "drop")

write.csv(mmd_summary, "data/transformed/mcap_temp_MMD.csv",
          row.names = FALSE)
message("Done — mcap_temp_MMD.csv written: ", nrow(mmd_summary), " days")

# -----------------------------------------------------------------------------
# 4. Summary statistics for reporting
# -----------------------------------------------------------------------------
cat("\n=== Bin Temperature Summary (daily mean) ===\n")
print(summary(daily_summary$bin_temp))

cat("\n=== Bin Temperature by Phase ===\n")
daily_summary %>%
  group_by(Phase) %>%
  summarise(
    mean = round(mean(bin_temp, na.rm = TRUE), 2),
    min  = round(min(bin_temp,  na.rm = TRUE), 2),
    max  = round(max(bin_temp,  na.rm = TRUE), 2)
  ) %>%
  print()

cat("\n=== MMD Summary ===\n")
print(summary(mmd_summary$MMD))

cat("\n=== MMD by Phase ===\n")
mmd_summary %>%
  group_by(Phase) %>%
  summarise(
    mean = round(mean(MMD, na.rm = TRUE), 2),
    min  = round(min(MMD,  na.rm = TRUE), 2),
    max  = round(max(MMD,  na.rm = TRUE), 2)
  ) %>%
  print()

# -----------------------------------------------------------------------------
# 5. Statistical test: bin_temp vs Jar_temp
# -----------------------------------------------------------------------------
data_test <- daily_summary %>% filter(!is.na(bin_temp) & !is.na(Jar_temp))
avg_diff  <- mean(data_test$bin_temp - data_test$Jar_temp, na.rm = TRUE)
cat("\n=== Bin vs Jar Temperature ===\n")
cat("Average difference (bin - jar):", round(avg_diff, 3), "°C\n")

shapiro_p <- shapiro.test(data_test$bin_temp - data_test$Jar_temp)$p.value
if (shapiro_p > 0.05) {
  test_result <- t.test(data_test$bin_temp, data_test$Jar_temp, paired = TRUE)
  cat("Test: Paired t-test | p =", round(test_result$p.value, 4), "\n")
} else {
  test_result <- wilcox.test(data_test$bin_temp, data_test$Jar_temp, paired = TRUE)
  cat("Test: Wilcoxon signed-rank | p =", round(test_result$p.value, 4), "\n")
}

# -----------------------------------------------------------------------------
# 6. Daily mean ± SE for plotting (from hourly data)
# -----------------------------------------------------------------------------
daily_plot <- mcap_temp %>%
  group_by(Date, Phase) %>%
  summarise(
    mean_temp = mean(bin_temp, na.rm = TRUE),
    se_temp   = sd(bin_temp, na.rm = TRUE) / sqrt(sum(!is.na(bin_temp))),
    .groups   = "drop"
  ) %>%
  mutate(Phase = factor(Phase, levels = c("acclimation", "heating", "recovery")))

# -----------------------------------------------------------------------------
# 7. Phase background rectangles
# -----------------------------------------------------------------------------
phase_bounds <- daily_plot %>%
  group_by(Phase) %>%
  summarise(xmin = min(Date), xmax = max(Date), .groups = "drop")

phase_colors <- c(
  "acclimation" = "#d4e8f0",
  "heating"     = "#f7d9c4",
  "recovery"    = "#d4ecd4"
)

phase_labels <- c(
  "acclimation" = "Acclimation",
  "heating"     = "Heating",
  "recovery"    = "Recovery"
)

# -----------------------------------------------------------------------------
# 8. Plot
# -----------------------------------------------------------------------------
p <- ggplot(daily_plot, aes(x = Date, y = mean_temp)) +
  geom_rect(
    data = phase_bounds,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Phase),
    alpha = 0.3, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = phase_colors, labels = phase_labels, name = "Phase") +
  geom_ribbon(
    aes(ymin = mean_temp - se_temp, ymax = mean_temp + se_temp),
    fill = "grey40", alpha = 0.25
  ) +
  geom_line(color = "grey20", linewidth = 0.8) +
  geom_text(
    data = phase_bounds %>%
      mutate(x_mid = xmin + (xmax - xmin) / 2),
    aes(x = x_mid, y = Inf, label = phase_labels[as.character(Phase)]),
    vjust = 1.5, size = 3.5, fontface = "italic", inherit.aes = FALSE
  ) +
  labs(x = "Date", y = "Mean daily temperature (°C)") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p)

if (!dir.exists("figs")) dir.create("figs")
ggsave("figs/fig_temperature.png", p,
       width = 9, height = 4, dpi = 300, bg = "white")