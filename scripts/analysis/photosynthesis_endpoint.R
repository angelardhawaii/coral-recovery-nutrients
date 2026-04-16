library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(cowplot)
library(multcompView)

# ==========================================
# Endpoint photosynthesis analysis (R100)
# ==========================================
setwd("C:/R_projects/src/photosynthesis")
# ------------------------------
# 1. Load and prepare data
# ------------------------------
photo_all <- read.csv("data/input/rlc_data.csv") %>%
  mutate(
    tp2 = as.character(tp2),
    tp2 = trimws(tp2),
    Treatment = factor(trt2, levels = c("CTL", "NH4", "NO3")),
    genotype = factor(genotype),
    dspi_pmax = dspi_pmax / 1e6,
    F0 = as.numeric(F0),
    Fm = as.numeric(Fm)
  )
photo <- photo_all %>%
  filter(
    phase == "recovery",
    trt2 %in% c("CTL", "NH4", "NO3")
  )

# ------------------------------
# 2. Metric setup
# ------------------------------
metrics <- c("F0", "Fm", "QuanYield", "deltaNPQ", "pmax", "ek_est", "hsat", "dspi_pmax")

metric_names <- list(
  F0 = expression(F[0]),
  Fm = expression(F[m]),
  QuanYield = expression(F[v] / F[m]),
  deltaNPQ = expression(Delta * NPQ),
  pmax = expression(P[max]),
  ek_est = expression(E[k]),
  hsat = expression(H[sat]),
  dspi_pmax = expression(DSPI)
)
metric_units <- list(
  F0 = "AU",
  Fm = "AU",
  QuanYield = "AU",
  deltaNPQ = "AU",
  pmax = expression(mu * mol ~ electrons ~ m^{-2} ~ s^{-1}),
  ek_est = expression(mu * mol ~ photons ~ m^{-2} ~ s^{-1}),
  hsat = "minutes",
  dspi_pmax = expression(mol ~ m^{-2} ~ d^{-1})
)
trt_cols <- c(
  "CTL" = "dodgerblue3",
  "NH4" = "firebrick3",
  "NO3" = "goldenrod2"
)
# ------------------------------
# 3. Legend
# ------------------------------
legend_plot <- ggplot(
  photo %>% filter(tp2 == "100"),
  aes(x = Treatment, y = as.numeric(pmax), color = Treatment)
) +
  geom_point(size = 3) +
  scale_color_manual(
    values = trt_cols,
    labels = c("Control", expression(NH[4]^"+"), expression(NO[3]^"-")),
    name = ""
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

g <- ggplotGrob(legend_plot)
legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
legend_grob <- g$grobs[[legend_index]]

# ------------------------------
# 4. Loop through metrics
# ------------------------------
plots_day100_nox <- list()

for (m in metrics) {
  
  rec_m <- photo %>%
    mutate(value = as.numeric(.data[[m]]))
  
  end_100 <- rec_m %>%
    filter(tp2 == "100" & !is.na(value))
  
  if (nrow(end_100) == 0) next
  
  # --------------------------
  # Model
  # --------------------------
  fit <- lmer(
    value ~ Treatment + (1 | genotype),
    data = end_100,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  # --------------------------
  # Post-hoc letters
  # Keep your original a,b,c style
  # --------------------------
  emm <- multcomp::cld(emmeans(fit, ~ Treatment), adjust = "holm")
  
  posthoc_letters <- emm %>%
    as.data.frame() %>%
    dplyr::select(Treatment, .group) %>%
    mutate(.group = trimws(.group)) %>%
    distinct() %>%
    mutate(group_letter = tolower(letters[match(.group, unique(.group))]))
  
  y_max <- max(end_100$value, na.rm = TRUE) * 1.1
  
  summary_pos <- data.frame(
    Treatment = levels(end_100$Treatment),
    y_pos = y_max
  ) %>%
    left_join(posthoc_letters, by = "Treatment")
  
  # --------------------------
  # Print model output
  # --------------------------
  cat("\n==============================\n")
  cat("### Metric:", m, "###\n")
  cat("==============================\n")
  
  cat("\n--- LMM Summary ---\n")
  print(summary(fit))
  
  cat("\n--- Estimated Marginal Means ---\n")
  print(emmeans(fit, ~ Treatment))
  
  cat("\n--- Pairwise Comparisons (Holm Adjusted) ---\n")
  print(pairs(emmeans(fit, ~ Treatment), adjust = "holm"))
  
  cat("\n--- Post-hoc Group Letters ---\n")
  print(emm)
  
  # --------------------------
  # Pre-bleaching baseline mean
  # --------------------------
  rec_all_m <- photo_all %>%
    mutate(value = as.numeric(.data[[m]]))
  
  mean_i <- rec_all_m %>%
    filter(tp2 == "i" & !is.na(value)) %>%
    summarise(mean_value = mean(value, na.rm = TRUE)) %>%
    pull(mean_value)
  
  # --------------------------
  # Plot
  # --------------------------
  p_box_nox <- ggplot(end_100, aes(x = Treatment, y = value, color = Treatment)) +
    geom_boxplot(fill = NA, linewidth = 0.7, outlier.shape = NA) +
    geom_jitter(
      aes(fill = Treatment),
      width = 0.15,
      shape = 21,
      size = 2,
      alpha = 0.7,
      stroke = 0.6
    )
  
  if (!is.na(mean_i)) {
    p_box_nox <- p_box_nox +
      geom_hline(
        yintercept = mean_i,
        linetype = "dashed",
        color = "black",
        linewidth = 0.8
      )
  }
  
  p_box_nox <- p_box_nox +
    geom_text(
      data = summary_pos,
      aes(x = Treatment, y = y_pos, label = group_letter),
      inherit.aes = FALSE,
      size = 5,
      color = "black"
    ) +
    labs(
      title = metric_names[[m]],
      y = metric_units[[m]],
      x = ""
    ) +
    scale_color_manual(
      values = trt_cols
    ) +
    scale_fill_manual(
      values = trt_cols
    ) +
    scale_x_discrete(
      labels = c(
        "CTL" = "Control",
        "NH4" = expression(NH[4]^"+"),
        "NO3" = expression(NO[3]^"-")
      )
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 12)
    )
  
  if (m == "hsat") {
    p_box_nox <- p_box_nox + coord_cartesian(ylim = c(150, 750))
  }
  
  plots_day100_nox[[m]] <- p_box_nox
}

# ------------------------------
# 5. Combine figure
# ------------------------------
figure_combined_nox <- plot_grid(
  plot_grid(plotlist = plots_day100_nox, ncol = 4, labels = "auto"),
  legend_grob,
  ncol = 1,
  rel_heights = c(1, 0.1)
)
print(figure_combined_nox)