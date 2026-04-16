#########################################################
# GAMM + Derivative Analysis for RLC Metrics
# - Fits GAMMs with treatment-specific smooths
# - Computes R²m / R²c from the mixed-effects part
# - Uses simple LMMs for overall treatment posthoc tests
# - Computes derivatives (relative rate of change)
# - 8-panel curve and derivative figures
#########################################################

# -----------------------------
# Packages
# -----------------------------
library(tidyverse)
library(mgcv)
library(gamm4)
library(gratia)
library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)
library(cowplot)
library(performance)
library(rlang)

`%>%` <- magrittr::`%>%`

setwd("C:/R_projects/src/photosynthesis")
# -----------------------------
# Read and prepare data
# -----------------------------
raw <- read.csv("data/input/rlc_data.csv")

df_base <- raw %>%
  filter(
    phase == "recovery",
    trt2 %in% c("CTL", "NH4", "NO3")
  ) %>%
  mutate(
    tp2          = as.numeric(gsub("[^0-9]", "", tp2)),
    tp2_shifted  = tp2,
    tp2_centered = tp2 - mean(tp2, na.rm = TRUE),
    Treatment    = factor(trt2, levels = c("CTL", "NH4", "NO3")),
    genotype     = factor(genotype)
  )
# -----------------------------
# Metrics, labels, colours
# -----------------------------
metrics <- c(
  "F0", "Fm", "QuanYield", "deltaNPQ",
  "pmax", "ek_est", "hsat", "dspi_pmax"
)
metric_names <- list(
  pmax      = expression(P[max]),
  deltaNPQ  = expression(Delta * NPQ),
  QuanYield = expression(F[v] / F[m]),
  F0        = expression(F[0]),
  Fm        = expression(F[m]),
  ek_est    = expression(E[k]),
  hsat      = expression(H[sat]),
  dspi_pmax = expression(DSPI)
)
metric_units <- list(
  pmax      = expression(mu * mol ~ electrons ~ m^{-2} ~ s^{-1}),
  deltaNPQ  = "AU",
  QuanYield = "AU",
  F0        = "AU",
  Fm        = "AU",
  ek_est    = expression(mu * mol ~ photons ~ m^{-2} ~ s^{-1}),
  hsat      = "minutes",
  dspi_pmax = expression(mol ~ m^{-2} ~ d^{-1})
)
color_set <- c(
  "CTL" = "dodgerblue3",
  "NH4" = "firebrick3",
  "NO3" = "goldenrod2"
)

deriv_y_label <- "Relative rate of change (per day)"

# -----------------------------
# Helper: R² label
# -----------------------------
get_r2_label <- function(model) {
  
  if (inherits(model, "lm")) {
    r2 <- tryCatch(summary(model)$r.squared, error = function(e) NA_real_)
    
    if (!is.na(r2)) {
      r2_fmt <- formatC(r2, digits = 3, format = "f")
      return(paste0("R² = ", r2_fmt, "\nR²c = NA"))
    } else {
      return("R² = NA\nR²c = NA")
    }
  }
  out <- tryCatch(r2_nakagawa(model), error = function(e) NULL)
  if (is.null(out)) return("R²m = NA\nR²c = NA")
  
  fmt_val <- function(v) {
    if (is.null(v) || length(v) == 0 || !is.numeric(v) || all(is.na(v))) {
      "NA"
    } else {
      formatC(v[1], digits = 3, format = "f")
    }
  }
  if (is.list(out) && (!is.null(out$R2_marginal) || !is.null(out$R2_conditional))) {
    r2m <- fmt_val(out$R2_marginal)
    r2c <- fmt_val(out$R2_conditional)
    return(paste0("R²m = ", r2m, "\nR²c = ", r2c))
  }
  if (is.numeric(out)) {
    r2 <- fmt_val(out[1])
    return(paste0("R² = ", r2, "\nR²c = NA"))
  }
  "R²m = NA\nR²c = NA"
}
# -----------------------------
# Containers
# -----------------------------
curve_plots   <- vector("list", length(metrics)); names(curve_plots)  <- metrics
deriv_plots   <- vector("list", length(metrics)); names(deriv_plots)  <- metrics
gamm_models   <- vector("list", length(metrics)); names(gamm_models)  <- metrics
deriv_models  <- vector("list", length(metrics)); names(deriv_models) <- metrics
posthoc_tests <- list()

smooth_tables      <- list()
posthoc_level_tbls <- list()
posthoc_deriv_tbls <- list()

# ======================================================
# MAIN LOOP
# ======================================================
for (metric in metrics) {
  # --------------------------
  # Metric-specific dataframe
  # --------------------------
  df <- df_base %>%
    transmute(
      tp2_shifted,
      tp2_centered,
      Treatment,
      genotype,
      value = if (metric == "dspi_pmax") !!sym(metric) * 1e-6 else !!sym(metric),
      log_value = log1p(value)
    ) %>%
    filter(!is.na(value))
  # --------------------------
  # GAMM fit
  # --------------------------
  mod <- gamm4(
    log_value ~ s(tp2_centered, by = Treatment, k = 7) + Treatment,
    random = ~(1 | genotype),
    data = df
  )
  gamm_models[[metric]] <- mod
  # --------------------------
  # Smooth-term summary
  # --------------------------
  st <- summary(mod$gam)$s.table %>%
    as.data.frame() %>%
    rownames_to_column("smooth_term") %>%
    mutate(
      metric = metric,
      p_value = `p-value`,
      .before = 1
    )
  
  smooth_tables[[metric]] <- st
  
  cat("\n--- GAMM smooth-term summary for", metric, "---\n")
  print(st)
  
  # --------------------------
  # Predictions for curve plot
  # --------------------------
  newdat <- expand.grid(
    tp2_centered = seq(min(df$tp2_centered), max(df$tp2_centered), length.out = 100),
    Treatment = levels(df$Treatment)
  )
  
  pred <- predict(mod$gam, newdata = newdat, se.fit = TRUE)
  
  newdat <- newdat %>%
    mutate(
      fit = expm1(pred$fit),
      lower = expm1(pred$fit - 1.96 * pred$se.fit),
      upper = expm1(pred$fit + 1.96 * pred$se.fit),
      tp2_shifted = tp2_centered + mean(df$tp2_shifted)
    )
  # --------------------------
  # Baseline from first acclimation day
  # --------------------------
  df_accl <- raw %>%
    filter(
      phase == "acclimation",
      trt2 %in% c("CTL", "NH4", "NO3")
    ) %>%
    mutate(tp2_num = as.numeric(gsub("[^0-9]", "", tp2)))
  
  first_accl_day <- suppressWarnings(min(df_accl$tp2_num, na.rm = TRUE))
  baseline_val <- NA_real_
  
  if (is.finite(first_accl_day)) {
    baseline_val <- df_accl %>%
      filter(tp2_num == first_accl_day) %>%
      summarise(
        baseline = mean(
          if (metric == "dspi_pmax") !!sym(metric) * 1e-6 else !!sym(metric),
          na.rm = TRUE
        )
      ) %>%
      pull(baseline)
  }
  
  if (!is.finite(baseline_val)) {
    first_rec_day <- min(df$tp2_shifted, na.rm = TRUE)
    
    baseline_val <- df %>%
      filter(abs(tp2_shifted - first_rec_day) < 1e-9) %>%
      summarise(baseline = mean(value, na.rm = TRUE)) %>%
      pull(baseline)
  }
  # --------------------------
  # Simple LMM for overall treatment posthoc
  # --------------------------
  lmm_simple <- lmer(log_value ~ Treatment + (1 | genotype), data = df)
  emm_simple <- emmeans(lmm_simple, ~ Treatment)
  comp_simple <- pairs(emm_simple, adjust = "holm")
  
  posthoc_tests[[metric]] <- comp_simple
  
  cat("\n--- LMM (log_value) ANOVA for", metric, "---\n")
  print(anova(lmm_simple))
  
  cat("\n--- LMM pairwise treatment contrasts (Holm) for", metric, "---\n")
  print(comp_simple)
  
  level_df <- as.data.frame(comp_simple) %>%
    mutate(metric = metric, .before = 1)
  
  posthoc_level_tbls[[metric]] <- level_df
  
  # --------------------------
  # R² label for curve plot
  # --------------------------
  r2_label_curve <- get_r2_label(mod$mer)
  
  x_pos <- min(df$tp2_shifted) +
    0.8 * (max(df$tp2_shifted) - min(df$tp2_shifted))
  
  y_pos <- min(df$value, na.rm = TRUE) +
    0.9 * (max(df$value, na.rm = TRUE) - min(df$value, na.rm = TRUE))
  
  # --------------------------
  # Curve plot
  # --------------------------
  p_curve <- ggplot(df, aes(x = tp2_shifted, y = value, color = Treatment)) +
    geom_point(alpha = 0.3) +
    geom_line(
      data = newdat,
      aes(x = tp2_shifted, y = fit, color = Treatment),
      linewidth = 1.2
    ) +
    geom_ribbon(
      data = newdat,
      aes(x = tp2_shifted, ymin = lower, ymax = upper, fill = Treatment),
      alpha = 0.1,
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = baseline_val,
      linetype = "dashed",
      color = "black",
      linewidth = 0.7
    ) +
    annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = r2_label_curve,
      size = 3.5,
      hjust = 1,
      vjust = 1,
      fontface = "plain"
    ) +
    labs(
      title = metric_names[[metric]],
      y = metric_units[[metric]],
      x = "Days"
    ) +
    scale_color_manual(values = color_set, name = NULL) +
    scale_fill_manual(values = color_set, name = NULL) +
    guides(
      color = guide_legend(direction = "horizontal"),
      fill = guide_legend(direction = "horizontal")
    ) +
    theme_minimal()
  
  curve_plots[[metric]] <- p_curve
  
  # --------------------------
  # Derivatives from GAMM
  # --------------------------
  deriv <- derivatives(
    mod$gam,
    select = "tp2_centered",
    partial_match = TRUE
  ) %>%
    mutate(
      Recovery_Day = tp2_centered + mean(df$tp2_shifted),
      rate_back  = exp(.derivative) - 1,
      lower_back = exp(.lower_ci) - 1,
      upper_back = exp(.upper_ci) - 1,
      Treatment = case_when(
        grepl("NH4", .smooth) ~ "NH4",
        grepl("NO3", .smooth) ~ "NO3",
        TRUE ~ "CTL"
      )
    ) %>%
    mutate(DayFactor = factor(round(Recovery_Day)))
  
  # --------------------------
  # Derivative model
  # --------------------------
  rate_lmm <- tryCatch(
    {
      m <- lmer(rate_back ~ Treatment + (1 | DayFactor), data = deriv)
      if (performance::check_singularity(m)) stop("singular")
      m
    },
    error = function(e) {
      lm(rate_back ~ Treatment, data = deriv)
    }
  )
  deriv_models[[metric]] <- rate_lmm
  
  posthoc_deriv <- emmeans(rate_lmm, pairwise ~ Treatment, adjust = "holm")
  posthoc_tests[[paste0(metric, "_deriv")]] <- posthoc_deriv
  
  cat("\n--- Derivative model ANOVA for", metric, "---\n")
  print(anova(rate_lmm))
  
  cat("\n--- Derivative pairwise treatment contrasts (Holm) for", metric, "---\n")
  print(posthoc_deriv)
  
  deriv_df <- as.data.frame(posthoc_deriv$contrasts) %>%
    mutate(metric = metric, .before = 1)
  
  posthoc_deriv_tbls[[metric]] <- deriv_df
  
  # --------------------------
  # Compact letter display for derivatives
  # --------------------------
  letters_obj <- cld(
    posthoc_deriv$emmeans,
    Letters = letters,
    adjust = "holm"
  )
  letters_df <- as.data.frame(letters_obj)
  letters_deriv <- letters_df[, c("Treatment", ".group")]
  
  last_day <- deriv %>%
    group_by(Treatment) %>%
    filter(Recovery_Day == max(Recovery_Day)) %>%
    summarise(
      x = max(Recovery_Day),
      y = mean(rate_back),
      .groups = "drop"
    ) %>%
    left_join(letters_deriv, by = "Treatment")
  
  r2_label_deriv <- get_r2_label(rate_lmm)
  
  # --------------------------
  # Derivative plot
  # --------------------------
  p_deriv <- ggplot(deriv, aes(x = Recovery_Day, y = rate_back, color = Treatment)) +
    geom_point(alpha = 0.3) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(
      aes(ymin = lower_back, ymax = upper_back, fill = Treatment),
      alpha = 0.05
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = r2_label_deriv,
      size = 3,
      hjust = 1.3,
      vjust = 1.15
    ) +
    geom_text(
      data = last_day,
      aes(x = x, y = y, label = .group, color = Treatment),
      vjust = -2,
      fontface = "bold",
      size = 4,
      inherit.aes = FALSE
    ) +
    labs(
      title = metric_names[[metric]],
      x = "Days",
      y = deriv_y_label
    ) +
    scale_color_manual(values = color_set, name = NULL) +
    scale_fill_manual(values = color_set, name = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
  
  deriv_plots[[metric]] <- p_deriv
}

# ======================================================
# Combine plots + legend
# ======================================================
curve_nolegend <- lapply(curve_plots, function(p) p + theme(legend.position = "none"))
deriv_nolegend <- lapply(deriv_plots, function(p) p + theme(legend.position = "none"))

p_with_legend <- curve_plots[[1]] + theme(legend.position = "bottom")
g <- ggplotGrob(p_with_legend)
legend_index <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
legend <- g$grobs[[legend_index]]

curve_grid <- cowplot::plot_grid(
  plotlist = curve_nolegend,
  ncol = 4
)
final_gamm_plot <- cowplot::plot_grid(
  curve_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)
deriv_grid <- cowplot::plot_grid(
  plotlist = deriv_nolegend,
  ncol = 4
)
final_deriv_plot <- cowplot::plot_grid(
  deriv_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)
print(final_gamm_plot)
print(final_deriv_plot)

# ======================================================
# Save outputs
# ======================================================
if (!dir.exists("figs")) dir.create("figs")
if (!dir.exists("output")) dir.create("output")

ggsave(
  "figs/fig_gamm_timeseries.png",
  final_gamm_plot,
  width = 10,
  height = 5.5,
  dpi = 600
)
ggsave(
  "figs/fig_gamm_derivatives.png",
  final_deriv_plot,
  width = 10,
  height = 5.5,
  dpi = 600
)
smooth_all <- bind_rows(smooth_tables)
write.csv(
  smooth_all,
  "output/gamm_smooth_tables_all_metrics.csv",
  row.names = FALSE
)
level_all <- bind_rows(posthoc_level_tbls)
write.csv(
  level_all,
  "output/gamm_level_posthoc_all_metrics.csv",
  row.names = FALSE
)
deriv_all <- bind_rows(posthoc_deriv_tbls)
write.csv(
  deriv_all,
  "output/gamm_derivative_posthoc_all_metrics.csv",
  row.names = FALSE
)
