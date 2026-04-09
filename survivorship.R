library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

# =========================
# 1. Load and clean data
# =========================
data_surv <- read.csv("data/input/surv_bw.csv") %>%
  filter(
    !is.na(surv_days),
    !is.na(event),
    !is.na(trt2),
    surv_days > 0
  ) %>%
  mutate(
    trt2 = factor(trt2, levels = c("CTL", "NH4", "NO3")),
    event = as.integer(event)
  )
# =========================
# 2. Descriptive statistics
# =========================
desc_stats <- data_surv %>%
  group_by(trt2) %>%
  summarise(
    n_total = n(),
    n_events = sum(event == 1),
    n_censored = sum(event == 0),
    mortality_rate = round(100 * mean(event == 1), 1),
    mean_survival_days = round(mean(surv_days), 1),
    median_survival_days = round(median(surv_days), 1),
    .groups = "drop"
  )
print(desc_stats)

# =========================
# 3. Survival object
# =========================
surv_object <- Surv(time = data_surv$surv_days, event = data_surv$event)

# =========================
# 4. Kaplan-Meier model
# =========================
km_fit <- survfit(surv_object ~ trt2, data = data_surv)

# Overall log-rank test
logrank_test <- survdiff(surv_object ~ trt2, data = data_surv)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)

print(logrank_test)
cat("Overall log-rank p-value:", signif(logrank_p, 3), "\n")

# Pairwise log-rank tests
pairwise_results <- pairwise_survdiff(
  Surv(surv_days, event) ~ trt2,
  data = data_surv,
  p.adjust.method = "holm"
)
print(pairwise_results)

# =========================
# 5. Cox proportional hazards model
# =========================
cox_model <- coxph(surv_object ~ trt2, data = data_surv)
cox_summary <- summary(cox_model)
print(cox_summary)

# Proportional hazards assumption
ph_test <- cox.zph(cox_model)
print(ph_test)

# Hazard ratio table
hr_table <- data.frame(
  term = rownames(cox_summary$coefficients),
  HR = round(cox_summary$conf.int[, "exp(coef)"], 2),
  lower_CI = round(cox_summary$conf.int[, "lower .95"], 2),
  upper_CI = round(cox_summary$conf.int[, "upper .95"], 2),
  p_value = signif(cox_summary$coefficients[, "Pr(>|z|)"], 3),
  row.names = NULL
)
print(hr_table)

# Median survival
median_survival <- surv_median(km_fit)
print(median_survival)

# =========================
# 6. Kaplan-Meier plot with risk table
# =========================
km_plot <- ggsurvplot(
  fit = km_fit,
  data = data_surv,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.12,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = TRUE,
  risk.table.col = "strata",
  censor = TRUE,
  xlim = c(0, 100),
  break.time.by = 20,
  palette = c("dodgerblue3", "salmon", "lightgoldenrod"),
  legend.title = "",
  legend.labs = c("Control", "NH4+", "NO3-"),
  xlab = "Recovery day",
  ylab = "Survival probability",
  ggtheme = theme_classic(base_size = 16),
  tables.theme = theme_cleantable(base_size = 12)
)

# Main panel formatting
km_plot$plot <- km_plot$plot +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13)
  )

# Risk table formatting
km_plot$table <- km_plot$table +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12)
  )
print(km_plot)

# =========================


