# mcap_combine_surv_bw.R
# Joins transformed BW growth data with survivorship to produce surv_bw.csv,
# the final input for growth.R.
#
# ── Inputs ────────────────────────────────────────────────────────────────────
#   data/transformed/bw_growth.csv   output of mcap_populate_bw.R
#     columns: genotype, nubbin_num, trt2, bw1–bw5, bw5_dead,
#              calc_acclimation, calc_heating
#
#   data/input/raw_surv.csv
#     columns: genotype, nubbin_num, trt2, surv_days, event
#     surv_days : days survived from 2024-08-16 (recovery start)
#     event     : 0 = survived to day 100, 1 = died during recovery,
#                 "heat" = died during heating (before recovery)
#
# ── calc_recovery formula ─────────────────────────────────────────────────────
#   rate   = ((bw5 - bw4) * 1000 / bw4) / surv_days
#   growth = sign(rate) * |rate|^(1/3)
#   "D"    where bw5 was missing (coral died before final weighing)
#
# ── Output ────────────────────────────────────────────────────────────────────
#   data/input/surv_bw.csv   ready for growth.R
#
# ── growth.R filter (applied downstream, not here) ───────────────────────────
#   filter(!calc_recovery %in% "D", !trt2 %in% "dead", surv_days >= 10)
#   -> retains 27 nubbins for final analysis

library(dplyr)

# ── Helper ────────────────────────────────────────────────────────────────────
signed_cbrt <- function(x) {
  ifelse(is.na(x), NA_real_, sign(x) * abs(x)^(1/3))
}

# ── Load ──────────────────────────────────────────────────────────────────────
bw_growth <- read.csv("data/transformed/bw_growth.csv", sep = ",",
                      colClasses = "character") %>%
  mutate(
    nubbin_num = as.integer(nubbin_num),
    bw4        = as.numeric(bw4),
    bw5_num    = suppressWarnings(as.numeric(bw5)),   # numeric bw5 for calc
    bw5_dead   = as.logical(bw5_dead)
  )

raw_surv <- read.csv("data/input/raw_surv.csv", sep = ",") %>%
  mutate(surv_days = as.numeric(surv_days)) %>%
  dplyr::select(nubbin_num, surv_days, event)

# ── Join ──────────────────────────────────────────────────────────────────────
combined <- bw_growth %>%
  left_join(raw_surv, by = "nubbin_num")

# ── calc_recovery ─────────────────────────────────────────────────────────────
combined <- combined %>%
  mutate(
    rec_rate     = ifelse(
      is.na(bw5_num) | is.na(bw4) | is.na(surv_days) | surv_days == 0,
      NA_real_,
      ((bw5_num - bw4) * 1000 / bw4) / surv_days
    ),
    calc_recovery = ifelse(
      bw5_dead,
      "D",
      as.character(signed_cbrt(rec_rate))
    )
  )

# ── Final output columns (matching surv_bw.csv structure) ────────────────────
surv_bw <- combined %>%
  dplyr::select(genotype, nubbin_num, trt2,
         bw1, bw2, bw3, bw4, bw5,
         surv_days, event,
         calc_acclimation, calc_heating, calc_recovery)

# ── Export ────────────────────────────────────────────────────────────────────
write.csv(surv_bw, "data/input/surv_bw.csv", row.names = FALSE)

# ── Summary ───────────────────────────────────────────────────────────────────
n_analysis <- surv_bw %>%
  filter(calc_recovery != "D",
         trt2 != "dead",
         as.numeric(surv_days) >= 10) %>%
  nrow()