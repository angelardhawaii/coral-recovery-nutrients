# mcap_populate_bw.R
# Transforms raw buoyant weight (BW) measurements into per-phase growth rates.
#
# ── Input ─────────────────────────────────────────────────────────────────────
#   data/input/raw_bw.csv
#     columns : genotype, nubbin_num, trt2,
#               2024-06-25, 2024-07-03, 2024-07-28, 2024-08-16, 2024-11-23
#     "D"     : coral had already died before that weighing date
#
# ── Measurement dates ─────────────────────────────────────────────────────────
#   bw1  2024-06-25  acclimation start
#   bw2  2024-07-03  acclimation end / heating start
#   bw3  2024-07-28  heating end / recovery start
#   bw4  2024-08-16  recovery mid-point
#   bw5  2024-11-23  recovery end
#
# ── Conversion formula ────────────────────────────────────────────────────────
#   Step 1 — Dry skeletal weight:  DW = BW * 1.54
#            (aragonite density correction; factor cancels in the growth ratio)
#   Step 2 — Per-day specific growth rate (mg g⁻¹ day⁻¹):
#            rate = ((BW_final - BW_initial) * 1000 / BW_initial) / days
#   Step 3 — Signed cube-root transform (preserves sign of negative growth):
#            growth = sign(rate) * |rate|^(1/3)
#
# ── Phase durations ───────────────────────────────────────────────────────────
#   calc_acclimation : bw1 -> bw2  /  10 days
#   calc_heating     : bw2 -> bw3  /  26 days
#   calc_recovery    : bw4 -> bw5  /  individual surv_days (joined in next script)
#
# NOTE: calc_recovery is NOT calculated here because it requires surv_days,
#       which comes from raw_surv.csv. It is computed in mcap_combine_surv_bw.R.
#
# ── Output ────────────────────────────────────────────────────────────────────
#   data/transformed/bw_growth.csv

library(dplyr)

# ── Phase durations (days) ────────────────────────────────────────────────────
DAYS_ACC  <- 10   # bw1 -> bw2
DAYS_HEAT <- 26   # bw2 -> bw3

# ── Helpers ───────────────────────────────────────────────────────────────────
signed_cbrt <- function(x) {
  ifelse(is.na(x), NA_real_, sign(x) * abs(x)^(1/3))
}

calc_growth <- function(bw_f, bw_i, days) {
  ifelse(
    is.na(bw_f) | is.na(bw_i),
    NA_real_,
    ((bw_f - bw_i) * 1000 / bw_i) / days
  )
}

# ── Load ──────────────────────────────────────────────────────────────────────
bw_raw <- read.csv("data/input/raw_bw.csv", sep = ",",
                   check.names = FALSE, colClasses = "character") %>%
  mutate(nubbin_num = as.integer(nubbin_num))

# Track "D" flags before converting to numeric
bw2_dead <- bw_raw[["2024-07-03"]]  == "D"
bw3_dead <- bw_raw[["2024-07-28"]]  == "D"
bw5_dead <- bw_raw[["2024-11-23"]]  == "D"

# Convert BW columns to numeric ("D" -> NA)
bw_growth <- bw_raw %>%
  mutate(across(c(`2024-06-25`, `2024-07-03`, `2024-07-28`,
                  `2024-08-16`, `2024-11-23`),
                ~ suppressWarnings(as.numeric(.x)))) %>%
  rename(bw1 = `2024-06-25`,
         bw2 = `2024-07-03`,
         bw3 = `2024-07-28`,
         bw4 = `2024-08-16`,
         bw5 = `2024-11-23`) %>%
  mutate(
    # Per-day growth rates (mg g⁻¹ day⁻¹), cube-root transformed
    calc_acclimation = signed_cbrt(calc_growth(bw2, bw1, DAYS_ACC)),
    calc_heating     = signed_cbrt(calc_growth(bw3, bw2, DAYS_HEAT)),
    # calc_recovery computed in mcap_combine_surv_bw.R (needs surv_days)

    # Restore "D" strings where bw was missing (for downstream filtering)
    bw2 = ifelse(bw2_dead, "D", as.character(bw2)),
    bw3 = ifelse(bw3_dead, "D", as.character(bw3)),
    bw5 = ifelse(bw5_dead, "D", as.character(bw5)),
    calc_acclimation = ifelse(bw2_dead, "D", as.character(calc_acclimation)),
    calc_heating     = ifelse(bw3_dead, "D", as.character(calc_heating)),

    # Carry bw5_dead flag for calc_recovery in next script
    bw5_dead = bw5_dead
  ) %>%
  dplyr::select(genotype, nubbin_num, trt2,
         bw1, bw2, bw3, bw4, bw5, bw5_dead,
         calc_acclimation, calc_heating)

# ── Export ────────────────────────────────────────────────────────────────────
write.csv(bw_growth, "data/transformed/bw_growth.csv", row.names = FALSE)