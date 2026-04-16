# =============================================================================
# build_rlc_data.R
# Joins mcap_ek.csv (long format PE parameters) with mcap_mean_hsaturation.csv
# (wide format hsat metrics) to produce rlc_data.csv for final analyses.
#
# Author : Ji Hoon Han
# Updated: 2026-04-12
#
# Inputs:
#   data/transformed/mcap_ek.csv
#   output/mcap_mean_hsaturation.csv
#
# Output:
#   data/input/rlc_data.csv
#
# Run order:
#   1. mc_phytotools.R            -> mcap_ek_alpha_normalized_2024.csv
#   2. mcap_populate_irradiance.R -> mcap_irrad.csv
#   3. mcap_populate_ek.R         -> mcap_ek.csv
#   4. mcap_hsaturation_final.R   -> mcap_mean_hsaturation.csv
#   5. THIS SCRIPT                -> rlc_data.csv
# =============================================================================

library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load
# -----------------------------------------------------------------------------
ek   <- read.csv("data/transformed/mcap_ek.csv",
                 sep = ",", stringsAsFactors = FALSE)
hsat <- read.csv("output/mcap_mean_hsaturation.csv",
                 sep = ",", stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 2. Pivot hsat wide -> long
# -----------------------------------------------------------------------------
rlc_days <- c(1, 9, 13, 19, 27, 34, 40, 53, 62, 67,
              72, 80, 86, 91, 94, 100, 112, 123, 136, 151)

hsat_long <- lapply(rlc_days, function(d) {
  data.frame(
    nubbinID   = as.character(hsat$nubbinID),
    rlc_day    = as.character(d),
    hsat       = hsat[[paste0("hsat_p",      d)]],
    day_length = hsat[[paste0("day_length_p", d)]],
    rel_hsat   = hsat[[paste0("rel_hsat_p",  d)]],
    dspi_pmax  = hsat[[paste0("dspi_pmax",   d)]],
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# -----------------------------------------------------------------------------
# 3. Prepare mcap_ek
#    - drop short-code 'phase' column first, then rename phase_full -> phase
#    - rename rlc_end_time -> rlc_time
#    - lowercase trt1 to match original (ctl, nut, dead)
#    - drop rank
# -----------------------------------------------------------------------------
ek_prep <- ek %>%
  dplyr::select(-any_of(c("phase", "rank"))) %>%   # drop short phase and rank first
  rename(
    rlc_time = rlc_end_time,
    phase    = phase_full                    # now rename phase_full -> phase
  ) %>%
  mutate(
    nubbinID = as.character(nubbinID),
    rlc_day  = as.character(rlc_day),
    trt1     = tolower(trt1)                 # ctl, nut, dead
  )

# -----------------------------------------------------------------------------
# 4. Join hsat metrics
# -----------------------------------------------------------------------------
rlc_data <- ek_prep %>%
  left_join(hsat_long, by = c("nubbinID", "rlc_day"))

# -----------------------------------------------------------------------------
# 5. Select and reorder columns to match original rlc_data.csv
# -----------------------------------------------------------------------------
rlc_data <- rlc_data %>%
  dplyr::select(
    Date, uid, nubbinID, trt1, trt2, genotype,
    rlc_time, rlc_day, rlc_order, phase, tp2,
    pmax, pmax_min, deltaNPQ, maxNPQ_Ypoint1,
    alpha_est, alpha_st_err, alpha_t, alpha_p,
    ek_est, ek_st_err, ek_t, ek_p,
    QuanYield, F0, Fm,
    hsat, dspi_pmax, day_length, rel_hsat
  )

# -----------------------------------------------------------------------------
# 6. Write
# -----------------------------------------------------------------------------
write.csv(rlc_data, "data/input/rlc_data.csv", row.names = FALSE)
