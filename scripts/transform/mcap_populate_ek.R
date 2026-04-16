# =============================================================================
# mcap_populate_ek.R
# Reads mcap_ek_alpha_normalized_2024.csv, applies correct column types,
# adds pmax_min, and writes mcap_ek.csv for use in the hsat pipeline.
#
# Author : Ji Hoon Han
# Updated: 2026-04-12
#
# Input  (data/transformed/):
#   mcap_ek_alpha_normalized_2024.csv
#     Date format: YYYY-MM-DD (written by mc_phytotools.R)
#     rlc_end_time format: H:MM:SS AM/PM (12-hour, from PAM output)
#
# Output (data/transformed/):
#   mcap_ek.csv
#     Date format: YYYY-MM-DD (preserved as character)
#     rlc_end_time format: HH:MM:SS (24-hour, for hsat pipeline)
# =============================================================================

library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load
# -----------------------------------------------------------------------------
ek <- read.csv("data/transformed/mcap_ek_alpha_normalized_2024.csv",
               sep = ",", stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 2. Date — verify YYYY-MM-DD format, keep as character
#    IMPORTANT: do NOT use as.Date() — causes NA on write/read cycle
# -----------------------------------------------------------------------------
message("Date sample from normalized file: ", paste(head(ek$Date, 3), collapse = ", "))

if (any(is.na(ek$Date)) || !all(nchar(trimws(ek$Date)) == 10)) {
  stop("Date column has unexpected values. Check mcap_ek_alpha_normalized_2024.csv ",
       "was produced by the updated mc_phytotools.R")
}

ek$Date <- as.character(ek$Date)

# -----------------------------------------------------------------------------
# 3. rlc_end_time — convert from 12-hour AM/PM to 24-hour HH:MM:SS
#    Input:  "7:44:34 AM"  ->  Output: "07:44:34"
# -----------------------------------------------------------------------------
ek$rlc_end_time <- format(
  as.POSIXct(ek$rlc_end_time, format = "%I:%M:%S %p"),
  "%H:%M:%S"
)

# -----------------------------------------------------------------------------
# 4. Type conversions
# -----------------------------------------------------------------------------
ek$rlc_day    <- as.factor(ek$rlc_day)
ek$nubbinID   <- as.factor(ek$nubbinID)
ek$rlc_order  <- as.factor(ek$rlc_order)
ek$phase      <- as.factor(ek$phase)
ek$phase_full <- as.factor(ek$phase_full)
ek$tp2        <- as.character(ek$tp2)
ek$genotype   <- as.factor(ek$genotype)
ek$trt1       <- as.character(ek$trt1)
ek$trt2       <- as.character(ek$trt2)

# -----------------------------------------------------------------------------
# 5. Add pmax_min (pmax per second -> per minute)
# -----------------------------------------------------------------------------
ek$pmax_min <- ek$pmax * 60

# -----------------------------------------------------------------------------
# 6. Write
# -----------------------------------------------------------------------------
write.csv(ek, "data/transformed/mcap_ek.csv", row.names = FALSE)