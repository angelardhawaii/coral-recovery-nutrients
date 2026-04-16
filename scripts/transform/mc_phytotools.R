# =============================================================================
# mc_phytotools.R
# Fits P-E curves (Webb model) to PAM fluorometer RLC data using phytotools,
# then joins experimental metadata to produce mcap_ek_alpha_normalized_2024.csv
#
# Author : Ji Hoon Han
# Created: 2024-11-07
# Updated: 2026-04-12
#
# Inputs  (data/input/):
#   all_pam_clean3.csv       -- cleaned PAM RLC measurements
#   rlc_day.csv              -- date -> rlc_day / phase_full / tp2 lookup
#   mcap_2024_meta.csv       -- per-uid treatment metadata
#
# Output  (data/transformed/):
#   mcap_ek_alpha_normalized_2024.csv
#
# Dependencies:
#   phytotools 1.0   -- removed from CRAN 2022-06-20, installed from archive
#   insol      1.2.2 -- removed from CRAN 2023-02-27, installed from archive
#   raster           -- available on CRAN (required by insol)
#
# System requirement:
#   Rtools45 (R 4.5.x):
#   https://cran.rstudio.com/bin/windows/Rtools/rtools45/rtools.html
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Install archived dependencies
#    Install order matters: raster -> insol -> phytotools
# -----------------------------------------------------------------------------
if (!requireNamespace("raster", quietly = TRUE)) {
  install.packages("raster")
}
if (!requireNamespace("insol", quietly = TRUE)) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/insol/insol_1.2.2.tar.gz",
    repos = NULL, type = "source"
  )
}
if (!requireNamespace("phytotools", quietly = TRUE)) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/phytotools/phytotools_1.0.tar.gz",
    repos = NULL, type = "source"
  )
}

# -----------------------------------------------------------------------------
# 2. Libraries
# -----------------------------------------------------------------------------
library(phytotools)
library(hash)
library(dplyr)
library(lubridate)

# -----------------------------------------------------------------------------
# 3. Load data
# -----------------------------------------------------------------------------
mc_phyto_nona  <- read.csv("data/input/pam_clean.csv",
                           sep = ",", stringsAsFactors = FALSE)
rlc_day_assign <- read.csv("data/input/rlc_day.csv",
                           sep = ",", strip.white = TRUE,
                           stringsAsFactors = FALSE)
meta_data      <- read.csv("data/input/mcap_2024_meta.csv",
                           sep = ",", stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 4. Prepare PAM data
#    Input dates are M/D/YYYY — convert to YYYY-MM-DD using mdy()
# -----------------------------------------------------------------------------
mc_phyto_nona$Date <- format(mdy(mc_phyto_nona$Date), "%Y-%m-%d")
mc_phyto_nona$uid  <- paste(mc_phyto_nona$Date, mc_phyto_nona$ID, sep = "_")
mc_phyto_nona$NPQ  <- as.numeric(mc_phyto_nona$NPQ)

# Fv/Fm: Y(II) at Epar == 0
mc_phyto_nona <- transform(mc_phyto_nona,
                           QuanYield = ifelse(Epar == 0, Y.II., NA))

# Remove zero/negative rETR rows before curve fitting
selectedData <- subset(mc_phyto_nona, rETR > 0)

# -----------------------------------------------------------------------------
# 5. Build rlc_day / phase_full / tp2 lookup hashes
#    Input dates are M/D/YYYY — convert to YYYY-MM-DD for hash keys
# -----------------------------------------------------------------------------
rlc_day_assign$Date <- format(mdy(rlc_day_assign$Date), "%Y-%m-%d")

hash_rlc_day    <- hash(rlc_day_assign$Date, rlc_day_assign$rlc_day)
hash_phase_full <- hash(rlc_day_assign$Date, rlc_day_assign$phase_full)
hash_tp2        <- hash(rlc_day_assign$Date, rlc_day_assign$tp2)

# -----------------------------------------------------------------------------
# 6. Unique RLC identifiers and per-curve date lookups
# -----------------------------------------------------------------------------
uniqueIds <- unique(selectedData$uid)
n         <- length(uniqueIds)
dates     <- substr(uniqueIds, 1, 10)    # YYYY-MM-DD portion

rlc_days_by_date   <- array(NA_real_,      dim = n)
phase_full_by_date <- array(NA_character_, dim = n)
tp2_by_date        <- array(NA_character_, dim = n)

for (i in seq_len(n)) {
  rlc_days_by_date[i]   <- hash_rlc_day[[dates[i]]]
  phase_full_by_date[i] <- hash_phase_full[[dates[i]]]
  tp2_by_date[i]        <- hash_tp2[[dates[i]]]
}

# -----------------------------------------------------------------------------
# 7. Per-curve NPQ metrics and RLC end time
# -----------------------------------------------------------------------------
maxNPQ_ypoint1 <- array(NA_real_,      dim = n)
minNPQ         <- array(NA_real_,      dim = n)
rlc_end_times  <- array(NA_character_, dim = n)

for (i in seq_len(n)) {
  lc <- subset(selectedData, uid == uniqueIds[i])
  maxNPQ_ypoint1[i] <- as.numeric(max(subset(lc, Y.II. > 0.1)$NPQ, na.rm = TRUE))
  minNPQ[i]         <- as.numeric(max(subset(lc, Epar == 25)$NPQ))
  rlc_end_times[i]  <- max(lc$Time)
}

deltaNPQ_ypoint1 <- maxNPQ_ypoint1 - minNPQ

# -----------------------------------------------------------------------------
# 8. Fit Webb P-E model (normalized by Y(II))
# -----------------------------------------------------------------------------
alpha <- array(NA_real_, dim = c(n, 4))
colnames(alpha) <- c("est", "st_err", "t", "p")
ek    <- array(NA_real_, dim = c(n, 4))
colnames(ek)    <- c("est", "st_err", "t", "p")

for (i in seq_len(n)) {
  Epar  <- selectedData$Epar[selectedData$uid  == uniqueIds[i]]
  y.II  <- selectedData$Y.II.[selectedData$uid == uniqueIds[i]]
  myfit <- fitWebb(Epar, y.II, normalize = TRUE)
  alpha[i, ] <- myfit$alpha
  ek[i, ]    <- myfit$ek
}

pmax <- round(alpha[, 1] * ek[, 1], digits = 2)

# -----------------------------------------------------------------------------
# 9. Extract Fv/Fm (Epar == 0) and rlc_order (Epar == 685)
# -----------------------------------------------------------------------------
first_row_of_rlc <- subset(mc_phyto_nona, Epar == 0)
last_row_rlc     <- subset(mc_phyto_nona, Epar == 685)

# -----------------------------------------------------------------------------
# 10. Assemble result data frame
#     Date is YYYY-MM-DD character — safe for downstream scripts
# -----------------------------------------------------------------------------
result_df <- data.frame(
  Date           = substr(uniqueIds, 1, 10),
  rlc_end_time   = rlc_end_times,
  nubbinID       = substr(uniqueIds, 12, 13),
  uid            = uniqueIds,
  rlc_day        = rlc_days_by_date,
  phase_full     = phase_full_by_date,
  tp2            = tp2_by_date,
  pmax           = pmax,
  maxNPQ_Ypoint1 = maxNPQ_ypoint1,
  deltaNPQ       = deltaNPQ_ypoint1,
  alpha_est      = alpha[, 1],
  alpha_st_err   = alpha[, 2],
  alpha_t        = alpha[, 3],
  alpha_p        = alpha[, 4],
  ek_est         = ek[, 1],
  ek_st_err      = ek[, 2],
  ek_t           = ek[, 3],
  ek_p           = ek[, 4],
  QuanYield      = first_row_of_rlc$QuanYield,
  F0             = first_row_of_rlc$F,
  Fm             = first_row_of_rlc$Fm.,
  rlc_order      = last_row_rlc$rlc_order,
  stringsAsFactors = FALSE
)

write.csv(result_df,
          "data/transformed/mcap_ek_alpha_normalized_2024.csv",
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 11. Join treatment metadata (no pret1/pret2)
# -----------------------------------------------------------------------------
mcap_data <- read.csv("data/transformed/mcap_ek_alpha_normalized_2024.csv",
                      sep = ",", stringsAsFactors = FALSE)

meta <- meta_data %>%
  dplyr::select(uid, genotype, phase, trt1, trt2) %>%
  mutate(
    genotype = as.factor(genotype),
    phase    = as.factor(phase),
    trt1     = as.factor(trt1),
    trt2     = as.factor(trt2)
  )

mcap_data <- mcap_data %>%
  left_join(meta, by = "uid")

# Rank within day and group into ~15-min rlc_order blocks
mcap_data <- mcap_data %>%
  group_by(Date) %>%
  mutate(rank = rank(rlc_order)) %>%
  ungroup() %>%
  mutate(rlc_order = floor((rank - 1) / 3) + 1)

write.csv(mcap_data,
          "data/transformed/mcap_ek_alpha_normalized_2024.csv",
          row.names = FALSE)