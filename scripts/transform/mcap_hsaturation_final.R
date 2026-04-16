# =============================================================================
# mcap_hsaturation_final.R
# Calculates hsat, day length, relative hsat, and DSPI per nubbin per
# RLC day, using the individual nubbin Ek as the saturation threshold.
#
# Author : Ji Hoon Han
# Updated: 2026-04-12
#
# Inputs  (data/transformed/):
#   mcap_ek.csv      -- Date: YYYY-MM-DD, rlc_end_time: HH:MM:SS 24-hr
#   mcap_irrad.csv   -- date_time: YYYY-MM-DD HH:MM:SS UTC, NO index col
#
# Output (output/):
#   mcap_mean_hsaturation.csv
#
# Run order:
#   1. mc_phytotools.R            -> mcap_ek_alpha_normalized_2024.csv
#   2. mcap_populate_irradiance.R -> mcap_irrad.csv
#   3. mcap_populate_ek.R         -> mcap_ek.csv
#   4. THIS SCRIPT                -> mcap_mean_hsaturation.csv
# =============================================================================

library(dplyr)
library(purrr)

# -----------------------------------------------------------------------------
# 1. Load
# -----------------------------------------------------------------------------
ek <- read.csv("data/transformed/mcap_ek.csv",
               sep = ",", stringsAsFactors = FALSE)

irradiance <- read.csv("data/transformed/mcap_irrad.csv",
                       sep = ",", stringsAsFactors = FALSE)

# Drop row index column if present (old files written with row.names = TRUE)
if ("" %in% names(irradiance)) irradiance[["" ]] <- NULL
if ("X" %in% names(irradiance)) irradiance[["X"]] <- NULL

# Parse irradiance date_time
irradiance$date_time <- as.POSIXct(irradiance$date_time,
                                   format = "%Y-%m-%d %H:%M:%S")

# Sanity checks
message("ek Date sample: ",            paste(head(ek$Date, 2), collapse = ", "))
message("ek rlc_end_time sample: ",    paste(head(ek$rlc_end_time, 2), collapse = ", "))
message("irradiance columns: ",        paste(names(irradiance), collapse = ", "))
message("irradiance date_time range: ",
        as.character(min(irradiance$date_time)), " to ",
        as.character(max(irradiance$date_time)))

if (any(ek$Date == "NA" | is.na(ek$Date))) {
  stop("ek Date column has NAs — re-run mcap_populate_ek.R")
}
if (all(is.na(irradiance$date_time))) {
  stop("irradiance date_time is all NA — re-run mcap_populate_irradiance.R")
}

# Fixed experiment end (day after last RLC measurement 11/23/2024)
# Used as run_end for nubbins that did not survive to day 151
EXPERIMENT_END <- as.POSIXct("2024-11-24 00:00:00",
                             format = "%Y-%m-%d %H:%M:%S")

# -----------------------------------------------------------------------------
# 2. Pivot ek to wide format — one row per nubbin, columns per RLC day
# -----------------------------------------------------------------------------
hsat1 <- ek %>%
  filter(as.character(rlc_day) == "1") %>%
  dplyr::select(nubbinID, Date, trt1, trt2, genotype,
         ek_est, rlc_end_time, pmax, pmax_min) %>%
  rename(day1_ek = ek_est, day1_rlc_time = rlc_end_time,
         day1_pmax = pmax, day1_pmax_min = pmax_min)

hsat9 <- ek %>%
  filter(as.character(rlc_day) == "9") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day9_ek = ek_est, day9_rlc_time = rlc_end_time,
         day9_pmax = pmax, day9_pmax_min = pmax_min)

hsat13 <- ek %>%
  filter(as.character(rlc_day) == "13") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day13_ek = ek_est, day13_rlc_time = rlc_end_time,
         day13_pmax = pmax, day13_pmax_min = pmax_min)

hsat19 <- ek %>%
  filter(as.character(rlc_day) == "19") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day19_ek = ek_est, day19_rlc_time = rlc_end_time,
         day19_pmax = pmax, day19_pmax_min = pmax_min)

hsat27 <- ek %>%
  filter(as.character(rlc_day) == "27") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day27_ek = ek_est, day27_rlc_time = rlc_end_time,
         day27_pmax = pmax, day27_pmax_min = pmax_min)

hsat34 <- ek %>%
  filter(as.character(rlc_day) == "34") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day34_ek = ek_est, day34_rlc_time = rlc_end_time,
         day34_pmax = pmax, day34_pmax_min = pmax_min)

hsat40 <- ek %>%
  filter(as.character(rlc_day) == "40") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day40_ek = ek_est, day40_rlc_time = rlc_end_time,
         day40_pmax = pmax, day40_pmax_min = pmax_min)

hsat53 <- ek %>%
  filter(as.character(rlc_day) == "53") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day53_ek = ek_est, day53_rlc_time = rlc_end_time,
         day53_pmax = pmax, day53_pmax_min = pmax_min)

hsat62 <- ek %>%
  filter(as.character(rlc_day) == "62") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day62_ek = ek_est, day62_rlc_time = rlc_end_time,
         day62_pmax = pmax, day62_pmax_min = pmax_min)

hsat67 <- ek %>%
  filter(as.character(rlc_day) == "67") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day67_ek = ek_est, day67_rlc_time = rlc_end_time,
         day67_pmax = pmax, day67_pmax_min = pmax_min)

hsat72 <- ek %>%
  filter(as.character(rlc_day) == "72") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day72_ek = ek_est, day72_rlc_time = rlc_end_time,
         day72_pmax = pmax, day72_pmax_min = pmax_min)

hsat80 <- ek %>%
  filter(as.character(rlc_day) == "80") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day80_ek = ek_est, day80_rlc_time = rlc_end_time,
         day80_pmax = pmax, day80_pmax_min = pmax_min)

hsat86 <- ek %>%
  filter(as.character(rlc_day) == "86") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day86_ek = ek_est, day86_rlc_time = rlc_end_time,
         day86_pmax = pmax, day86_pmax_min = pmax_min)

hsat91 <- ek %>%
  filter(as.character(rlc_day) == "91") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day91_ek = ek_est, day91_rlc_time = rlc_end_time,
         day91_pmax = pmax, day91_pmax_min = pmax_min)

hsat94 <- ek %>%
  filter(as.character(rlc_day) == "94") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day94_ek = ek_est, day94_rlc_time = rlc_end_time,
         day94_pmax = pmax, day94_pmax_min = pmax_min)

hsat100 <- ek %>%
  filter(as.character(rlc_day) == "100") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day100_ek = ek_est, day100_rlc_time = rlc_end_time,
         day100_pmax = pmax, day100_pmax_min = pmax_min)

hsat112 <- ek %>%
  filter(as.character(rlc_day) == "112") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day112_ek = ek_est, day112_rlc_time = rlc_end_time,
         day112_pmax = pmax, day112_pmax_min = pmax_min)

hsat123 <- ek %>%
  filter(as.character(rlc_day) == "123") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day123_ek = ek_est, day123_rlc_time = rlc_end_time,
         day123_pmax = pmax, day123_pmax_min = pmax_min)

hsat136 <- ek %>%
  filter(as.character(rlc_day) == "136") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day136_ek = ek_est, day136_rlc_time = rlc_end_time,
         day136_pmax = pmax, day136_pmax_min = pmax_min)

hsat151 <- ek %>%
  filter(as.character(rlc_day) == "151") %>%
  dplyr::select(nubbinID, rlc_end_time, ek_est, pmax, pmax_min) %>%
  rename(day151_ek = ek_est, day151_rlc_time = rlc_end_time,
         day151_pmax = pmax, day151_pmax_min = pmax_min)

# -----------------------------------------------------------------------------
# 3. Join all RLC days — one wide row per nubbin
# -----------------------------------------------------------------------------
hsat_all <- hsat1 %>%
  left_join(hsat9,   by = "nubbinID") %>%
  left_join(hsat13,  by = "nubbinID") %>%
  left_join(hsat19,  by = "nubbinID") %>%
  left_join(hsat27,  by = "nubbinID") %>%
  left_join(hsat34,  by = "nubbinID") %>%
  left_join(hsat40,  by = "nubbinID") %>%
  left_join(hsat53,  by = "nubbinID") %>%
  left_join(hsat62,  by = "nubbinID") %>%
  left_join(hsat67,  by = "nubbinID") %>%
  left_join(hsat72,  by = "nubbinID") %>%
  left_join(hsat80,  by = "nubbinID") %>%
  left_join(hsat86,  by = "nubbinID") %>%
  left_join(hsat91,  by = "nubbinID") %>%
  left_join(hsat94,  by = "nubbinID") %>%
  left_join(hsat100, by = "nubbinID") %>%
  left_join(hsat112, by = "nubbinID") %>%
  left_join(hsat123, by = "nubbinID") %>%
  left_join(hsat136, by = "nubbinID") %>%
  left_join(hsat151, by = "nubbinID")

message("hsat_all rows: ", nrow(hsat_all))

# -----------------------------------------------------------------------------
# 4. Helper: calculate hsat over a window of days
# -----------------------------------------------------------------------------
calc_hsat_by_day <- function(run_irradiance, day1_date, days_to_consider,
                             threshold) {
  hsat_by_day      <- numeric(length(days_to_consider))
  daylight_minutes <- numeric(length(days_to_consider))
  rel_hsat_by_day  <- numeric(length(days_to_consider))
  
  base <- as.POSIXct(paste0(day1_date, " 00:00:00"),
                     format = "%Y-%m-%d %H:%M:%S")
  
  for (i in seq_along(days_to_consider)) {
    start <- base + (86400 * days_to_consider[i])
    end   <- start + 86399
    
    hsat_by_day[i] <- sum(
      run_irradiance$date_time > start &
        run_irradiance$date_time < end &
        run_irradiance$Epar > threshold
    )
    daylight_minutes[i] <- sum(
      run_irradiance$date_time > start &
        run_irradiance$date_time < end &
        run_irradiance$Epar >= 1
    )
    rel_hsat_by_day[i] <- ifelse(daylight_minutes[i] > 0,
                                 hsat_by_day[i] / daylight_minutes[i], NA)
  }
  list(hsat_by_day      = hsat_by_day,
       daylight_minutes = daylight_minutes,
       rel_hsat_by_day  = rel_hsat_by_day)
}

# -----------------------------------------------------------------------------
# 5. Main per-nubbin calculation
#    run_end uses EXPERIMENT_END if day151_rlc_time is NA (nubbin died early)
# -----------------------------------------------------------------------------
calculate_hsat <- function(day1_date, day1_rlc_time, day151_rlc_time,
                           day1_ek,  day9_ek,  day13_ek, day19_ek,
                           day27_ek, day34_ek, day40_ek, day53_ek,
                           day62_ek, day67_ek, day72_ek, day80_ek,
                           day86_ek, day91_ek, day94_ek, day100_ek,
                           day112_ek, day123_ek, day136_ek, day151_ek) {
  
  # Start from midnight of day1 so full p1 window is captured
  run_start <- as.POSIXct(paste0(day1_date, " 00:00:00"),
                          format = "%Y-%m-%d %H:%M:%S")
  
  # Use fixed experiment end if nubbin did not survive to day 151
  if (is.na(day151_rlc_time) || day151_rlc_time == "NA") {
    run_end <- EXPERIMENT_END
  } else {
    run_end <- as.POSIXct(paste(day1_date, day151_rlc_time),
                          format = "%Y-%m-%d %H:%M:%S") +
      152 * 86400
  }
  
  run_irr <- irradiance[irradiance$date_time >= run_start &
                          irradiance$date_time <  run_end, ]
  ivs <- list(
    p1   = calc_hsat_by_day(run_irr, day1_date, 0:3,     day1_ek),
    p9   = calc_hsat_by_day(run_irr, day1_date, 6:10,    day9_ek),
    p13  = calc_hsat_by_day(run_irr, day1_date, 10:14,   day13_ek),
    p19  = calc_hsat_by_day(run_irr, day1_date, 16:20,   day19_ek),
    p27  = calc_hsat_by_day(run_irr, day1_date, 24:28,   day27_ek),
    p34  = calc_hsat_by_day(run_irr, day1_date, 31:35,   day34_ek),
    p40  = calc_hsat_by_day(run_irr, day1_date, 38:42,   day40_ek),
    p53  = calc_hsat_by_day(run_irr, day1_date, 51:55,   day53_ek),
    p62  = calc_hsat_by_day(run_irr, day1_date, 60:64,   day62_ek),
    p67  = calc_hsat_by_day(run_irr, day1_date, 65:69,   day67_ek),
    p72  = calc_hsat_by_day(run_irr, day1_date, 70:74,   day72_ek),
    p80  = calc_hsat_by_day(run_irr, day1_date, 78:82,   day80_ek),
    p86  = calc_hsat_by_day(run_irr, day1_date, 84:88,   day86_ek),
    p91  = calc_hsat_by_day(run_irr, day1_date, 89:93,   day91_ek),
    p94  = calc_hsat_by_day(run_irr, day1_date, 92:96,   day94_ek),
    p100 = calc_hsat_by_day(run_irr, day1_date, 98:102,  day100_ek),
    p112 = calc_hsat_by_day(run_irr, day1_date, 110:114, day112_ek),
    p123 = calc_hsat_by_day(run_irr, day1_date, 121:125, day123_ek),
    p136 = calc_hsat_by_day(run_irr, day1_date, 134:138, day136_ek),
    p151 = calc_hsat_by_day(run_irr, day1_date, 148:151, day151_ek)
  )
  
  hv  <- lapply(ivs, function(x) mean(x$hsat_by_day,     na.rm = TRUE))
  dlv <- lapply(ivs, function(x) mean(x$daylight_minutes, na.rm = TRUE))
  rv  <- lapply(ivs, function(x) mean(x$rel_hsat_by_day,  na.rm = TRUE))
  
  list(
    individual_hsat       = hv,
    individual_day_length = dlv,
    individual_rel_hsat   = rv,
    acclimation_hsat = mean(c(hv$p1, hv$p9), na.rm = TRUE),
    heating_hsat     = mean(c(hv$p9, hv$p13, hv$p19, hv$p27, hv$p34),
                            na.rm = TRUE),
    recovery_hsat    = mean(c(hv$p34, hv$p40, hv$p53, hv$p62, hv$p67,
                              hv$p72, hv$p80, hv$p86, hv$p91, hv$p94,
                              hv$p100, hv$p112, hv$p123, hv$p136, hv$p151),
                            na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# 6. Run for all nubbins
# -----------------------------------------------------------------------------
message("Running hsat for ", nrow(hsat_all), " nubbins...")

results <- pmap(list(
  hsat_all$Date,          hsat_all$day1_rlc_time, hsat_all$day151_rlc_time,
  hsat_all$day1_ek,       hsat_all$day9_ek,       hsat_all$day13_ek,
  hsat_all$day19_ek,      hsat_all$day27_ek,      hsat_all$day34_ek,
  hsat_all$day40_ek,      hsat_all$day53_ek,      hsat_all$day62_ek,
  hsat_all$day67_ek,      hsat_all$day72_ek,      hsat_all$day80_ek,
  hsat_all$day86_ek,      hsat_all$day91_ek,      hsat_all$day94_ek,
  hsat_all$day100_ek,     hsat_all$day112_ek,     hsat_all$day123_ek,
  hsat_all$day136_ek,     hsat_all$day151_ek
), calculate_hsat)

# -----------------------------------------------------------------------------
# 7. Extract results
# -----------------------------------------------------------------------------
hsat_all <- hsat_all %>%
  mutate(
    hsat_p1         = map_dbl(results, ~ .x$individual_hsat[["p1"]]),
    day_length_p1   = map_dbl(results, ~ .x$individual_day_length[["p1"]]),
    rel_hsat_p1     = map_dbl(results, ~ .x$individual_rel_hsat[["p1"]]),
    hsat_p9         = map_dbl(results, ~ .x$individual_hsat[["p9"]]),
    day_length_p9   = map_dbl(results, ~ .x$individual_day_length[["p9"]]),
    rel_hsat_p9     = map_dbl(results, ~ .x$individual_rel_hsat[["p9"]]),
    hsat_p13        = map_dbl(results, ~ .x$individual_hsat[["p13"]]),
    day_length_p13  = map_dbl(results, ~ .x$individual_day_length[["p13"]]),
    rel_hsat_p13    = map_dbl(results, ~ .x$individual_rel_hsat[["p13"]]),
    hsat_p19        = map_dbl(results, ~ .x$individual_hsat[["p19"]]),
    day_length_p19  = map_dbl(results, ~ .x$individual_day_length[["p19"]]),
    rel_hsat_p19    = map_dbl(results, ~ .x$individual_rel_hsat[["p19"]]),
    hsat_p27        = map_dbl(results, ~ .x$individual_hsat[["p27"]]),
    day_length_p27  = map_dbl(results, ~ .x$individual_day_length[["p27"]]),
    rel_hsat_p27    = map_dbl(results, ~ .x$individual_rel_hsat[["p27"]]),
    hsat_p34        = map_dbl(results, ~ .x$individual_hsat[["p34"]]),
    day_length_p34  = map_dbl(results, ~ .x$individual_day_length[["p34"]]),
    rel_hsat_p34    = map_dbl(results, ~ .x$individual_rel_hsat[["p34"]]),
    hsat_p40        = map_dbl(results, ~ .x$individual_hsat[["p40"]]),
    day_length_p40  = map_dbl(results, ~ .x$individual_day_length[["p40"]]),
    rel_hsat_p40    = map_dbl(results, ~ .x$individual_rel_hsat[["p40"]]),
    hsat_p53        = map_dbl(results, ~ .x$individual_hsat[["p53"]]),
    day_length_p53  = map_dbl(results, ~ .x$individual_day_length[["p53"]]),
    rel_hsat_p53    = map_dbl(results, ~ .x$individual_rel_hsat[["p53"]]),
    hsat_p62        = map_dbl(results, ~ .x$individual_hsat[["p62"]]),
    day_length_p62  = map_dbl(results, ~ .x$individual_day_length[["p62"]]),
    rel_hsat_p62    = map_dbl(results, ~ .x$individual_rel_hsat[["p62"]]),
    hsat_p67        = map_dbl(results, ~ .x$individual_hsat[["p67"]]),
    day_length_p67  = map_dbl(results, ~ .x$individual_day_length[["p67"]]),
    rel_hsat_p67    = map_dbl(results, ~ .x$individual_rel_hsat[["p67"]]),
    hsat_p72        = map_dbl(results, ~ .x$individual_hsat[["p72"]]),
    day_length_p72  = map_dbl(results, ~ .x$individual_day_length[["p72"]]),
    rel_hsat_p72    = map_dbl(results, ~ .x$individual_rel_hsat[["p72"]]),
    hsat_p80        = map_dbl(results, ~ .x$individual_hsat[["p80"]]),
    day_length_p80  = map_dbl(results, ~ .x$individual_day_length[["p80"]]),
    rel_hsat_p80    = map_dbl(results, ~ .x$individual_rel_hsat[["p80"]]),
    hsat_p86        = map_dbl(results, ~ .x$individual_hsat[["p86"]]),
    day_length_p86  = map_dbl(results, ~ .x$individual_day_length[["p86"]]),
    rel_hsat_p86    = map_dbl(results, ~ .x$individual_rel_hsat[["p86"]]),
    hsat_p91        = map_dbl(results, ~ .x$individual_hsat[["p91"]]),
    day_length_p91  = map_dbl(results, ~ .x$individual_day_length[["p91"]]),
    rel_hsat_p91    = map_dbl(results, ~ .x$individual_rel_hsat[["p91"]]),
    hsat_p94        = map_dbl(results, ~ .x$individual_hsat[["p94"]]),
    day_length_p94  = map_dbl(results, ~ .x$individual_day_length[["p94"]]),
    rel_hsat_p94    = map_dbl(results, ~ .x$individual_rel_hsat[["p94"]]),
    hsat_p100       = map_dbl(results, ~ .x$individual_hsat[["p100"]]),
    day_length_p100 = map_dbl(results, ~ .x$individual_day_length[["p100"]]),
    rel_hsat_p100   = map_dbl(results, ~ .x$individual_rel_hsat[["p100"]]),
    hsat_p112       = map_dbl(results, ~ .x$individual_hsat[["p112"]]),
    day_length_p112 = map_dbl(results, ~ .x$individual_day_length[["p112"]]),
    rel_hsat_p112   = map_dbl(results, ~ .x$individual_rel_hsat[["p112"]]),
    hsat_p123       = map_dbl(results, ~ .x$individual_hsat[["p123"]]),
    day_length_p123 = map_dbl(results, ~ .x$individual_day_length[["p123"]]),
    rel_hsat_p123   = map_dbl(results, ~ .x$individual_rel_hsat[["p123"]]),
    hsat_p136       = map_dbl(results, ~ .x$individual_hsat[["p136"]]),
    day_length_p136 = map_dbl(results, ~ .x$individual_day_length[["p136"]]),
    rel_hsat_p136   = map_dbl(results, ~ .x$individual_rel_hsat[["p136"]]),
    hsat_p151       = map_dbl(results, ~ .x$individual_hsat[["p151"]]),
    day_length_p151 = map_dbl(results, ~ .x$individual_day_length[["p151"]]),
    rel_hsat_p151   = map_dbl(results, ~ .x$individual_rel_hsat[["p151"]])
  )

# -----------------------------------------------------------------------------
# 8. DSPI = hsat * pmax_min
# -----------------------------------------------------------------------------
hsat_all <- hsat_all %>%
  mutate(
    dspi_pmax1   = hsat_p1   * day1_pmax_min,
    dspi_pmax9   = hsat_p9   * day9_pmax_min,
    dspi_pmax13  = hsat_p13  * day13_pmax_min,
    dspi_pmax19  = hsat_p19  * day19_pmax_min,
    dspi_pmax27  = hsat_p27  * day27_pmax_min,
    dspi_pmax34  = hsat_p34  * day34_pmax_min,
    dspi_pmax40  = hsat_p40  * day40_pmax_min,
    dspi_pmax53  = hsat_p53  * day53_pmax_min,
    dspi_pmax62  = hsat_p62  * day62_pmax_min,
    dspi_pmax67  = hsat_p67  * day67_pmax_min,
    dspi_pmax72  = hsat_p72  * day72_pmax_min,
    dspi_pmax80  = hsat_p80  * day80_pmax_min,
    dspi_pmax86  = hsat_p86  * day86_pmax_min,
    dspi_pmax91  = hsat_p91  * day91_pmax_min,
    dspi_pmax94  = hsat_p94  * day94_pmax_min,
    dspi_pmax100 = hsat_p100 * day100_pmax_min,
    dspi_pmax112 = hsat_p112 * day112_pmax_min,
    dspi_pmax123 = hsat_p123 * day123_pmax_min,
    dspi_pmax136 = hsat_p136 * day136_pmax_min,
    dspi_pmax151 = hsat_p151 * day151_pmax_min
  )

# -----------------------------------------------------------------------------
# 9. Write
# -----------------------------------------------------------------------------
write.csv(hsat_all, "output/mcap_mean_hsaturation.csv", row.names = FALSE)