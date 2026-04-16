# =============================================================================
# mcap_populate_irradiance.R
# Reads irradiance_combined.csv, parses date/time, constructs POSIXct
# date_time column, and writes mcap_irrad.csv for the hsat pipeline.
#
# Author : Ji Hoon Han
# Updated: 2026-04-12
#
# Input  (data/input/):
#   irradiance_combined.csv    -- 1-minute Epar time series (135,774 rows)
#                                 Date format: YYYY-MM-DD
#                                 Time format: H:MM:SS (no AM/PM)
#
# Output (data/transformed/):
#   mcap_irrad.csv
# =============================================================================

library(hms)

# -----------------------------------------------------------------------------
# 1. Load
# -----------------------------------------------------------------------------
irradiance <- read.csv("data/input/irradiance_raw.csv",
                       sep = ",", stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 2. Date — convert from M/D/YYYY to YYYY-MM-DD character
# -----------------------------------------------------------------------------
library(lubridate)
irradiance$Date <- format(mdy(irradiance$Date), "%Y-%m-%d")

# -----------------------------------------------------------------------------
# 3. Time — parse as HMS
# -----------------------------------------------------------------------------
irradiance$Time <- as_hms(irradiance$Time)

# -----------------------------------------------------------------------------
# 4. date_time — combine Date + Time into POSIXct (local timezone)
# -----------------------------------------------------------------------------
irradiance$date_time <- as.POSIXct(
  paste(irradiance$Date, irradiance$Time),
  format = "%Y-%m-%d %H:%M:%S"
)

# -----------------------------------------------------------------------------
# 5. Write
# -----------------------------------------------------------------------------
write.csv(irradiance, "data/transformed/mcap_irrad.csv", row.names = FALSE)