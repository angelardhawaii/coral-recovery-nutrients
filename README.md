# Coral Recovery Under Nutrient Enrichment

Analysis scripts accompanying the publication:

> **Nutrient enrichment enhances post-bleaching survival with asynchronous symbiont recovery**
— Ji Hoon Han, Angela Richards Donà, Claire Lewis, Rob Toonan, Celia Smith —
*Proceedings of the Royal Society B*
2026. DOI: [DOI]

---

## Overview

This repository contains the R analysis scripts used to evaluate the effects of
nutrient enrichment on coral recovery following a thermal stress event. Corals
were exposed to three nutrient treatments — control (CTL), ammonium (NH4), and
nitrate (NO3) — and monitored across three experimental phases: acclimation,
heating, and recovery.

---

## Scripts

| File | Description |
|------|-------------|
| `growth.R` | Buoyant weight growth rates across acclimation, heating, and recovery phases; linear mixed models with post-hoc comparisons |
| `survivorship.R` | Kaplan-Meier survival curves and Cox proportional hazards models by treatment |
| `photosynthesis_endpoint.R` | Endpoint photosynthetic metrics (F₀, Fₘ, quantum yield, NPQ, Pmax, Ek, Hsat, DSPI) during the recovery phase; linear mixed models |
| `photosynthesis_trend(GAMM).R` | Generalized additive mixed models (GAMMs) of photosynthetic metric trajectories over time during recovery; derivative analysis of rate of change |

---

## Data

Input data files are located in `data/input/`:

| File | Description |
|------|-------------|
| `surv_bw.csv` | Buoyant weight measurements and survival records by individual, treatment, and experimental phase |
| `rlc_data.csv` | Rapid light curve (RLC) photosynthetic parameters by individual, treatment, genotype, and time point |

---

## Dependencies

All scripts are written in R. Required packages:

```
dplyr, tidyr, ggplot2, lme4, lmerTest, emmeans, multcomp, multcompView,
survival, survminer, mgcv, gamm4, gratia, cowplot, gridExtra, performance
```

---

## Contact

**Ji Hoon (Justin) Han** — han9@hawaii.edu
Or
**Angela Richards Donà** — angelard@hawaii.edu
