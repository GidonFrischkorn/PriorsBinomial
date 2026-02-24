# PriorsBinomial

Simulation project establishing principled default priors for Bayes Factor
estimation in binomial GLMs (logistic and probit regression).

## Repository Structure

```
PriorsBinomial/
├── R/                        # Helper functions (source at top of scripts)
│   ├── link_functions.R      # logit/ilogit, probit/iprobit, analytical densities
│   ├── prior_sampling.R      # sample_prior() — draw from normal/logistic/cauchy
│   ├── plotting.R            # ggplot2 figure functions
│   └── bf_helpers.R          # fit_and_get_bf() for Goal 2 (brms + bridgesampling)
├── scripts/                  # Numbered pipeline (run in order)
│   ├── 01_PriorPredictive.R          # Goal 1: SimDesign prior predictive simulation
│   ├── 02_SaveData_PriorPredictive.R # Goal 1: post-process results
│   ├── 03_BFCalibration.R            # Goal 2: SimDesign BF calibration
│   └── 04_SaveData_BFCalibration.R   # Goal 2: post-process results
├── output/                   # gitignored — per-condition .rds from SimDesign
├── figures/                  # gitignored — rendered figures
├── reports/                  # Quarto documents
│   ├── prior_predictive_analysis.qmd
│   └── bf_calibration_analysis.qmd
└── local/                    # gitignored — working files (DOCX, notes)
```

## Key Design Decisions

- **Contrast coding:** All simulations use sum-to-zero ±1 coding. This gives
  symmetric prior predictive distributions for delta_p.
- **Intercept prior:** Both matched (Logistic for logit, Normal for probit) and
  misfit (e.g., Normal for logit) are simulated to quantify the difference.
- **Matched prior principle:** Logistic(0,1) on logit scale → Uniform(0,1) on
  probability scale (exact analytical result by probability integral transform);
  same for Normal(0,1) on probit scale.
- **SimDesign:** Goal 1 uses 10 replications × 10,000 draws per condition.
  Goal 2 uses 500 replications per condition.
- **Reference project:** See ComputationalValidity repo for SimDesign patterns.

## Coding Conventions

- snake_case throughout
- Roxygen2 `@param`/`@return` doc in all `R/` files
- `Attach(condition)` inside Generate/Analyse (SimDesign convention)
- Extra config via `fixed_objects = list(...)` in runSimulation calls
- `here::here()` for all file paths
- Guard pattern: `if (!file.exists(...)) { runSimulation(...) }`
- `parallel::detectCores() - 2` for ncores

## Packages Required

- SimDesign
- ggplot2, ggridges, viridis
- here
- brms, bridgesampling (Goal 2 only)
