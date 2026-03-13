# PriorsBinomial

> How to Set Priors for Hypothesis Testing in Generalized Linear Models: A Three-Step Workflow with an Application to Binomial Models

**Authors:** Gidon T. Frischkorn¹², Joscha Dutli¹, Philipp Musfeld¹³, Klaus Oberauer¹

¹ Department of Psychology, University of Zurich
² Faculty of Behavioral Sciences and Psychology, University of Lucerne
³ Department of Psychology, University of Amsterdam

**OSF:** Simulation output files and supplementary materials are shared at <https://osf.io/yjgkt/overview>

---

## Overview

Choosing default priors for regression coefficients in binomial GLMs is non-trivial because the link function (logit or probit) mediates how a prior on the linear predictor translates to the probability scale. This project uses large-scale prior predictive simulation to answer:

1. **Goal 1 — Prior Predictive Analysis:** What do priors of different families and widths imply about baseline success probabilities and effect sizes ($\delta p = p_1 - p_2$) on the probability scale?
2. **Goal 2 — BF Calibration:** How well does Bayesian model comparison (Bayes factors via bridge sampling) recover the ground-truth model under these priors?

A central result is the **matched prior principle**: a Logistic(0, 1) prior on the logit scale induces a Uniform(0, 1) distribution on the probability scale (and analogously, Normal(0, 1) on the probit scale). This provides a natural, non-informative starting point for intercept priors.

All simulations use **sum-to-zero (±1) contrast coding** so that the intercept represents the grand mean and the slope encodes the half-difference between conditions.

---

## Repository Structure

```
PriorsBinomial/
├── R/                              # Helper functions (sourced by scripts & reports)
│   ├── link_functions.R            # logit/ilogit, probit/iprobit, analytical densities
│   ├── prior_sampling.R            # sample_prior() — draw from normal/logistic/cauchy
│   ├── plotting.R                  # ggplot2 figure functions
│   └── bf_helpers.R                # fit_and_get_bf() — brms + bridgesampling (Goal 2)
│
├── scripts/                        # Analysis pipeline — run in numbered order
│   ├── 01_PriorPredictive.R        # Goal 1: Prior predictive (fixed effects)
│   ├── 01b_PriorPredictive_RE.R    # Goal 1: Prior predictive (random effects)
│   ├── 02_SaveData_PriorPredictive.R       # Goal 1: post-process FE results
│   ├── 02b_SaveData_PriorPredictive_RE.R   # Goal 1: post-process RE results
│   ├── 03_BFCalibration.R          # Goal 2: BF calibration (fixed effects)
│   ├── 03b_BFCalibration_RE.R      # Goal 2: BF calibration (random effects)
│   ├── 04_SaveData_BFCalibration.R         # Goal 2: post-process FE BF results
│   ├── 04b_SaveData_BFCalibration_RE.R     # Goal 2: post-process RE BF results
│   └── 05_BFValidation_BridgeSampling.R    # Validate Savage-Dickey BFs vs bridge sampling
│
├── reports/                        # Quarto documents (apaquarto format)
│   ├── priors_binomial_glm_ampps.qmd           # Main manuscript
│   ├── priors_binomial_glm_ampps_supplement.qmd # Supplementary materials
│   └── references.bib              # BibTeX bibliography
│
├── output/                         # Simulation outputs (gitignored; shared on OSF)
├── figures/                        # Rendered figures (gitignored)
└── local/                          # Working files, notes (gitignored)
```

---

## Simulation Design

### Goal 1 — Prior Predictive (Fixed Effects)

| Factor | Levels |
|--------|--------|
| Link function | `logit`, `probit` |
| Intercept prior family (`dist_b0`) | `logistic`, `normal`, `cauchy` |
| Intercept prior SD (`sd_b0`) | 0.50, 0.75, 1.00, 1.50 |
| Slope prior family (`dist_b1`) | `normal`, `logistic`, `cauchy` |
| Slope prior SD (`sd_b1`) | 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, 2.00 |

**Total conditions:** 2 × 3 × 4 × 3 × 10 = **720**
**Draws per condition:** 10 replications × 10,000 draws = 100,000

**Model:**

$$
g(p_{ij}) = \beta_0 + \beta_1 x_i, \quad x_i \in \{+1, -1\}
$$

$$
\delta p = p_1 - p_2 = g^{-1}(\beta_0 + \beta_1) - g^{-1}(\beta_0 - \beta_1)
$$

### Goal 1b — Prior Predictive (Random Effects)

Extends Goal 1 to mixed models with subject-level random intercepts and/or slopes, sweeping hyperprior widths for the random-effect standard deviations.

### Goal 2 — BF Calibration

For each prior configuration, 500 simulated datasets are generated under the null ($\beta_1 = 0$) and alternative ($\beta_1 \neq 0$). Models are fit with `brms` and Bayes factors are computed via `bridgesampling` to assess false-positive and true-positive rates.

### BF Validation (Script 05)

Validates Savage-Dickey density ratio Bayes factors against bridge sampling on a targeted subset of conditions to confirm equivalence of the two BF computation methods.

---

## Quickstart

### 1. Install required packages

```r
install.packages(c(
  "SimDesign", "ggplot2", "ggridges", "viridis",
  "scales", "dplyr", "tidyr", "here", "knitr"
))

# Goal 2 only
install.packages(c("brms", "bridgesampling"))
```

### 2. Run the pipeline

Scripts must be run in order. Each script guards against re-running if output already exists.

```r
# Goal 1 — fixed-effects prior predictive
source("scripts/01_PriorPredictive.R")
source("scripts/02_SaveData_PriorPredictive.R")

# Goal 1b — random-effects prior predictive
source("scripts/01b_PriorPredictive_RE.R")
source("scripts/02b_SaveData_PriorPredictive_RE.R")

# Goal 2 — BF calibration (long-running; recommended on a compute server)
source("scripts/03_BFCalibration.R")
source("scripts/03b_BFCalibration_RE.R")
source("scripts/04_SaveData_BFCalibration.R")
source("scripts/04b_SaveData_BFCalibration_RE.R")

# BF validation
source("scripts/05_BFValidation_BridgeSampling.R")
```

> **Runtime:** Goal 1 (~5–15 min with parallel workers). Goal 2 is substantially longer due to repeated `brms` model fitting.

### 3. Render the manuscript

```r
quarto::quarto_render("reports/priors_binomial_glm_ampps.qmd")
quarto::quarto_render("reports/priors_binomial_glm_ampps_supplement.qmd")
```

Or from the terminal:

```bash
quarto render reports/priors_binomial_glm_ampps.qmd
quarto render reports/priors_binomial_glm_ampps_supplement.qmd
```

The manuscript renders to HTML, PDF, and DOCX via the apaquarto format.

---

## Key Design Decisions

- **Matched prior principle:** Logistic(0, 1) on the logit scale → Uniform(0, 1) on the probability scale (exact by the probability integral transform). Normal(0, 1) on the probit scale yields the same result. These serve as canonical non-informative intercept priors.
- **Contrast coding:** All simulations use sum-to-zero ±1 coding. The intercept ($\beta_0$) represents the grand-mean log-odds; $\beta_1$ encodes the half-difference. This gives symmetric prior predictive distributions for $\delta p$.
- **Intercept vs. slope priors:** Both matched and mismatched prior families are simulated (e.g., Normal prior on the logit scale) to quantify the cost of prior misspecification.
- **Parallel execution:** Scripts use `parallel::detectCores() - 2` workers via SimDesign's `ncores` argument.

---

## Recommended Priors (Summary)

Based on the prior predictive analysis:

| Parameter | Link | Recommended Prior | Notes |
|-----------|------|-------------------|-------|
| Intercept $\beta_0$ | logit | Logistic(0, 1) | Uniform(0,1) on probability scale |
| Intercept $\beta_0$ | probit | Normal(0, 1) | Uniform(0,1) on probability scale |
| Slope $\beta_1$ | logit / probit | Normal(0, σ), σ ∈ [0.10, 0.50] | Encodes small-to-moderate effects; avoid Cauchy |

Ready-to-use `brms` prior specifications are provided in the manuscript's Recommendations section.

---

## Data Availability

Simulation output files are shared on OSF: <https://osf.io/yjgkt/overview>

The `output/` directory is gitignored due to file size. To reproduce locally, run the pipeline scripts or download the output files from OSF.

---

## Citation

If you use this simulation framework, please cite:

- Frischkorn, G. T., Dutli, J., Musfeld, P., & Oberauer, K. (in preparation). *How to set priors for hypothesis testing in generalized linear models: A three-step workflow with an application to binomial models.*
- Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software*, *80*(1), 1–28.

---

## License

This project is for research purposes. See individual file headers for function-level documentation.
