# PriorsBinomial

> Principled Default Priors for Bayesian Hypothesis Testing in Binomial GLMs: A Matched-Prior Framework with Empirical Evaluation

**Authors:** Gidon T. Frischkorn¹, Joscha Dutli¹, Philipp Musfeld²¹, Klaus Oberauer¹

¹ Department of Psychology, University of Zurich
² Department of Psychology, University of Amsterdam

**Status:** Submitted to *Psychological Methods* (July 2026). An earlier, tutorial-oriented version is available as a [preprint](https://doi.org/10.31234/osf.io/q7byw_v1).

**OSF:** Pre-computed simulation outputs are archived at <https://doi.org/10.17605/OSF.IO/YJGKT>

---

## Overview

Choosing default priors for regression coefficients in binomial GLMs is non-trivial because the link function (logit or probit) mediates how a prior on the linear predictor translates to the probability scale. The paper proposes and validates a three-step workflow:

1. **Step 1 — Matched priors:** Identify prior families whose CDF matches the inverse link function. By the probability integral transform, a Logistic(0, 1) prior on the logit scale induces a Uniform(0, 1) distribution on the probability scale (analogously, Normal(0, 1) for probit) — a maximally uninformative anchor.
2. **Step 2 — Prior predictive simulation:** Evaluate what priors of different families and scales imply for baseline success probabilities and effect sizes ($\delta p = p_1 - p_2$) on the probability scale.
3. **Step 3 — Bayes factor calibration:** Verify via simulation that the recommended priors yield well-calibrated Bayes factors (Savage-Dickey density ratio, validated against bridge sampling), including a direct comparison to the JZS/Cauchy default and formal consistency checks.

A central result of the calibration analysis is that prior scale acts as a **decision criterion, not a sensitivity parameter**: it shifts the Bayes factor's operating point on the ROC curve but leaves discriminability unchanged. The advantage of matched priors is principled scale selection on the probability scale, where effect sizes are directly interpretable.

All simulations use **sum-to-zero (±1) contrast coding** so that the intercept represents the grand mean and the slope encodes the half-difference between conditions on the link scale.

---

## Repository Structure

```text
PriorsBinomial/
├── R/                              # Helper functions (sourced by scripts & reports)
│   ├── link_functions.R            # logit/ilogit, probit/iprobit, analytical densities
│   ├── prior_sampling.R            # sample_prior() — draw from normal/logistic/cauchy
│   ├── plotting.R                  # ggplot2 figure functions
│   └── bf_helpers.R                # fit_and_get_bf() — brms + bridgesampling
│
├── scripts/                        # Analysis pipeline — run in numbered order
│   ├── 01_PriorPredictive.R        # Step 2: Prior predictive (fixed effects)
│   ├── 01b_PriorPredictive_RE.R    # Step 2: Prior predictive (random-effects hyperpriors)
│   ├── 02_SaveData_PriorPredictive.R       # Post-process FE results
│   ├── 02b_SaveData_PriorPredictive_RE.R   # Post-process RE results
│   ├── 03_BFCalibration.R          # Step 3: BF calibration (fixed effects)
│   ├── 03b_BFCalibration_RE.R      # Step 3: BF calibration (RE hyperpriors)
│   ├── 03c_BFCalibration_Consistency.R     # Step 3: Monotone-growth & Lindley consistency
│   ├── 03d_BFCalibration_JZS.R     # Step 3: JZS/Cauchy comparison (scale × family)
│   ├── 03e_BFCalibration_Probit.R  # Step 3: Probit confirmation (backup simulation)
│   ├── 04_SaveData_BFCalibration.R         # Post-process FE BF results
│   ├── 04b_SaveData_BFCalibration_RE.R     # Post-process RE BF results
│   ├── 04c_SaveData_BFCalibration_Consistency.R
│   ├── 04d_SaveData_BFCalibration_JZS.R
│   ├── 05_BFValidation_BridgeSampling.R    # Validate Savage-Dickey BFs vs bridge sampling
│   ├── 06_BFValidation_Figures.R   # Figures for the BF validation
│   ├── 07_AppliedExample.R         # Applied example (Oberauer, 2019, simple span)
│   ├── 08_Figure_JZS_ROC.R         # Manuscript Figure: dual ROC comparison
│   ├── 09_Figure_JZS_AUC.R         # Manuscript Figure: AUC comparison
│   └── 10_Figure_JZS_logBF.R       # Supplement Figure: logBF distributions
│
├── reports/                        # Quarto documents (apaquarto format)
│   ├── priors_binomial_glm_psychmet.qmd                       # Main manuscript
│   ├── priors_binomial_glm_psychmet_supplement.qmd            # Supplementary materials
│   ├── priors_binomial_glm_psychmet_translational_abstract.qmd # Translational abstract
│   └── references.bib              # BibTeX bibliography
│
├── output/                         # Simulation outputs (gitignored; archived on OSF)
├── figures/                        # Rendered figures (gitignored; regenerate via scripts 06, 08–10)
└── local/                          # Working files, notes (gitignored)
```

---

## Simulation Design

### Step 2 — Prior Predictive (Fixed Effects)

| Factor | Levels |
| -------- | -------- |
| Link function | `logit`, `probit` |
| Intercept prior family (`dist_b0`) | `logistic`, `normal`, `cauchy` |
| Intercept prior SD (`sd_b0`) | 0.50, 0.75, 1.00, 1.50 |
| Effect prior family (`dist_b1`) | `normal`, `logistic`, `cauchy` |
| Effect prior SD (`sd_b1`) | 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, 2.00 |

**Total conditions:** 2 × 3 × 4 × 3 × 10 = **720**, with 10 replications × 10,000 draws each.

**Model:**

$$
g(p_{ij}) = \beta_0 + \beta_1 x_i, \quad x_i \in \{+1, -1\}
$$

$$
\delta p = p_1 - p_2 = g^{-1}(\beta_0 + \beta_1) - g^{-1}(\beta_0 - \beta_1)
$$

### Step 2 — Prior Predictive (Random-Effects Hyperpriors)

Crosses four hyperprior families (Exponential, Gamma(shape = 2), Half-Normal, Half-Student-*t*(df = 3)) with five hyperprior means (0.15, 0.25, 0.50, 1.00, 2.50) and two random-effects structures (intercept only; intercept + slope), holding the fixed-effects priors at the matched configuration — 40 conditions.

### Step 3 — BF Calibration

Data are simulated under the logit link with true effects $b_1 \in \{0, 0.05, 0.10, 0.20, 0.50\}$, crossing sample size (30, 60, 100), trials per subject (20, 50), effect prior family (Normal, Logistic, Cauchy) at scales 0.25 and 0.50, and matched vs. misfit intercept priors — 360 conditions with 100 replications each. Bayes factors are computed via the Savage-Dickey density ratio; a separate simulation varies only the random-effects hyperprior.

Additional Step-3 simulations:

- **Consistency (03c):** Monotone evidence growth with effect size and sample-size (Lindley) consistency for the matched Logistic(0, 0.25) prior (bridge sampling).
- **JZS/Cauchy comparison (03d):** Two prior scales (0.25 vs. √2/2 ≈ 0.707) as design factors; within every replication, both families (Logistic, Cauchy) are fitted to the same datasets — 36 conditions × 300 replications. Basis of the ROC/AUC analysis in the manuscript.
- **Probit confirmation (03e):** Backup simulation confirming well-behaved BF calibration for the matched Normal prior under the probit link.

### BF Validation (Scripts 05–06)

Validates Savage-Dickey density ratio Bayes factors against bridge sampling on a targeted subset of conditions.

---

## Quickstart

### 1. Install required packages

```r
install.packages(c(
  "SimDesign", "ggplot2", "ggridges", "viridis",
  "scales", "dplyr", "tidyr", "here", "knitr"
))

# BF calibration (Step 3) additionally requires
install.packages(c("brms", "bridgesampling"))
```

### 2. Run the pipeline (or download outputs from OSF)

The compute-heavy simulation outputs are archived on OSF (<https://doi.org/10.17605/OSF.IO/YJGKT>). To reproduce the manuscript without re-running the simulations, download the `output/` folder from OSF into the repository root and skip to step 3.

Scripts must be run in order; each guards against re-running if output already exists.

```r
# Step 2 — prior predictive
source("scripts/01_PriorPredictive.R")
source("scripts/01b_PriorPredictive_RE.R")
source("scripts/02_SaveData_PriorPredictive.R")
source("scripts/02b_SaveData_PriorPredictive_RE.R")

# Step 3 — BF calibration (long-running; recommended on a compute server)
source("scripts/03_BFCalibration.R")
source("scripts/03b_BFCalibration_RE.R")
source("scripts/03c_BFCalibration_Consistency.R")
source("scripts/03d_BFCalibration_JZS.R")
source("scripts/03e_BFCalibration_Probit.R")
source("scripts/04_SaveData_BFCalibration.R")
source("scripts/04b_SaveData_BFCalibration_RE.R")
source("scripts/04c_SaveData_BFCalibration_Consistency.R")
source("scripts/04d_SaveData_BFCalibration_JZS.R")

# BF validation
source("scripts/05_BFValidation_BridgeSampling.R")
```

### 3. Generate the pre-computed figures

The manuscript and supplement embed three figures generated outside the Quarto render:

```r
source("scripts/06_BFValidation_Figures.R")
source("scripts/08_Figure_JZS_ROC.R")
source("scripts/09_Figure_JZS_AUC.R")
source("scripts/10_Figure_JZS_logBF.R")
```

### 4. Render the manuscript

```bash
quarto render reports/priors_binomial_glm_psychmet.qmd
quarto render reports/priors_binomial_glm_psychmet_supplement.qmd
```

The documents render to HTML, PDF, and DOCX via the apaquarto format.

---

## Key Design Decisions

- **Matched prior principle:** Logistic(0, 1) on the logit scale → Uniform(0, 1) on the probability scale (exact by the probability integral transform); Normal(0, 1) on the probit scale yields the same result. Scaling below 1 concentrates mass away from floor and ceiling.
- **Contrast coding:** All simulations use sum-to-zero ±1 coding. The intercept ($\beta_0$) represents the grand-mean on the link scale; $\beta_1$ encodes the half-difference. This gives symmetric prior predictive distributions for $\delta p$ and makes the BF invariant to group labeling.
- **Intercept vs. effect priors:** Both matched and misfit prior families are simulated to quantify the consequences of prior-family misspecification.
- **BF estimation:** Savage-Dickey density ratio for the calibration grids (validated against bridge sampling in scripts 05–06); bridge sampling for the consistency analysis.
- **Parallel execution:** Scripts use `parallel::detectCores() - 2` workers via SimDesign's `ncores` argument.

---

## Recommended Priors (Summary)

Default recommendations from the manuscript (sum-to-zero ±1 coding, predictors scaled to SD = 1):

| Parameter | Link | Recommended prior | Key property |
| ----------- | ------ | ------------------- | -------------- |
| Intercept $b_0$ | logit | Logistic(0, 0.75) | < 10% prior mass in floor/ceiling region |
| Intercept $b_0$ | probit | Normal(0, 0.75) | same |
| Effect $b_1$ | logit | Logistic(0, s), s = 0.25–0.30 | 90th percentile of $\lvert\delta p\rvert$ ≈ 0.29–0.34 |
| Effect $b_1$ | probit | Normal(0, s), s = 0.25–0.30 | 90th percentile of $\lvert\delta p\rvert$ ≈ 0.27–0.32 |
| RE hyperprior $\sigma_u$ | both | Exponential(4) | mean = 0.25 link-scale units |

Rescale the effect prior when the predictor SD differs from 1 (see the manuscript's *From Defaults to Informed Priors* section for scaling formulas and orthonormal-contrast lookup tables), and report BFs across $\mathrm{sd}_{b_1} \in \{0.15, 0.25, 0.35\}$ as a sensitivity check. Ready-to-use `brms` prior specifications are provided in the manuscript's Recommendations section.

---

## Data Availability

Pre-computed simulation outputs are archived on OSF under a CC-BY 4.0 license: <https://doi.org/10.17605/OSF.IO/YJGKT>

The `output/` and `figures/` directories are gitignored due to file size. To reproduce the manuscript locally, download `output/` from OSF (or re-run the pipeline), regenerate `figures/` via scripts 06 and 08–10, and render the Quarto documents.

---

## Citation

- Frischkorn, G. T., Dutli, J., Musfeld, P., & Oberauer, K. (2026). *Principled default priors for Bayesian hypothesis testing in binomial GLMs: A matched-prior framework with empirical evaluation.* Manuscript submitted for publication.
- An earlier, tutorial-oriented version is available as a preprint: Frischkorn, G. T., Dutli, J., Musfeld, P., & Oberauer, K. (2026). *How to set priors for hypothesis testing in generalized linear models: A three-step workflow with an application to binomial models.* PsyArXiv. <https://doi.org/10.31234/osf.io/q7byw_v1>

---

## License

Simulation outputs and materials are archived on OSF under CC-BY 4.0. See individual file headers in `R/` and `scripts/` for function-level documentation.
