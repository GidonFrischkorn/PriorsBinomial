# ==============================================================================
# Script 06: Figures for BF Validation (Savage-Dickey vs Bridge Sampling)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)

source(here("R", "plotting.R"))

# --- Load replication-level data from per-condition files ---------------------

results_dir  <- here("output", "Simulation_BFValidation")
n_conditions <- 18

raw_rows <- vector("list", n_conditions)
for (i in seq_len(n_conditions)) {
  cond_file <- file.path(results_dir, paste0("BFValid_Cond-", i, ".rds"))
  if (!file.exists(cond_file)) {
    warning("Missing: ", cond_file)
    next
  }
  cond_data    <- readRDS(cond_file)
  raw_rows[[i]] <- cbind(
    as.data.frame(as.list(cond_data$condition)),
    as.data.frame(cond_data$results)
  )
}

bf_raw <- bind_rows(raw_rows)
message("Loaded ", nrow(bf_raw), " replications from ",
        sum(!vapply(raw_rows, is.null, logical(1))), " / ", n_conditions,
        " conditions")

# --- Prepare plotting data ----------------------------------------------------

bf_plot <- bf_raw |>
  filter(!is.na(BF10_SD), !is.na(BF10_BS), BF10_SD > 0, BF10_BS > 0) |>
  mutate(
    log10_SD   = log10(BF10_SD),
    log10_BS   = log10(BF10_BS),
    mean_log10 = (log10_SD + log10_BS) / 2,
    diff_log10 = log10_SD - log10_BS,
    dist_b1    = factor(dist_b1,
                        levels = c("cauchy", "normal", "logistic"),
                        labels = c("Cauchy", "Normal", "Logistic")),
    true_b1_label = factor(paste0("b[1] == ", true_b1)),
    sd_b1_label   = factor(paste0("sd[b[1]] == ", sd_b1)),
    sd_b1_factor  = factor(sd_b1)
  )

message("Valid replications for plotting: ", nrow(bf_plot))

# --- Figure S1: Scatter plot --------------------------------------------------

dir.create(here("figures"), showWarnings = FALSE)

fig_scatter <- plot_bf_validation_scatter(bf_plot)
ggsave(here("figures", "bf_validation_scatter.pdf"), fig_scatter,
       width = 12, height = 8)
message("Saved: figures/bf_validation_scatter.pdf")

# --- Figure S2: Bland-Altman plot ---------------------------------------------

fig_ba <- plot_bf_validation_bland_altman(bf_plot)
ggsave(here("figures", "bf_validation_bland_altman.pdf"), fig_ba,
       width = 12, height = 5)
message("Saved: figures/bf_validation_bland_altman.pdf")
