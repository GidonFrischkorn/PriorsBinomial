# ==============================================================================
# Script 15b: Figures for the random-slopes BF calibration (after Script 15)
# ==============================================================================
# Reads the Script 15 outputs and writes PNGs to figures/ (gitignored), one set
# per arm (arms are never pooled):
#   figures/rs_violin_<arm>.png          log10 BF10 by family x scale, facets
#                                        n_subjects (rows) x true_b1 (columns),
#                                        n_trials = 20 and slope SD = 0.30
#                                        (the cells all arms share)
#   figures/rs_violin_<arm>_all.png      logit_factor only: all n_trials x
#                                        slope SD cells (rows = n_trials/slope SD)
#   figures/rs_roc_<arm>_b1_<x>.png      dual ROC (H1 / H0 detection) at true
#                                        effect x, columns = n_subjects, BF = 3
#                                        operating points
#   figures/rs_auc_<arm>.png             AUC vs n_subjects
#
# Input suffix: "" for the merged run, "_smoke" when Script 15 ran on the
# smoke object (auto-detected), or "_" + PB_RS_TAG if that variable is set.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
})
source(here("R", "plotting.R"))

tag <- if (nzchar(Sys.getenv("PB_RS_TAG"))) paste0("_", Sys.getenv("PB_RS_TAG")) else
  if (file.exists(here("output", "bf_calibration_rs_rates.rds"))) "" else "_smoke"
load_out <- function(name) {
  f <- here("output", paste0(name, tag, ".rds"))
  if (!file.exists(f)) stop("Missing ", f, "; run scripts/15_SaveData_BFCalibration_RS.R first")
  readRDS(f)
}
reps   <- load_out("bf_calibration_rs_sd_decomp_reps")   # replication x prior, log10_bf10
curves <- load_out("bf_calibration_rs_roc_curves")
auc    <- load_out("bf_calibration_rs_roc_auc")
message("Figure input suffix: '", tag, "' (source: ", attr(reps, "input_source"), ")")

dir.create(here("figures"), showWarnings = FALSE)
fig_tag <- tag   # figures carry the same suffix so smoke figures are recognisable
save_png <- function(p, name, width, height) {
  f <- here("figures", paste0(name, fig_tag, ".png"))
  ggsave(f, p, width = width, height = height, dpi = 150, bg = "white")
  message("Saved ", f)
}

arms <- unique(reps$arm)
for (a in arms) {
  r <- filter(reps, arm == a)

  # Violin: the shared cells (n_trials = 20, slope SD = 0.30); the logit_factor
  # arm has n_trials = 50 and other slope SDs in addition
  shared <- filter(r, n_trials == 20, true_sd_slope == 0.30)
  if (nrow(shared) == 0) shared <- r   # smoke: logit_factor only has n_trials = 50
  save_png(plot_rs_violin(shared, facet_rows = "n_subjects") +
             labs(title = a, subtitle = "n_trials = 20, true slope SD = 0.30"),
           paste0("rs_violin_", a), width = 11, height = 8)
  if (a == "logit_factor" && n_distinct(paste(r$n_trials, r$true_sd_slope)) > 1) {
    r_all <- mutate(r, cell_label = paste0("T = ", n_trials, ", slope SD = ", true_sd_slope))
    save_png(plot_rs_violin(r_all, facet_rows = "cell_label") + labs(title = a),
             paste0("rs_violin_", a, "_all"), width = 11, height = 12)
  }

  # ROC per true effect; columns = n_subjects; shared cells only
  cv <- filter(curves, arm == a, n_trials == 20 | a == "logit_factor" & n_trials == 50,
               true_sd_slope == 0.30)
  for (b1 in sort(unique(cv$true_b1))) {
    save_png(plot_rs_roc(filter(cv, true_b1 == b1), op_threshold = 3,
                         facet_cols = "n_subjects") +
               labs(title = sprintf("%s, true b1 = %s", a, b1)),
             sprintf("rs_roc_%s_b1_%s", a, b1), width = 12, height = 8)
  }

  # AUC
  au <- filter(auc, arm == a, !is.na(auc))
  if (nrow(au) > 0) {
    save_png(plot_rs_auc(au) + labs(title = a), paste0("rs_auc_", a), width = 10, height = 7)
  } else {
    message("No AUC rows for arm ", a, " (no cell with a matching null cell); figure skipped")
  }
}
message("DONE: script 15b")
