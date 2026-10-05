# ==============================================================================
# Script 14: Merge chunked outputs of Script 13 (random-slopes BF calibration)
# ==============================================================================
# Reads every <root>/rows_*/res_rows_*.rds written by run_chunk() in Script 13,
# checks that each expected Design row is present exactly once and that its
# design columns match make_design_rs(), and saves
#   <root>/res_bf_calibration_rs_merged.rds   (summaries, one row per Design row)
#   <root>/res_bf_calibration_rs_raw.rds      (all replications)
#
# Usage:
#   Rscript scripts/14_Merge_BFCalibration_RS.R [root] [expected_rows]
#   root           default output/Simulation_BFCalibration_RS
#                  (relative paths are taken relative to the project root)
#   expected_rows  default "all" (e.g. "1-2" for a test merge)
# Stops with an error if a row is missing or duplicated.
# ==============================================================================

library(SimDesign)
library(here)

source(here("R", "bf_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) > 0) args[1] else file.path("output", "Simulation_BFCalibration_RS")
if (!grepl("^/", root)) root <- here(root)
expected_spec <- if (length(args) > 1) args[2] else "all"

Design   <- make_design_rs()
expected <- parse_rows_spec(expected_spec, nrow(Design))

files <- list.files(root, pattern = "^res_rows_.*\\.rds$", recursive = TRUE,
                    full.names = TRUE)
if (length(files) == 0) stop("No chunk results (res_rows_*.rds) under ", root)
message("Found ", length(files), " chunk file(s) under ", root)

chunks <- lapply(files, readRDS)
summ <- do.call(rbind, lapply(chunks, function(x) as.data.frame(x)))
raw  <- do.call(rbind, lapply(chunks, function(x) as.data.frame(SimResults(x))))

# ---- Completeness checks -----------------------------------------------------
tab      <- table(factor(summ$design_row, levels = seq_len(nrow(Design))))
missing  <- setdiff(expected, summ$design_row)
dupes    <- as.integer(names(tab)[tab > 1])
extra    <- setdiff(unique(summ$design_row), expected)

design_cols <- c("arm", "true_b1", "n_subjects", "n_trials", "true_sd_slope")
ref      <- as.data.frame(Design)[match(summ$design_row, Design$design_row), design_cols]
mismatch <- which(!apply(summ[, design_cols] == ref, 1, all))

message(sprintf("Rows expected: %d | present: %d | missing: %d | duplicated: %d | unexpected: %d | design-column mismatches: %d",
                length(expected), length(unique(summ$design_row)), length(missing),
                length(dupes), length(extra), length(mismatch)))
if (length(missing))  message("Missing rows: ", paste(missing, collapse = ","))
if (length(dupes))    message("Duplicated rows: ", paste(dupes, collapse = ","))
if (length(extra))    message("Unexpected rows: ", paste(extra, collapse = ","))
if (length(missing) || length(dupes) || length(mismatch)) {
  stop("Merge check failed (see messages above); nothing written.")
}

# Replications per row in the raw results
reps_per_row <- table(raw$design_row)
message("Replications per row: min ", min(reps_per_row), ", max ", max(reps_per_row))

summ <- summ[order(summ$design_row), ]
raw  <- raw[order(raw$design_row), ]
saveRDS(summ, file.path(root, "res_bf_calibration_rs_merged.rds"))
saveRDS(raw,  file.path(root, "res_bf_calibration_rs_raw.rds"))
message("MERGE OK: ", nrow(summ), " rows, ", nrow(raw), " replications -> ",
        file.path(root, "res_bf_calibration_rs_merged.rds"))
