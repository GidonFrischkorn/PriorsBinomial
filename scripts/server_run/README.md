# Random-slopes BF calibration on a remote VM (UZH ScienceCloud, Ubuntu)

Runs `scripts/13_BFCalibration_RS.R` in mode `full`: the approved grid
(grid C + logit slope SD 0.50), 72 Design rows, 200 replications, bridge
sampling with 4 x 5,000 post-warmup draws. Written 2026-10-05. Nothing in this
file has been run on a VM yet.

## 0. Before you start (on the laptop)

- **Commit and push first.** `scripts/13_BFCalibration_RS.R`,
  `scripts/14_Merge_BFCalibration_RS.R`, `R/bf_helpers.R` and
  `R/link_functions.R` have uncommitted changes, and `git clone` does not
  include uncommitted work.
- This folder (`scripts/server_run/`) is tracked, so the clone includes it.
- `output/` is gitignored, so nothing large has to be copied to the VM.

## 1. VM setup (Ubuntu; not tested on a VM)

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran git rsync \
  libcurl4-openssl-dev libssl-dev libxml2-dev libv8-dev
# R 4.6.x from the CRAN Ubuntu repository (https://cran.r-project.org/bin/linux/ubuntu/)
git clone -b revision-psychmet-r1 <repo-url> ~/PriorsBinomial
```

R packages: these versions come from `sessionInfo()` of the bridge smoke run
(`logs/t3b_smoke_bridge.log`, macOS, R 4.6.1). Install the same or newer:

| Package | Version | Package | Version |
|---|---|---|---|
| SimDesign | 2.27 | brms | 2.23.0 |
| rstan | 2.32.7 | StanHeaders | 2.39.1 |
| cmdstanr | 0.9.0 (from stan-dev r-universe) | posterior | 1.7.0 |
| bridgesampling | 1.2-1 | logspline | 2.1.22 |
| withr | 3.0.3 | here | 1.0.2 |
| qs2 (SimDesign temp files) | 0.3.1 | Rcpp | 1.1.2 |
| dplyr | 1.2.1 | tibble | 3.3.1 |

CmdStan **2.40.0**.

```r
install.packages(c("SimDesign", "brms", "rstan", "posterior", "bridgesampling",
                   "logspline", "withr", "here", "qs2"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan(version = "2.40.0", cores = 4)
```

If CmdStan is installed somewhere non-standard, set `export CMDSTAN=/path/to/cmdstan-2.40.0`.

## 2. Pre-flight on the VM (about 25 min)

```bash
cd ~/PriorsBinomial && mkdir -p logs
timeout 2h Rscript scripts/13_BFCalibration_RS.R smoke_bridge > logs/vm_smoke_bridge.log 2>&1
timeout 2h Rscript scripts/13_BFCalibration_RS.R chunkcheck   > logs/vm_chunkcheck.log 2>&1
grep -E "SMOKE_BRIDGE|CHUNKCHECK|EXIT|Error" logs/vm_*.log
```

Compare the `fit_wall mean` values with this laptop's:

| Row | Laptop `fit_wall mean` |
|---|---|
| probit ns50/nt20 | 11.6 s |
| logit_continuous ns50/nt20 | 104.5 s |
| logit ns50/nt50 | 12.0 s |

Scale the sizing table below by the ratio. `chunkcheck` must print
`identical data: TRUE | identical BFs: TRUE`.

## 3. Launch (resumable, safe to re-launch)

```bash
cd ~/PriorsBinomial
PB_NCORES=30 bash scripts/server_run/run_rs_local.sh          # all 72 rows, log: logs/t5_full.log
# or split across two VMs (balanced by estimated CPU-h, see below):
PB_NCORES=30 bash scripts/server_run/run_rs_local.sh "1,3,4,..."   # log: logs/t5_full_rows_<rows>.log
tail -f logs/t5_full.log      # done when a line "==== EXIT 0" appears
```

The launcher is detached (`nohup`), so it survives an SSH disconnect.

**Re-launching:**
- Running the same command again is safe. A pid file stops a duplicate run on the same rows.
- After a crash or kill, SimDesign resumes from `SIMDESIGN-TEMPFILE_rows_<label>` in the chunk folder, and finished Design rows are skipped. This was tested on the laptop: the run was killed after row 1 of 2 and re-launched; the log shows "Resuming simulation from ..." and the run then finished (`logs/t3b_resume_launcher.log`).
- A finished chunk, i.e. one where `res_rows_<label>.rds` exists, is never re-run. Templates still compile first, which takes about 2.5 min.
- To redo a chunk from scratch, delete its folder `output/Simulation_BFCalibration_RS/rows_<label>/`.

**Environment variables:**

| Variable | Default | Meaning |
|---|---|---|
| `PB_NCORES` | `detectCores()-2` | Workers; one core each, chains run sequentially. |
| `PB_REPS` | 200 | Replications per Design row. |
| `PB_SEED` | 20261005 | Global seed. |
| `PB_TIMEOUT` | 21d | Hard timeout for the run. |
| `PB_OUT_ROOT` | `output/Simulation_BFCalibration_RS` | Output root folder. |

Per-row seeds are derived from the global seed for the full 72-row Design, so
results do not depend on how the rows are split across chunks or VMs. This was
verified: rows 1-2 run as one chunk and as two one-row chunks gave identical
data and all 24 identical BFs (`logs/t3b_chunkcheck.log`).

## 4. Copy results back and merge

```bash
# on the laptop
rsync -avz <user>@<vm>:~/PriorsBinomial/output/Simulation_BFCalibration_RS ./output/
rsync -avz <user>@<vm>:~/PriorsBinomial/logs/t5_full* ./logs/
Rscript scripts/14_Merge_BFCalibration_RS.R      # checks all 72 rows present exactly once
```

With two VMs, rsync both into the same `output/Simulation_BFCalibration_RS/`.
The chunk folders have distinct names. The merge writes
`res_bf_calibration_rs_merged.rds` (summaries) and
`res_bf_calibration_rs_raw.rds` (all replications). It stops if any row is
missing or duplicated.

## 5. Sizing (EXTRAPOLATION from laptop timings, `logs/t3b_sizing.log`)

The CPU-hours come from per-fit wall times measured on this laptop (Apple
Silicon, 1 core per worker).
- **Measured directly:** logit ns30/nt50 and ns50/nt50, probit ns50/nt20, continuous ns50/nt20. These are 24 of 72 rows and 46% of the CPU-h.
- **Derived for the other rows:** the cheaper Savage-Dickey timing of that row's cell × the measured bridge/Savage-Dickey ratio.
- **Spot check:** row 1 of the full grid (logit ns20/nt20) took 151 s per replication in the resume test, against 155 s predicted.

Totals:
- 1,027 CPU-h for 200 replications: logit_factor 582, logit_continuous 411, probit_factor 34.
- The largest single row is design_row 39 (continuous, ns = 50), at about 76 CPU-h.

| VM | Workers | Ideal wall-clock | With last-batch idling* |
|---|---|---|---|
| 16 vCPU | 14 | 73 h (3.1 d) | ~77 h (3.2 d) |
| 32 vCPU | 30 | 34 h (1.4 d) | ~36 h (1.5 d) |
| 64 vCPU | 62 | 17 h (0.7 d) | ~21 h (0.9 d) |

\*SimDesign parallelises over the 200 replications within a row. With w
workers, the last of ceil(200/w) batches is partly idle; the adjustment assumes
equal replication times.

**VM speed is unknown.** A vCPU is often a hyperthread and can be markedly
slower than a laptop performance core. Rescale using the pre-flight smoke run.
RAM was not measured on the workers; the SimDesign main process reported about
285 MB.

**Balanced 2-VM split** (about 515 / 512 CPU-h):

```
VM1 PB_ROWS="1,3,4,6,8,12,14,16,18,20,21,22,24,26,29,31,32,37,39,41,43,45,47,49,51,53,55,57,58,60,62,64,66,68,70,72"
VM2 PB_ROWS="2,5,7,9,10,11,13,15,17,19,23,25,27,28,30,33,34,35,36,38,40,42,44,46,48,50,52,54,56,59,61,63,65,67,69,71"
```

Contiguous halves `1-36` / `37-72` are unbalanced (455 / 572 CPU-h).

## 6. Consistency re-run with bridge sampling (after T5; Script 16)

`scripts/16_BFCalibration_Consistency_Bridge.R` repeats the consistency
simulation of Script 03c (between-subjects, Logistic(0, 0.25) effect prior,
Logistic(0, 0.75) intercept, random intercepts with Exponential(4), 6 true
effects x 5 sample sizes, n_trials = 20) with the bridge engine of
`R/bf_helpers.R`: two rstan models compiled once (H1, H0), refits by
`update()`, 4 x 5,000 post-warmup draws, 50 replications per cell (3,000
brms fits). It is independent of the random-slopes run and needs only the
packages of section 1.

```bash
cd ~/PriorsBinomial && mkdir -p logs
# pre-flight: one cell, one replication (about 5 min; expects a line SMOKE_CONSISTENCY_BRIDGE)
timeout 30m Rscript scripts/16_BFCalibration_Consistency_Bridge.R smoke > logs/t7_consistency_smoke_vm.log 2>&1; echo "EXIT $?" >> logs/t7_consistency_smoke_vm.log
# full run, detached
PB_NCORES=30 nohup bash -c 'timeout 7d Rscript scripts/16_BFCalibration_Consistency_Bridge.R full; echo "==== EXIT $?"' > logs/t7_consistency_full.log 2>&1 < /dev/null &
tail -f logs/t7_consistency_full.log      # done when "CONSISTENCY BRIDGE COMPLETE" and "==== EXIT 0" appear
```

Environment variables: `PB_NCORES` (workers), `PB_REPS` (default 50),
`PB_SEED` (default 20261104; per-cell seeds via `genSeeds()` on the Design),
`PB_N_SAMPLING` (post-warmup draws per chain, default 5000).

Re-launching after a crash is safe: SimDesign resumes from
`output/Simulation_BFCalibration_Consistency_Bridge/SIMDESIGN-TEMPFILE_consistency_bridge`,
and a finished run (`res_bf_calibration_consistency_bridge.rds` present) is
skipped; the flat files are then rebuilt from it. Copy back and inspect:

```bash
# on the laptop
rsync -avz <user>@<vm>:~/PriorsBinomial/output/Simulation_BFCalibration_Consistency_Bridge ./output/
rsync -avz <user>@<vm>:~/PriorsBinomial/output/bf_calibration_consistency_bridge_*.rds ./output/
rsync -avz <user>@<vm>:~/PriorsBinomial/logs/t7_consistency_* ./logs/
```

Outputs: `output/bf_calibration_consistency_bridge_raw.rds` (one row per
replication; `bf10_bridge`, `log10_bf10`, diagnostics, `n_fail`) and
`output/bf_calibration_consistency_bridge_summaries.rds` (one row per cell).
Failed fits stay NA (`allow_na = TRUE`) and are counted in `n_fail`.
