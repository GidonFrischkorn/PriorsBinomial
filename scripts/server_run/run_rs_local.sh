#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Launch the random-slopes BF calibration (scripts/13_BFCalibration_RS.R, mode
# "full") on a single machine, detached from the terminal.
#
#   bash scripts/server_run/run_rs_local.sh            # all 72 Design rows
#   bash scripts/server_run/run_rs_local.sh 1-36       # a subset (e.g. VM 1 of 2)
#
# Environment (optional): PB_NCORES (default detectCores()-2), PB_REPS (200),
# PB_SEED (20261005), PB_TIMEOUT (default 21d), CMDSTAN (CmdStan path).
#
# Safe to re-launch after a disconnect or crash: a pid file prevents two runs
# on the same rows, and SimDesign resumes from its temp file in the chunk
# directory (finished Design rows are skipped). A finished chunk is never
# re-run (res_rows_*.rds guard).
# ------------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

ROWS="${1:-${PB_ROWS:-all}}"
TAG="$(echo "$ROWS" | tr ',' '_')"
TIMEOUT="${PB_TIMEOUT:-21d}"
mkdir -p logs
LOG="logs/t5_full.log"
[ "$ROWS" != "all" ] && LOG="logs/t5_full_rows_${TAG}.log"
PIDFILE="logs/t5_full_rows_${TAG}.pid"

if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "Already running (pid $(cat "$PIDFILE")) for rows '$ROWS'; log: $LOG"
  exit 0
fi

{
  echo "==== launch $(date '+%Y-%m-%dT%H:%M:%S%z') on $(hostname) | rows=$ROWS | timeout=$TIMEOUT"
  uptime
} >> "$LOG"

nohup bash -c "timeout '$TIMEOUT' Rscript scripts/13_BFCalibration_RS.R full '$ROWS'; \
  echo \"==== EXIT \$? \$(date '+%Y-%m-%dT%H:%M:%S%z')\"" >> "$LOG" 2>&1 < /dev/null &
echo $! > "$PIDFILE"
echo "Started pid $(cat "$PIDFILE"); log: $LOG"
echo "Follow with: tail -f $LOG   (finished when a line '==== EXIT 0' appears)"
