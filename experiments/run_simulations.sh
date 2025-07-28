#!/usr/bin/env bash

# ./run_simulations.sh log      0.16253
# ./run_simulations.sh identity 0.16253

set -euo pipefail

# Usage check
if [ $# -ne 2 ]; then
  echo "Usage: $0 <link_type> <baseline>"
  echo "  <link_type> : log or identity"
  echo "  <baseline>  : target p_base (e.g. 0.16253)"
  exit 1
fi

link_type=$1
baseline=$2

echo "Starting simulation run at $(date)"
echo "  link_type = $link_type"
echo "  baseline  = $baseline"

cd /home/rstudio/teh_detection/experiments

# Timestamped logfile
logfile="trials_${link_type}_$(date +%Y%m%d_%H%M%S).log"

# Launch in background
nohup Rscript run_trials_parallel.R \
  --link_type "$link_type" \
  --baseline  "$baseline" \
  > "$logfile" 2>&1 &

echo "Job submitted with nohup."
echo "  Logs → simulations/$logfile"

