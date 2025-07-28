#!/usr/bin/env bash
set -euo pipefail

# Usage: run_simulations.sh <link_type>
if [ $# -ne 1 ]; then
  echo "Usage: $0 <link_type>"
  echo "  <link_type> : log or identity"
  exit 1
fi

link_type=$1

echo "Starting simulation run at $(date)"
echo "  link_type = $link_type"

cd /home/rstudio/teh_detection/experiments

logfile="trials_${link_type}_$(date +%Y%m%d_%H%M%S).log"

nohup Rscript run_trials_parallel.R "$link_type" > "$logfile" 2>&1 &

echo "Job submitted with nohup."
echo "  Logs → $PWD/$logfile"

