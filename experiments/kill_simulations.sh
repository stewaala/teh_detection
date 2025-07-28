#!/bin/bash
#
# kill_trials.sh
#  - Kill all R processes running run_trials_parallel.R (whatever cwd)

set -euo pipefail

SCRIPT_NAME="run_trials_parallel\\.R"

echo "🔍 Killing any R processes invo­ked with “$SCRIPT_NAME”…"

pkill -TERM -f "$SCRIPT_NAME" \
  && echo "🛑 Sent TERM to all matching processes." \
  || echo "✔ No matching processes found."

echo "✅ Done."

