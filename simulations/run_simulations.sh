#!/bin/bash

# Run simulations in the background using nohup
echo "Starting simulation run at $(date)"
cd /home/rstudio/teh_detection/simulations
nohup Rscript run_trials_parallel.R > trials.log 2>&1 &
echo "Job submitted with nohup. Output will go to simulations/trials.log"

