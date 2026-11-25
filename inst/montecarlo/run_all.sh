#!/bin/bash
# Runs every Monte Carlo design in order of priority for Paper 1.
cd "$(dirname "$0")"
export MC_MIRROR=${MC_MIRROR:-./results/raw}
for d in main_n2000 main_n500 misspec weakiv j3 main_large; do
  echo "=== $d $(date -u +%H:%M:%S)" >> mc_run.log
  Rscript mc_copula.R $d 4 >> mc_run.log 2>&1
done
echo "=== ALL DONE $(date -u +%H:%M:%S)" >> mc_run.log
