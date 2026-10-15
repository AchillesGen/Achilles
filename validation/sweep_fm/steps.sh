#!/bin/bash
cd "$(dirname "$0")/../.."
for h in 0.5 1.0; do
  python validation/sweep_20260928/driver.py final_fm_h$h qe 150 5000 \
    "{\"Class\":\"VegasFlowMatchingFlow\",\"hidden\":64,\"step_size\":$h,\"vegas_iterations\":10}" 2>&1 | grep RESULT
done
