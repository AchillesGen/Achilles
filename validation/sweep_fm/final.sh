#!/bin/bash
# Interleaved, sequential final comparison: VegasFlowMatchingFlow (defaults) vs VegasSplineFlow (qe_p2).
cd "$(dirname "$0")/../.."
for i in 1 2 3; do
  python validation/sweep_20260928/driver.py final_fm_$i qe 150 5000 \
    '{"Class":"VegasFlowMatchingFlow","hidden":64,"step_size":0.05,"vegas_iterations":10}' 2>&1 | grep RESULT
  python validation/sweep_20260928/driver.py final_spline_$i qe 150 5000 \
    '{"Class":"VegasSplineFlow","hidden":32,"bins":16,"vegas_iterations":10,"loss_power":2.0}' 2>&1 | grep RESULT
done
echo FINAL_DONE
