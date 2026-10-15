#!/bin/bash
cd "$(dirname "$0")/../.."
for mode in qe res; do
  for acc in 1.0e-3 1.0e-2; do
    # Comparable timing needs a quiet machine (the final_* flow runs ran at load ~1).
    while (( $(cut -d' ' -f1 /proc/loadavg | cut -d. -f1) >= 3 )); do sleep 60; done
    echo "start mc_${mode}_acc${acc} load=$(cut -d' ' -f1 /proc/loadavg)"
    python validation/sweep_20260928/driver.py mc_${mode}_acc${acc} $mode 150 5000 \
      "{\"Kind\":\"MultiChannel\",\"Accuracy\":\"$acc\"}" 2>&1 | grep -E 'RESULT|Error|Traceback'
  done
done
echo MC_DONE
