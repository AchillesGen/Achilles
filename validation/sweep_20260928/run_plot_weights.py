# Runs plot_weights.py unmodified and prints the numbers it computes.
import os, sys, runpy
import matplotlib
matplotlib.use("Agg")
script, events, cutoff, outdir = sys.argv[1:5]
events = os.path.abspath(events)
os.makedirs(outdir, exist_ok=True)
os.chdir(outdir)
sys.argv = [script, events, "--cutoff", cutoff]
ns = runpy.run_path(script, run_name="__main__")
w = ns["wgt"]
print(f"{os.path.basename(events)} cutoff={cutoff}% n={len(w)} mean={ns['mean']:.6e} "
      f"max_wgt={ns['max_wgt']:.6e} eff={ns['eff']:.3f}% relstd={w.std()/w.mean():.3f}")
