import pyhepmc
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="A weight plotter")
parser.add_argument("filename", type=str, help="Event file to plot")
args = parser.parse_args()

stem = Path(args.filename).stem

weights = []
with pyhepmc.open(args.filename) as f:
    for event in f:
        weights.append(event.weight("CV"))

wgt = np.array(weights)
mean = np.mean(wgt)
max_wgt = np.percentile(wgt, 99)
eff = mean/max_wgt*100

fig, ax = plt.subplots()
plt.hist(np.log10(wgt), bins=500, histtype='step')
plt.axvline(x=np.log10(mean), color='red', ls='--')
plt.axvline(x=np.log10(max_wgt), color='green', ls='--')
plt.text(0.05, 0.95, f"Eff = {eff:.2f}%", transform=ax.transAxes)
plt.ylabel('# Events')
plt.xlabel(r'$log_{10}(wgt)$')
plt.yscale('log')
plt.savefig(f'weights_{stem}.png')
plt.show()
