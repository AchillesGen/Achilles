"""Write comparison configuration files for achilles-cascade.

  <out>/12C_mixed_configs.out.gz    same rho_p, rho_n as QMC, no correlations
  <out>/12C_jastrow_configs.out.gz  same rho_p, rho_n plus the QMC g(r), nothing else

Run: python3 scripts/cascade_studies/make_configs.py <outdir>
"""

import pathlib
import sys

import numpy as np

import configs

rng = np.random.default_rng(7)
out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".")
out.mkdir(parents=True, exist_ok=True)

pos, isp, w = configs.load("QMC_configs.out.gz")
pos, isp = configs.unweight(pos, isp, w, rng)
# Order nucleons in every configuration as protons first so slots share isospin
order = np.argsort(~isp, axis=1, kind="stable")
pos = np.take_along_axis(pos, order[..., None], axis=1)
isp = np.take_along_axis(isp, order, axis=1)
assert (isp == isp[0]).all()

nout = 36000
r, g = configs.g_of_r(pos, rng)
configs.write(out / "12C_mixed_configs.out.gz", configs.mixed(pos, rng, nout, isp), isp[0])
configs.write(out / "12C_jastrow_configs.out.gz", configs.jastrow(pos, r, g, rng, nout, isp=isp), isp[0])
print("unweighted QMC configs:", len(pos))
