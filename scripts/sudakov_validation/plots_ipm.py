"""Figure 8: the IPM-corrected coherent veto.
Usage: python3 plots_ipm.py <val dir> <density file> <fig dir>"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chi2_table import regroup, load, EL, EH, ET
import ipm_veto

VAL, DENS, FIG = sys.argv[1:4]
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 8.5,
    "legend.frameon": False, "lines.linewidth": 1.5, "axes.titlesize": 9, "axes.titleweight": "bold"})
J = {k: json.load(open(os.path.join(VAL, p))) for k, p in {
    "nasa_old": "data_nasa/data_runs_nasa.json", "nasa_sud": "data_nasa/data_runs_nasa_sudakov.json",
    "nasa_new": "data_nasa/data_runs_nasa_step2.json", "nasa_ipm": "data_ipm/ipm_nasa.json",
    "nn_old": "data/data_runs.json", "nn_new": "data/data_runs_nn_step2.json", "nn_ipm": "data_ipm_nn/ipm_nn.json",
    "nn_acc": "data_acc_nn/acc_nn.json"}.items()}
car, tra = load(os.path.join(VAL, "carlson.csv")), load(os.path.join(VAL, "transparency.csv"))


def sig(src, a):
    r = J[src]
    return np.vstack([regroup(r[f"sigR_low_{a}"], EL, "s"), regroup(r[f"sigR_high_{a}"], EH, "s")])


def draw(ax, g, lab, col, ls="-"):
    ax.errorbar(g[:, 0], g[:, 1], yerr=g[:, 2], color=col, ls=ls, marker="o", ms=2.5, lw=1.4, label=lab)


def sig_axes(ax, title):
    ax.errorbar(car[:, 0], car[:, 1], yerr=car[:, 2], fmt="s", color=INK, ms=2.8, mfc="white", mew=0.8,
                lw=0.7, zorder=5, label="data: Carlson (1996), digitised")
    ax.set_xscale("log"); ax.set_ylim(0, 650)
    ax.set_xlabel("proton kinetic energy $T_p$ [MeV]"); ax.set_ylabel("$\\sigma_R$ [mb]")
    ax.set_title(title); ax.legend(fontsize=6.6, loc="upper right")


fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.0))
# (a) response functions
a = ax[0, 0]
d = np.loadtxt(DENS, comments="#"); r, rho = d[:, 0], d[:, 1]
rms = np.sqrt(np.trapezoid(rho * r**4, r) / np.trapezoid(rho * r**2, r)); b = rms / np.sqrt(26 / 12)
qs = np.linspace(5, 600, 120)
S = np.array([ipm_veto.s_inc(q / ipm_veto.HBARC, b)[0] for q in qs])
F = np.array([ipm_veto.s_inc(q / ipm_veto.HBARC, b)[1] for q in qs])
L = np.array([ipm_veto.lfg_allowed(q / ipm_veto.HBARC, r, rho) for q in qs])
a.plot(qs, S, color=INK, label="exact IPM: nucleus left excited, $S_{\\rm inc}(q)$")
a.plot(qs, L, color=C[0], ls="--", label="cascade LFG Pauli: allowed fraction $L(q)$")
a.plot(qs, L * (1 - F**2), color=C[1], ls="--", label="LFG $\\times\\,(1-F^2)$: the $F^2$ rule")
a.plot(qs, F**2, color=C[1], ls=":", label="$F^2(q)$ veto")
a.plot(qs, np.clip(1 - S / L, 0, 1), color=C[2], label="IPM veto $v(q)=1-S_{\\rm inc}/L$")
a.set_xlabel("momentum transfer $q$ [MeV/c]"); a.set_ylabel("probability"); a.set_ylim(0, 1.02)
a.set_title("(a) $^{12}$C: what the Pauli factor already removes"); a.legend(fontsize=6.6, loc="center right")
# (b) NASA sigma_R
a = ax[0, 1]
draw(a, sig("nasa_old", "Base"), "Achilles Base, Gaussian", INK2, ":")
draw(a, sig("nasa_sud", "A_fold_finiteA"), "A: no coherent veto", C[0], ":")
draw(a, sig("nasa_new", "A_fold_finiteA_coh"), "A: $F^2$ veto", C[0], "--")
draw(a, sig("nasa_ipm", "A_fold_finiteA_ipm"), "A: IPM veto", C[0], "-")
draw(a, sig("nasa_sud", "B_eik_Gauss_single"), "B: no coherent veto", C[1], ":")
draw(a, sig("nasa_ipm", "B_eik_Gauss_single_ipm"), "B: IPM veto", C[1], "-")
sig_axes(a, "(b) $\\sigma_R$, NASA (Tot) NN input")
# (c) NN set sigma_R
a = ax[1, 0]
draw(a, sig("nn_old", "Continuous"), "Achilles Continuous", INK2, "--")
draw(a, sig("nn_new", "Base"), "Achilles Base, Gaussian", INK2, ":")
draw(a, sig("nn_old", "A_fold_finiteA"), "A: no coherent veto", C[0], ":")
draw(a, sig("nn_ipm", "A_fold_finiteA_ipm"), "A: IPM veto", C[0], "-")
draw(a, sig("nn_ipm", "B_eik_Gauss_single_ipm"), "B: IPM veto", C[1], "-")
sig_axes(a, "(c) $\\sigma_R$, Achilles default NN set")
# (d) transparency, NN set, in acceptance
a = ax[1, 1]
for src, key, lab, col, ls in [("nn_acc", "Continuous", "Achilles Continuous", INK2, "--"),
                               ("nn_acc", "Base", "Achilles Base", INK2, ":"),
                               ("nn_ipm", "A_fold_finiteA_ipm", "A: IPM veto", C[0], "-"),
                               ("nn_ipm", "B_eik_Gauss_single_ipm", "B: IPM veto", C[1], "-")]:
    draw(a, regroup(J[src][f"transp_{key}"], ET, "a"), lab, col, ls)
a.errorbar(tra[:, 0], tra[:, 1], yerr=tra[:, 2], fmt="s", color=INK, ms=3, mfc="white", mew=0.8, lw=0.7,
           zorder=5, label="data: $^{12}$C(e,e$'$p), digitised")
a.set_ylim(0.3, 1.0); a.set_xlabel("proton kinetic energy $T_p$ [MeV]"); a.set_ylabel("transparency (in acceptance)")
a.set_title("(d) transparency, default NN set"); a.legend(fontsize=6.6, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig8_ipm_veto.png"), dpi=200)
print("ok")
