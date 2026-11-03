"""Figure 7: Pauli blocking on/off, the coherent F^2(q) veto, the in-medium correction and the
acceptance-defined transparency, NASA ("Tot") NN cross sections.
Usage: python3 plots_step2.py <val dir> <fig dir>
<val dir> holds carlson.csv, transparency.csv, data_nasa/*.json and data_acc/acc_nasa.json."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chi2_table import regroup, load, EL, EH, ET

VAL, FIG = sys.argv[1:3]
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 8.5,
    "legend.frameon": False, "lines.linewidth": 1.5, "axes.titlesize": 9, "axes.titleweight": "bold"})

J = {k: json.load(open(os.path.join(VAL, p))) for k, p in {
    "old": "data_nasa/data_runs_nasa.json", "sud": "data_nasa/data_runs_nasa_sudakov.json",
    "new": "data_nasa/data_runs_nasa_step2.json", "acc": "data_acc/acc_nasa.json"}.items()}
car, tra = load(os.path.join(VAL, "carlson.csv")), load(os.path.join(VAL, "transparency.csv"))


def sig(src, a):
    r = J[src]
    return np.vstack([regroup(r[f"sigR_low_{a}"], EL, "s"), regroup(r[f"sigR_high_{a}"], EH, "s")])


def draw(ax, g, lab, col, ls="-"):
    ax.errorbar(g[:, 0], g[:, 1], yerr=g[:, 2], color=col, ls=ls, marker="o", ms=2.5, lw=1.4, label=lab)


def sig_axes(ax, title, ymax=700):
    ax.errorbar(car[:, 0], car[:, 1], yerr=car[:, 2], fmt="s", color=INK, ms=2.8, mfc="white", mew=0.8,
                lw=0.7, zorder=5, label="data: Carlson (1996), digitised")
    ax.set_xscale("log"); ax.set_ylim(0, ymax)
    ax.set_xlabel("proton kinetic energy $T_p$ [MeV]"); ax.set_ylabel("$\\sigma_R$ [mb]")
    ax.set_title(title); ax.legend(fontsize=6.6, loc="upper right")


fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.0))
# (a) Pauli on/off
a = ax[0, 0]
draw(a, sig("sud", "B_eik_Gauss_single"), "B eikonal Gauss, Pauli (single attempt)", C[1])
draw(a, sig("new", "B_eik_Gauss_nopauli"), "B eikonal Gauss, no Pauli", C[1], "--")
draw(a, sig("sud", "A_fold_finiteA"), "A fold + finite-A, Pauli", C[0])
draw(a, sig("new", "A_fold_finiteA_nopauli"), "A fold + finite-A, no Pauli", C[0], "--")
sig_axes(a, "(a) Pauli blocking on / off", ymax=1150)
# (b) coherent veto and in-medium, Track B
a = ax[0, 1]
draw(a, sig("old", "Base"), "Achilles Base, Gaussian", INK2, ":")
draw(a, sig("sud", "B_eik_Gauss_single"), "B eikonal Gauss (single)", C[1])
draw(a, sig("new", "B_eik_Gauss_single_coh"), "  + coherent $F^2$ veto", C[2])
draw(a, sig("new", "B_eik_Gauss_single_coh_med"), "  + coherent + in-medium", C[6])
draw(a, sig("new", "B_k1_Box_single_coh"), "B $\\kappa$=1 box + coherent", C[4], "--")
sig_axes(a, "(b) Track B: coherent veto, in-medium $\\sigma$")
# (c) Track A and Base
a = ax[1, 0]
draw(a, sig("old", "Base"), "Achilles Base, Gaussian", INK2, ":")
draw(a, sig("new", "Base_med"), "Achilles Base, Gaussian, in-medium", INK2, "-.")
draw(a, sig("sud", "A_fold_finiteA"), "A fold + finite-A", C[0])
draw(a, sig("new", "A_fold_finiteA_coh"), "  + coherent $F^2$ veto", C[2])
draw(a, sig("new", "A_fold_finiteA_coh_med"), "  + coherent + in-medium", C[6])
sig_axes(a, "(c) Track A and Achilles Base")
# (d) transparency: any interaction vs acceptance
a = ax[1, 1]
for (key, lab, col) in [("Base", "Achilles Base, Gaussian", INK2), ("B_eik_Gauss_single", "B eikonal Gauss (single)", C[1]),
                        ("B_eik_Gauss_single_coh", "B + coherent", C[2]), ("A_fold_finiteA_coh", "A + coherent", C[0])]:
    tr = J["acc"][f"transp_{key}"]
    draw(a, regroup(tr, ET, "t"), lab + ": no interaction", col, ":")
    draw(a, regroup(tr, ET, "a"), lab + ": in acceptance", col, "-")
a.errorbar(tra[:, 0], tra[:, 1], yerr=tra[:, 2], fmt="s", color=INK, ms=3, mfc="white", mew=0.8, lw=0.7,
           zorder=5, label="data: $^{12}$C(e,e$'$p), digitised")
a.set_ylim(0.3, 1.0); a.set_xlabel("proton kinetic energy $T_p$ [MeV]"); a.set_ylabel("transparency")
a.set_title("(d) transparency: no interaction vs in acceptance")
a.legend(fontsize=6.2, loc="lower left", ncol=2)
a.text(0.98, 0.97, "acceptance: $\\Delta E<80$ MeV, $|\\Delta p|<300$ MeV/c", transform=a.transAxes, ha="right", va="top", fontsize=7, color=INK2)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig7_pauli_coherent_medium.png"), dpi=200)
print("ok")
