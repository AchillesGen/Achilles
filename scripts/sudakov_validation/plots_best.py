"""Figure 9: Achilles Base vs the best Track A and Track B set-ups, Achilles default NN input.
Figure 10: the same with NASA (Tot) input, with the published Base curve of arXiv:2007.15570.
sigma_R from the 14 fm beam-radius runs (Base is not converged at 7 fm below ~40 MeV);
any-collision hit definition. Transparency counted in acceptance.
Usage: python3 plots_best.py <val dir> <fig dir> [<paper Fig. 6 QMC Gauss Tot curve csv>]"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chi2_table import regroup, load, chi2, EL, EH, ET

VAL, FIG = sys.argv[1:3]
PAPER = sys.argv[3] if len(sys.argv) > 3 else None
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "legend.frameon": False, "axes.titlesize": 11, "axes.titleweight": "bold"})
car, tra = load(os.path.join(VAL, "carlson.csv")), load(os.path.join(VAL, "transparency.csv"))
low = car[:, 0] < 100
SERIES = [  # key, label, colour, style
    ("Base", "Achilles Base", INK2, "--"),
    ("B_eik_Gauss_single_ipm", "Track B (best)", ORANGE, "-"),
    ("A_fold_finiteA_ipm", "Track A (best)", BLUE, "-"),
]


def jl(p):
    return json.load(open(os.path.join(VAL, p)))


def figure(sig_json, tr_jsons, nn_name, out, paper=None):
    S = jl(sig_json)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    if paper is not None:
        pc = np.loadtxt(paper, delimiter=",", comments="#")
        pc = pc[(pc[:, 1] > 0) & (pc[:, 0] > 8)]
        ax[0].plot(pc[:, 0], pc[:, 1], color=INK2, ls=":", lw=1.6,
                   label="published Base curve (Fig. 6; 2020 code, escape rule)")
    for key, lab, col, ls in SERIES:
        g = np.vstack([regroup(S[f"sigR_low_{key}"], EL, "s"), regroup(S[f"sigR_high_{key}"], EH, "s")])
        c_lo, c_hi = chi2(car, g, low)[0], chi2(car, g, ~low)[0]
        ax[0].plot(g[:, 0], g[:, 1], color=col, ls=ls, lw=2,
                   label=f"{lab}   $\\chi^2/N$: {c_lo:.0f} (<100 MeV), {c_hi:.0f} (>100)")
        ax[0].fill_between(g[:, 0], g[:, 1] - g[:, 2], g[:, 1] + g[:, 2], color=col, alpha=0.15, lw=0)
        T = jl(tr_jsons[key])
        t = regroup(T[f"transp_{key}"], ET, "a")
        c_t = chi2(tra, t, np.ones(len(tra), bool), logx=False)[0]
        ax[1].plot(t[:, 0], t[:, 1], color=col, ls=ls, lw=2, label=f"{lab}   $\\chi^2/N$: {c_t:.1f}")
        ax[1].fill_between(t[:, 0], t[:, 1] - t[:, 2], t[:, 1] + t[:, 2], color=col, alpha=0.15, lw=0)
    ax[0].errorbar(car[:, 0], car[:, 1], yerr=car[:, 2], fmt="o", color=INK, ms=3.2, mfc="white", mew=0.9,
                   lw=0.7, zorder=5, label="data: Carlson (1996)")
    ax[1].errorbar(tra[:, 0], tra[:, 1], yerr=tra[:, 2], fmt="o", color=INK, ms=3.8, mfc="white", mew=0.9,
                   lw=0.7, zorder=5, label="data: $^{12}$C(e,e$'$p)")
    ax[0].set_xscale("log"); ax[0].set_ylim(0, 700)
    ax[0].set_xlabel("proton kinetic energy $T_p$ [MeV]"); ax[0].set_ylabel("$\\sigma_R$ [mb]")
    ax[0].set_title("p + $^{12}$C reaction cross section")
    ax[1].set_ylim(0.4, 0.9); ax[1].set_xlabel("proton kinetic energy $T_p$ [MeV]"); ax[1].set_ylabel("transparency")
    ax[1].set_title("$^{12}$C transparency (in acceptance)")
    ax[0].legend(fontsize=7.8, loc="upper right"); ax[1].legend(fontsize=8, loc="upper right")
    fig.text(0.5, 0.005, f"All use {nn_name} NN cross sections; $\\sigma_R$ counts any collision, beam radius 14 fm. "
             "Track A: density rate + fold + finite-A + IPM coherent veto. "
             "Track B: eikonal Gaussian pairs, single Pauli attempt, IPM veto. Bands: MC statistics.",
             ha="center", fontsize=7.2, color=INK2)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(os.path.join(FIG, out), dpi=200)


figure("data_r14_nn/r14_nn.json", {"Base": "data/data_runs_nn_step2.json", "B_eik_Gauss_single_ipm": "data_ipm_nn/ipm_nn.json",
       "A_fold_finiteA_ipm": "data_ipm_nn/ipm_nn.json"}, "Achilles' default", "fig9_best_vs_base.png")
figure("data_r14/r14_nasa.json", {"Base": "data_acc/acc_nasa.json", "B_eik_Gauss_single_ipm": "data_ipm/ipm_nasa.json",
       "A_fold_finiteA_ipm": "data_ipm/ipm_nasa.json"}, "NASA (Tot)", "fig10_best_vs_base_nasa.png", PAPER)
print("ok")
