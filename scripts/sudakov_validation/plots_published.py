"""Figure 6: the published Achilles cascade set-ups (Isaacson et al., PRC 103, 015502):
Base algorithm with SAID/GEANT4 elastic ("El") and NASA total ("Tot") NN cross sections,
Gaussian and cylinder probabilities, against the Sudakov tracks run with the same input.
Usage: python3 plots_published.py <nasa.json> <nasa_sudakov.json> <geant.json> <carlson.csv> <transparency.csv> <fig dir>"""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NASA, NASA_S, GEANT, CARL, TRANS, FIG = sys.argv[1:7]
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
    "legend.frameon": False, "lines.linewidth": 1.6, "axes.titlesize": 9.5, "axes.titleweight": "bold"})
MN = 938.92
def T_of_p(p): return np.sqrt(p**2 + MN**2) - MN

def regroup(e, edges, mode):
    lo, hi = np.array(e["lo"]), np.array(e["hi"]); h, n = np.array(e["hits"]), np.array(e["nohits"])
    Tc = T_of_p(0.5*(lo+hi)); out = []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (Tc >= a) & (Tc < b); H, N = h[s].sum(), n[s].sum(); t = H+N
        if t < 50: continue
        f = H/t; er = math.sqrt(f*(1-f)/t); Tm = np.average(Tc[s], weights=(h+n)[s])
        out.append((Tm, 10*math.pi*49*f, 10*math.pi*49*er) if mode == "s" else (Tm, 1-f, er))
    return np.array(out)

def load_csv(fn):
    return np.array([[float(x) for x in l.split(",")] for l in open(fn) if not (l.startswith("#") or l.startswith("Tp"))])

el = [8, 12, 16, 20, 25, 30, 36, 43, 52, 62, 75, 90, 105]
eh = [105, 130, 160, 200, 250, 310, 380, 460, 560, 680, 820, 980]
et = [140, 220, 300, 380, 460, 560, 680, 820, 980, 1150, 1350, 1550, 1800, 1950]
J = {"nasa": json.load(open(NASA)), "nasa_s": json.load(open(NASA_S)), "geant": json.load(open(GEANT))}
SERIES = [  # (json, alg key, label, colour, linestyle)
    ("geant", "Base", "Achilles Base, Gaussian, El (SAID/GEANT4)", INK2, ":"),
    ("nasa", "Base", "Achilles Base, Gaussian, Tot (NASA)", INK, "--"),
    ("nasa", "Base_cyl", "Achilles Base, cylinder, Tot (NASA)", INK, "-."),
    ("nasa_s", "B_eik_Gauss", "B eikonal Gauss, Pauli retry, Tot", C[1], "-"),
    ("nasa_s", "B_eik_Gauss_single", "B eikonal Gauss, Pauli single attempt, Tot", C[2], "-"),
    ("nasa_s", "A_fold_finiteA", "A, $\\rho\\otimes\\hat p$ + finite-A, Tot", C[0], "-"),
]
car, tra = load_csv(CARL), load_csv(TRANS)
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.3))
for src, a, lab, col, ls in SERIES:
    r = J[src]
    g = np.vstack([regroup(r[f"sigR_low_{a}"], el, "s"), regroup(r[f"sigR_high_{a}"], eh, "s")])
    ax[0].errorbar(g[:, 0], g[:, 1], yerr=g[:, 2], color=col, ls=ls, marker="o", ms=2.8, lw=1.5, label=lab)
    g = regroup(r[f"transp_{a}"], et, "t")
    ax[1].errorbar(g[:, 0], g[:, 1], yerr=g[:, 2], color=col, ls=ls, marker="o", ms=2.8, lw=1.5, label=lab)
ax[0].errorbar(car[:, 0], car[:, 1], yerr=car[:, 2], fmt="s", color=INK, ms=3, mfc="white", mew=0.9, lw=0.8, zorder=5,
               label="data: Carlson (1996), digitised")
ax[1].errorbar(tra[:, 0], tra[:, 1], yerr=tra[:, 2], fmt="s", color=INK, ms=3.5, mfc="white", mew=0.9, lw=0.8, zorder=5,
               label="data: $^{12}$C(e,e$'$p), digitised")
ax[0].set_xscale("log"); ax[0].set_ylim(0, 700)
ax[0].set_xlabel("proton kinetic energy $T_p$ [MeV]"); ax[0].set_ylabel("$\\sigma_R$ [mb]")
ax[1].set_xlabel("proton kinetic energy $T_p$ [MeV]"); ax[1].set_ylabel("transparency"); ax[1].set_ylim(0.3, 1.0)
ax[0].set_title("(a) p + $^{12}$C $\\sigma_R$: published set-ups vs Sudakov tracks")
ax[1].set_title("(b) $^{12}$C transparency")
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.3, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.14, 1, 1))
fig.savefig(os.path.join(FIG, "fig6_published_setups.png"), dpi=200)
print("ok")
