"""Figures for the Sudakov validation (Step 4 a-f) and the data comparison (Step 4 h).
Usage: python3 plots.py <forward out dir> <compare_NN.json> <data_runs.json or -> <carlson.csv> <transparency.csv> <fig dir>"""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import special, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
FWD, CMPJ, DATAJ, CARL, TRANS, FIG = sys.argv[1:7]
os.makedirs(FIG, exist_ok=True)

# reference palette (light), fixed order; data in primary ink
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
    "legend.frameon": False, "lines.linewidth": 1.6, "lines.markersize": 5,
    "axes.titlesize": 9.5, "axes.titleweight": "bold",
})
MN = 938.92

def T_of_p(p): return np.sqrt(p**2 + MN**2) - MN
def load(name): return np.genfromtxt(os.path.join(FWD, name + ".csv"), delimiter=",", names=True)

A12, RMS = 12, 2.33
R = math.sqrt(5/3)*RMS; rho0 = A12/(4/3*math.pi*R**3); mu40 = 4.0*rho0
def ext_Pn(n, mu):
    y = 2*mu*R; return (n+1)*special.gammainc(n+2, y)/(2*mu**2*R**2)
def int_T(mu):
    y = 2*mu*R
    return 3/(4*R**3)*(R**2*(1-np.exp(-y))/mu - (2-np.exp(-y)*(y*y+2*y+2))/(4*mu**3))

# ---------------------------------------------------------------- Figure 1: Step 4(a), (b)
def fig_closed_forms():
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
    d = load("a_ext_A"); N = len(d)
    n = np.arange(0, 9); cnt = np.bincount(d["nint"].astype(int), minlength=9)[:9]/N
    ex = np.array([ext_Pn(k, mu40) for k in n])
    ax[0].bar(n, ex, width=0.6, color=C[0], alpha=0.35, label="closed form", zorder=2)
    ax[0].errorbar(n, cnt, yerr=np.sqrt(cnt*(1-cnt)/N), fmt="o", color=INK, ms=4, label="engine (Track A)", zorder=3)
    ax[0].set_xlabel("number of collisions $N$"); ax[0].set_ylabel("$P(N)$, $b<R$")
    ax[0].set_title("(a) collision number, uniform sphere, 40 mb"); ax[0].legend()

    dep = d["depth"][d["nint"] > 0]
    Lg = np.linspace(0, 2*R, 4001); fL = Lg/(2*R**2); pint = 1-ext_Pn(0, mu40)
    edges = np.linspace(0, 2*R, 31); cen = 0.5*(edges[1:]+edges[:-1]); wbin = edges[1]-edges[0]
    h, _ = np.histogram(dep, bins=edges)
    dens = h/len(dep)/wbin; err = np.sqrt(h)/len(dep)/wbin
    xs = np.linspace(0, 2*R, 300)
    pdf = [np.trapezoid(fL*(Lg > x)*mu40*np.exp(-mu40*x), Lg)/pint for x in xs]
    ax[1].plot(xs, pdf, color=C[0], lw=2.0, label="closed form")
    ax[1].errorbar(cen, dens, yerr=err, fmt="o", color=INK, ms=3.5, label="engine")
    ax[1].set_xlabel("depth of first collision [fm]"); ax[1].set_ylabel("probability density [fm$^{-1}$]")
    ax[1].set_title("(b) first-interaction point (KS p = 0.41)"); ax[1].legend(loc="upper right")

    # toy A=3: first-collision probabilities per profile
    from analyse_helpers import toy_exact, TOYS
    profs = ["Atom", "Box", "Gauss"]; w = 0.25
    for k, prof in enumerate(profs):
        dd = load(f"b_toy3_{prof}"); N = len(dd)
        P, _ = toy_exact(TOYS[3], prof)
        obs = np.array([np.mean(dd["first"] == j) for j in range(3)])
        x = np.arange(3) + (k-1)*w
        ax[2].bar(x, P, width=w*0.9, color=C[k], alpha=0.45, label=f"{prof}: exact", zorder=2)
        ax[2].errorbar(x, obs, yerr=np.sqrt(obs*(1-obs)/N), fmt="o", color=INK, ms=3.5, zorder=3,
                       label="engine" if k == 0 else None)
    ax[2].set_xticks(range(3)); ax[2].set_xticklabels(["nucleon 1", "nucleon 2", "nucleon 3"])
    ax[2].set_ylabel("P(first collision with j)"); ax[2].set_yscale("log")
    ax[2].set_title("(c) toy A = 3, Track B, three time profiles"); ax[2].legend(fontsize=7.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_closed_forms.png"), dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 2: Step 4(c), (d)
def fig_majorant_weights():
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.3))
    s = np.linspace(20, 60, 200)
    ax[0].plot(s, [int_T(x*0.1*rho0) for x in s], color=C[0], lw=2.0, label="closed form $T(\\sigma)$")
    for c, mk, col in (("1.25", "s", C[1]), ("2", "D", C[2]), ("4", "^", C[6])):
        d = load(f"d_int_A_c{c}"); N = len(d)
        for k, f in ((0, 0.8), (1, 1.2)):
            wcol = d[f"w{k}"]*(d["nint"] == 0)
            ax[0].errorbar(40*f + (float(c)-2)*0.25, wcol.mean(), yerr=wcol.std()/math.sqrt(N), fmt=mk,
                           color=col, ms=5, label=f"reweighted from 40 mb, headroom {c}" if k == 0 else None)
    d = load("a_int_A"); T = np.mean(d["nint"] == 0)
    ax[0].errorbar(40, T, yerr=math.sqrt(T*(1-T)/len(d)), fmt="o", color=INK, label="generated at 40 mb")
    d = load("c_int_A_h08"); Tb = np.mean(d["nint"] == 0)
    ax[0].errorbar(40, Tb, yerr=math.sqrt(Tb*(1-Tb)/len(d)), fmt="x", color=C[7], ms=7, mew=1.6,
                   label="majorant 0.8x (counter fires)")
    ax[0].set_xlabel("$\\sigma$ [mb]"); ax[0].set_ylabel("transparency, internal source")
    ax[0].set_title("(a) reweighting closure and violated majorant"); ax[0].legend(fontsize=7)

    # ESS vs headroom
    Lg = np.linspace(0, 2*R, 4001); fL = 3/(4*R**3)*(R**2-Lg**2/4)
    cs = np.linspace(1.05, 6, 120)
    for f, col in ((0.8, C[2]),):
        pred = []
        for c in cs:
            a = 1/(1.01*1.2*c); ap = f*a
            pred.append(1/np.trapezoid(fL*np.exp(1.01*1.2*c*mu40*Lg*(ap-a)**2/(a*(1-a))), Lg))
        ax[1].plot(cs, pred, color=col, lw=1.8, label="predicted $1/E[W^2]$ ($\\sigma\\times0.8$ and $\\times1.2$ coincide)")
    for c in ("1.25", "2", "4"):
        d = load(f"d_int_A_c{c}")
        for k, col in ((0, C[1]), (1, C[2])):
            w = d[f"w{k}"]; ax[1].plot(float(c), w.sum()**2/(w**2).sum()/len(w), "o", color=INK, ms=5)
    ax[1].plot([], [], "o", color=INK, label="engine")
    ax[1].set_xlabel("headroom $c=\\hat\\Gamma/\\Gamma_{\\rm nominal}$"); ax[1].set_ylabel("ESS / N")
    ax[1].set_title("(b) weight variance vs majorant headroom"); ax[1].legend(fontsize=7.5, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_majorant_weights.png"), dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 3: Step 4(e)
def fig_qmc():
    rows = [("A on $\\rho$ (optical)", "e_ext_A_rho", 229.6),
            ("A on $\\rho\\otimes\\hat p$", "e_ext_A_fold", 266.4),
            ("A on $\\rho\\otimes\\hat p$ + finite-A", "e_ext_A_foldFA", 271.9),
            ("B, $\\kappa$=1, decorrelated", "e_ext_B_k1_dec1", 271.9),
            ("B, $\\kappa$=1, QMC", "e_ext_B_k1_dec0", 275.6),
            ("B, eikonal, decorrelated", "e_ext_B_eik_dec1", 266.3),
            ("B, eikonal, QMC", "e_ext_B_eik_dec0", 270.0)]
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    for i, (lab, name, ref) in enumerate(rows):
        d = load(name); p = np.mean(d["nint"] > 0); N = len(d)
        s = 10*math.pi*64*p; e = 10*math.pi*64*math.sqrt(p*(1-p)/N)
        col = C[0] if lab.startswith("A") else C[1]
        ax.errorbar(s, i, xerr=e, fmt="o", color=col, ms=5)
        ax.plot(ref, i, "|", color=INK, ms=11, mew=1.6)
    ax.plot([], [], "o", color=C[0], label="engine, Track A"); ax.plot([], [], "o", color=C[1], label="engine, Track B")
    ax.plot([], [], "|", color=INK, ms=10, mew=1.6, label="reference (design doc §2.E)")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows]); ax.invert_yaxis()
    ax.set_xlabel("$\\sigma_R$ [mb], p + $^{12}$C, $\\sigma_{NN}$ = 40 mb, no deflection")
    ax.set_title("QMC configurations: Track A vs Track B"); ax.legend(fontsize=7.5, loc="center")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_qmc_AvsB.png"), dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 4: Step 4(f)
ALGS = [("Veto", "Achilles Veto", INK2, "--"), ("Continuous", "Achilles Continuous", INK, "-"),
        ("B_k1_Atom", "B, $\\kappa$=1, atom (Veto law and profile)", C[4], ":"),
        ("B_k1_Box", "B, $\\kappa$=1, box (Continuous law and profile)", C[1], "-"),
        ("B_eik_Gauss", "B, eikonal, Gaussian (proposed)", C[2], "-"),
        ("A_rho", "A, $\\rho$", C[3], "-"), ("A_fold_finiteA", "A, $\\rho\\otimes\\hat p$ + finite-A", C[0], "-")]
def fig_achilles():
    r = json.load(open(CMPJ))
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.6))
    for j, (mode, key, lab) in enumerate((("CrossSection", "sigmaR", "$\\sigma_R$ [mb]"), ("Transparency", "T", "transparency"))):
        for a, name, col, ls in ALGS:
            res = r[f"NN_{mode}_{a}"]["result"]
            p = np.array([x["p"] for x in res]); y = np.array([x[key] for x in res]); e = np.array([x["err"] for x in res])
            ax[j].errorbar(p, y, yerr=e, color=col, ls=ls, marker="o", ms=3.5, label=name, capsize=0)
        ax[j].set_xlabel("proton momentum [MeV/c]"); ax[j].set_ylabel(lab)
    ax[0].set_title("(a) p + $^{12}$C reaction cross section"); ax[1].set_title("(b) $^{12}$C transparency (internal source)")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.12, 1, 1)); fig.savefig(os.path.join(FIG, "fig4_achilles_comparison.png"), dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 5: Step 4(h)
DALGS = [("Continuous", "Achilles Continuous", INK2, "--"),
         ("B_k1_Box", "B, $\\kappa$=1, box", C[1], "-"),
         ("B_eik_Gauss", "B, eikonal, Gaussian", C[2], "-"),
         ("A_fold_finiteA", "A, $\\rho\\otimes\\hat p$ + finite-A", C[0], "-"),
         ("A_rho", "A, $\\rho$", C[3], "-")]
def regroup(entry, edges_T, mode):
    lo, hi = np.array(entry["lo"]), np.array(entry["hi"])
    h, n = np.array(entry["hits"]), np.array(entry["nohits"])
    Tc = T_of_p(0.5*(lo+hi))
    out = []
    for a, b in zip(edges_T[:-1], edges_T[1:]):
        s = (Tc >= a) & (Tc < b)
        H, Nn = h[s].sum(), n[s].sum(); tot = H+Nn
        if tot < 50: out.append((np.nan, np.nan, np.nan)); continue
        f = H/tot; e = math.sqrt(f*(1-f)/tot)
        Tm = np.average(Tc[s], weights=(h+n)[s])
        if mode == "CrossSection":
            area = math.pi*7.0**2*10; out.append((Tm, area*f, area*e))
        else:
            out.append((Tm, 1-f, e))
    return np.array(out)

def load_csv(fn):
    rows = []
    for line in open(fn):
        if line.startswith("#") or line.startswith("Tp"): continue
        rows.append([float(x) for x in line.split(",")])
    return np.array(rows)

def fig_data():
    if DATAJ == "-" or not os.path.exists(DATAJ): return None
    r = json.load(open(DATAJ)); car = load_csv(CARL); tra = load_csv(TRANS)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    edges_low = [8, 12, 16, 20, 25, 30, 36, 43, 52, 62, 75, 90, 105]
    edges_high = [105, 130, 160, 200, 250, 310, 380, 460, 560, 680, 820, 980]
    chi = {}
    for a, name, col, ls in DALGS:
        lowk, highk = f"sigR_low_{a}", f"sigR_high_{a}"
        if lowk not in r or "lo" not in r[lowk] or "lo" not in r.get(highk, {}): continue
        g = np.vstack([regroup(r[lowk], edges_low, "CrossSection"), regroup(r[highk], edges_high, "CrossSection")])
        g = g[~np.isnan(g[:, 0])]
        ax[0].errorbar(g[:, 0], g[:, 1], yerr=g[:, 2], color=col, ls=ls, marker="o", ms=3, lw=1.5, label=name)
        # chi2/N against data (model interpolated in log T)
        model = np.interp(np.log(car[:, 0]), np.log(g[:, 0]), g[:, 1], left=np.nan, right=np.nan)
        merr = np.interp(np.log(car[:, 0]), np.log(g[:, 0]), g[:, 2])
        ok = ~np.isnan(model)
        for lab, sel in (("10-100 MeV", (car[:, 0] < 100)), ("100-860 MeV", (car[:, 0] >= 100))):
            m = ok & sel
            c2 = np.sum((car[m, 1]-model[m])**2/(car[m, 2]**2 + merr[m]**2)); chi[(a, "sigR", lab)] = (c2/m.sum(), m.sum(),
                np.mean(model[m]/car[m, 1]))
    ax[0].errorbar(car[:, 0], car[:, 1], yerr=car[:, 2], fmt="s", color=INK, ms=3, mfc="white", mew=0.9, lw=0.8,
                   label="data: Carlson (1996), digitised", zorder=5)
    ax[0].set_xscale("log"); ax[0].set_xlabel("proton kinetic energy $T_p$ [MeV]"); ax[0].set_ylabel("$\\sigma_R$ [mb]")
    ax[0].set_title("(a) p + $^{12}$C reaction cross section"); ax[0].set_ylim(0, 700)
    edges_t = [140, 220, 300, 380, 460, 560, 680, 820, 980, 1150, 1350, 1550, 1800, 1950]
    for a, name, col, ls in DALGS:
        k = f"transp_{a}"
        if k not in r or "lo" not in r[k]: continue
        g = regroup(r[k], edges_t, "Transparency"); g = g[~np.isnan(g[:, 0])]
        ax[1].errorbar(g[:, 0], g[:, 1], yerr=g[:, 2], color=col, ls=ls, marker="o", ms=3, lw=1.5, label=name)
        model = np.interp(tra[:, 0], g[:, 0], g[:, 1], left=np.nan, right=np.nan)
        merr = np.interp(tra[:, 0], g[:, 0], g[:, 2]); m = ~np.isnan(model)
        c2 = np.sum((tra[m, 1]-model[m])**2/(tra[m, 2]**2+merr[m]**2)); chi[(a, "T", "all")] = (c2/m.sum(), m.sum(), np.mean(model[m]/tra[m, 1]))
    ax[1].errorbar(tra[:, 0], tra[:, 1], yerr=tra[:, 2], fmt="s", color=INK, ms=3.5, mfc="white", mew=0.9, lw=0.8,
                   label="data: $^{12}$C(e,e$'$p), digitised", zorder=5)
    ax[1].set_xlabel("proton kinetic energy $T_p$ [MeV]  ($Q^2\\simeq 2m_NT_p$)"); ax[1].set_ylabel("transparency")
    ax[1].set_title("(b) $^{12}$C transparency (model: no-collision probability)"); ax[1].set_ylim(0.3, 1.0)
    ax[0].legend(fontsize=7, loc="upper right"); ax[1].legend(fontsize=7, loc="upper right")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig5_data.png"), dpi=200); plt.close(fig)
    with open(os.path.join(FIG, "chi2_data.txt"), "w") as f:
        for (a, obs, lab), (c2n, n, ratio) in chi.items():
            f.write(f"{a}\t{obs}\t{lab}\tchi2/N = {c2n:.2f}\tN = {n}\tmean model/data = {ratio:.3f}\n")
    return chi

fig_closed_forms(); fig_majorant_weights(); fig_qmc(); fig_achilles()
chi = fig_data()
if chi:
    for k, v in chi.items(): print(k, "chi2/N = %.2f  N = %d  model/data = %.3f" % v)
print("figures in", FIG)
