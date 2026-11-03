"""Oracle for the Sudakov validation suite (Step 4 a-e).
Reads the CSVs written by run_all.sh and compares them with closed forms, exact enumeration
and the reference numbers of the design document.  Usage: python3 analyse.py <out dir>"""
import sys, os, math, itertools
import numpy as np
from scipy import stats, special, integrate

O = sys.argv[1]
A12, RMS = 12, 2.33
R = math.sqrt(5/3)*RMS
rho0 = A12/(4/3*math.pi*R**3)
results = []          # (section, test, value, expected, z or p, pass)

def load(name):
    return np.genfromtxt(os.path.join(O, name + '.csv'), delimiter=',', names=True)

def diag(name):
    s = open(os.path.join(O, name + '.diag')).read().split()
    return {k: float(v) for k, v in (t.split('=') for t in s)}

def check_z(sec, test, val, err, exp, zmax=3.0, exp_err=0.0):
    z = (val - exp)/math.sqrt(err**2 + exp_err**2) if err > 0 else 0.0
    results.append((sec, test, f"{val:.5f} ± {err:.5f}", f"{exp:.5f}", f"z = {z:+.2f}", abs(z) < zmax))

def check_p(sec, test, stat, p, desc):
    results.append((sec, test, desc, "", f"p = {p:.3f}", p > 0.01))

def chi2_counts(obs, exp, minexp=5):
    obs = np.asarray(obs, float); exp = np.asarray(exp, float)
    # merge bins with small expectation into the last one
    o, e = [], []
    co = ce = 0.0
    for oi, ei in zip(obs, exp):
        co += oi; ce += ei
        if ce >= minexp:
            o.append(co); e.append(ce); co = ce = 0.0
    if ce > 0:
        o[-1] += co; e[-1] += ce
    o, e = np.array(o), np.array(e)
    c2 = np.sum((o-e)**2/e); dof = len(o)-1
    return c2, dof, stats.chi2.sf(c2, dof)

# ---------------- closed forms (uniform sphere) ----------------
def ext_Pn(n, mu):
    y = 2*mu*R; return (n+1)*special.gammainc(n+2, y)/(2*mu**2*R**2)
def int_T(mu):
    y = 2*mu*R
    return 3/(4*R**3)*(R**2*(1-math.exp(-y))/mu - (2-math.exp(-y)*(y*y+2*y+2))/(4*mu**3))
def int_Pn(n, mu):
    y = 2*mu*R
    return 3/(4*R**3)*(R**2/mu*special.gammainc(n+1, y) - (n+2)*(n+1)/(4*mu**3)*special.gammainc(n+3, y))
mu40 = 4.0*rho0

def section_a():
    for name, sig in (('a_ext_A', 40), ('a_ext_A20', 20), ('a_ext_A30', 30)):
        d = load(name); mu = sig*0.1*rho0; N = len(d)
        pint = np.mean(d['nint'] > 0); pex = 1-ext_Pn(0, mu)
        check_z('a', f'A ext {sig} mb: P_int', pint, math.sqrt(pex*(1-pex)/N), pex)
        cnt = np.bincount(d['nint'].astype(int), minlength=40)[:40]
        c2, dof, p = chi2_counts(cnt, [N*ext_Pn(n, mu) for n in range(40)])
        check_p('a', f'A ext {sig} mb: N_int chi2', c2, p, f"chi2 = {c2:.1f}/{dof}")
        if sig == 40:
            dep = d['depth'][d['nint'] > 0]
            Lg = np.linspace(0, 2*R, 4001); fL = Lg/(2*R**2)
            cdf = lambda x: np.trapezoid(fL*(1-np.exp(-mu*np.minimum(x, Lg))), Lg)/pex
            ks = stats.kstest(dep, np.vectorize(cdf))
            check_p('a', 'A ext 40 mb: KS first-interaction depth', ks.statistic, ks.pvalue, f"D = {ks.statistic:.4f}")
    d = load('a_int_A'); N = len(d); T = np.mean(d['nint'] == 0); Tex = int_T(mu40)
    check_z('a', 'A int 40 mb: transparency', T, math.sqrt(Tex*(1-Tex)/N), Tex)
    cnt = np.bincount(d['nint'].astype(int), minlength=40)[:40]
    c2, dof, p = chi2_counts(cnt, [N*int_Pn(n, mu40) for n in range(40)])
    check_p('a', 'A int 40 mb: N_int chi2', c2, p, f"chi2 = {c2:.1f}/{dof}")
    check_z('a', 'A int 40 mb: <N>', d['nint'].mean(), d['nint'].std()/math.sqrt(N), 3*mu40*R/4)
    # hazard thinning a = 0.5
    d = load('a_int_A_pauli'); T = np.mean(d['nint'] == 0); Tex = int_T(0.5*mu40)
    check_z('a', 'A int, factor a=0.5: T = closed form at a*mu', T, math.sqrt(Tex*(1-Tex)/len(d)), Tex)

# Track B on Poisson configurations: N | b ~ Poisson(m(b)), m(b) = rho0 int d2s L(s) q(b - s)
def make_m(qfun, nb=161):
    bs = np.linspace(0, R+4, nb)
    # polar quadrature over s: s in [0, R], phi in [0, pi]
    sg, wsg = np.polynomial.legendre.leggauss(160); s = 0.5*R*(sg+1); ws = 0.5*R*wsg
    pg, wpg = np.polynomial.legendre.leggauss(160); ph = 0.5*math.pi*(pg+1); wp = 0.5*math.pi*wpg
    L = 2*np.sqrt(np.maximum(R**2-s**2, 0))
    out = []
    for b in bs:
        d2 = b*b + s[:, None]**2 - 2*b*s[:, None]*np.cos(ph[None, :])
        out.append(2*rho0*np.sum(ws[:, None]*wp[None, :]*s[:, None]*L[:, None]*qfun(d2)))
    return lambda b: np.interp(b, bs, out)
sigf = 4.0
def q_k1(d2): return np.exp(-math.pi*d2/sigf)
def Ein(x): return special.exp1(x)+np.log(x)+np.euler_gamma
x0 = 2.6; s2 = sigf/(2*math.pi*Ein(x0))
def q_eik(d2): return 1-np.exp(-x0*np.exp(-d2/(2*s2)))

def check_poisson(sec, test, d, mfun):
    b = d['b']; m = mfun(b); N = len(d)
    pint = np.mean(d['nint'] > 0); pex = np.mean(1-np.exp(-m))
    check_z(sec, test + ': P_int', pint, math.sqrt(pex*(1-pex)/N), pex)
    nmax = 30
    exp = [np.sum(stats.poisson.pmf(n, m)) for n in range(nmax)]
    cnt = np.bincount(d['nint'].astype(int), minlength=nmax)[:nmax]
    c2, dof, p = chi2_counts(cnt, exp)
    check_p(sec, test + ': N_int chi2', c2, p, f"chi2 = {c2:.1f}/{dof}")
    return pex

def section_a_B():
    m_k1 = make_m(q_k1); m_eik = make_m(q_eik)
    for prof in ('Atom', 'Box', 'Gauss'):
        check_poisson('a', f'B Poisson configs, kappa=1, {prof}', load(f'a_ext_B_k1_{prof}'), m_k1)
        check_poisson('a', f'B Poisson configs, eikonal x0=2.6, {prof}', load(f'a_ext_B_eik_{prof}'), m_eik)
    m_th = make_m(lambda d2: 1-np.sqrt(1-q_eik(d2)))
    check_poisson('a', 'B Poisson configs, factor a=0.5 -> 1-(1-p)^a per pair', load('a_ext_B_pauli'), m_th)
    return m_eik

# ---------------- (b) toy enumeration ----------------
TOYS = {2: [(0.3, 0, -0.4), (0.5, 0.2, 0.4)],
        3: [(0.3, 0, -1.0), (0.8, 0, 0), (0.2, 0.1, 0.5)],
        4: [(0, 0.3, -1.5), (0.6, 0, -0.3), (0.1, -0.4, 0.3), (0.9, 0.5, 1.2)]}
def toy_exact(pts, profile, t0=-np.inf):
    pts = np.array(pts, float); b2 = pts[:, 0]**2 + pts[:, 1]**2; zc = pts[:, 2]
    p = np.exp(-math.pi*b2/sigf); X = -np.log1p(-p)
    if profile == 'Atom':
        P = np.zeros(len(p)); surv = 1.0
        for j in np.argsort(zc):
            if zc[j] <= t0: continue
            P[j] = surv*p[j]; surv *= 1-p[j]
        return P, p
    # hazards in z (v cancels); box half-width sqrt(3) s, Gaussian rms s, s^2 = sigma/2pi
    s = math.sqrt(sigf/(2*math.pi)); l = math.sqrt(3)*s
    zg = np.linspace(-25, 25, 500001); dz = zg[1]-zg[0]; zm = zg[:-1]+dz/2
    H = []
    for j in range(len(p)):
        if profile == 'Box':
            h = np.where(np.abs(zm-zc[j]) < l, X[j]/(2*l), 0.0)
        else:
            phi = np.where(np.abs(zm-zc[j]) <= 5*s, np.exp(-(zm-zc[j])**2/(2*s*s))/(math.sqrt(2*math.pi)*s*math.erf(5/math.sqrt(2))), 0.0)
            h = X[j]*phi
        H.append(np.where(zm >= t0, h, 0.0))
    H = np.array(H); Ht = H.sum(0)
    S = np.exp(-np.concatenate([[0], np.cumsum(Ht*dz)]))
    frac = np.divide(H, Ht, out=np.zeros_like(H), where=Ht > 0)
    P = (frac*(S[:-1]-S[1:])).sum(1)
    pfull = 1-np.exp(-H.sum(1)*dz)
    return P, pfull

def section_b():
    for t, pts in TOYS.items():
        for prof in ('Atom', 'Box', 'Gauss'):
            d = load(f'b_toy{t}_{prof}'); N = len(d)
            P, pj = toy_exact(pts, prof)
            # first-collision index (including none = -1)
            obs = [np.sum(d['first'] == j) for j in range(len(pts))] + [np.sum(d['first'] == -1)]
            exp = list(N*P) + [N*(1-P.sum())]
            c2, dof, p = chi2_counts(obs, exp, minexp=1)
            check_p('b', f'toy A={t} {prof}: first-collision index chi2', c2, p, f"chi2 = {c2:.1f}/{dof}")
            # subset probabilities (product law, any profile)
            masks = np.arange(2**len(pts))
            expm = [N*np.prod([pj[i] if (m >> i) & 1 else 1-pj[i] for i in range(len(pts))]) for m in masks]
            obsm = [np.sum(d['mask'] == m) for m in masks]
            c2, dof, p = chi2_counts(obsm, expm, minexp=1)
            check_p('b', f'toy A={t} {prof}: collision-subset chi2', c2, p, f"chi2 = {c2:.1f}/{dof}")
    for prof in ('Atom', 'Box', 'Gauss'):
        d = load(f'b_toy3int_{prof}'); N = len(d)
        P, _ = toy_exact(TOYS[3], prof, t0=0.0)
        pesc = 1-P.sum()
        check_z('b', f'toy A=3 internal source {prof}: P_escape', np.mean(d['first'] == -1),
                math.sqrt(pesc*(1-pesc)/N), pesc)

# ---------------- (c) violated majorant ----------------
def section_c(m_eik):
    eff = 0.8*1.01
    d = load('c_ext_A_h08'); pint = np.mean(d['nint'] > 0)
    pb = 1-ext_Pn(0, eff*mu40); pt = 1-ext_Pn(0, mu40)
    check_z('c', 'A ext, majorant 0.8x: P_int equals closed form at min(G, Ghat)', pint, math.sqrt(pb*(1-pb)/len(d)), pb)
    dg = diag('c_ext_A_h08')
    results.append(('c', 'A ext, majorant 0.8x: violation counter', f"{int(dg['violations'])}/{int(dg['trials'])} trials", '> 0',
                    f"bias on sigma_R {100*(pb/pt-1):+.1f}%", dg['violations'] > 0))
    d = load('c_int_A_h08'); T = np.mean(d['nint'] == 0); Tb = int_T(eff*mu40)
    check_z('c', 'A int, majorant 0.8x: T equals closed form at min(G, Ghat)', T, math.sqrt(Tb*(1-Tb)/len(d)), Tb)
    # realistic: rho_hat = 0.16 on the QMC density
    for kind, ref in (('ext', 'e_ext_A_rho'), ('int', 'e_int_A_rho')):
        db, dr = load(f'c_{kind}_A_rho016'), load(ref)
        if kind == 'ext':
            fb, fr = np.mean(db['nint'] > 0), np.mean(dr['nint'] > 0); lab = 'sigma_R'
        else:
            fb, fr = np.mean(db['nint'] == 0), np.mean(dr['nint'] == 0); lab = 'T'
        err = math.sqrt(fb*(1-fb)/len(db) + fr*(1-fr)/len(dr))
        dg = diag(f'c_{kind}_A_rho016')
        results.append(('c', f'A {kind}, rho_hat = 0.16 fm^-3 (QMC peak 0.20): {lab} biased vs unbiased',
                        f"{fb:.4f} vs {fr:.4f} ({100*(fb/fr-1):+.2f}% ± {100*err/fr:.2f}%)", '',
                        f"violations {int(dg['violations'])}/{int(dg['trials'])}", dg['violations'] > 0))
    # Track B with headroom 0.7: per-pair hazard min(h, 0.7h) -> p_eff = 1-(1-p)^0.7
    m_b = make_m(lambda d2: 1-(1-q_eik(d2))**0.7)
    d = load('c_ext_B_h07')
    pint = np.mean(d['nint'] > 0); pb = np.mean(1-np.exp(-m_b(d['b']))); pt = np.mean(1-np.exp(-m_eik(d['b'])))
    check_z('c', 'B, majorant 0.7x: P_int equals prediction with p_eff = 1-(1-p)^0.7', pint, math.sqrt(pb*(1-pb)/len(d)), pb)
    dg = diag('c_ext_B_h07')
    results.append(('c', 'B, majorant 0.7x: violation counter', f"{int(dg['violations'])}/{int(dg['trials'])} trials", '> 0',
                    f"bias on P_int {100*(pb/pt-1):+.1f}%", dg['violations'] > 0))

# ---------------- (d) reweighting closure ----------------
def section_d(m_eik):
    Lg = np.linspace(0, 2*R, 4001); fL = 3/(4*R**3)*(R**2-Lg**2/4)
    for c in ('1.25', '2', '4'):
        d = load(f'd_int_A_c{c}'); N = len(d)
        # the majorant covers the largest variation (x1.2) and carries the 1% shell safety
        a = 1/(1.01*1.2*float(c))
        for k, f in ((0, 0.8), (1, 1.2)):
            w = d[f'w{k}']
            est = np.mean(w*(d['nint'] == 0)); err = np.std(w*(d['nint'] == 0))/math.sqrt(N)
            check_z('d', f'A int, headroom {c}: reweighted T at sigma x{f}', est, err, int_T(f*mu40))
            ess = np.sum(w)**2/np.sum(w**2)/N
            muhat = 1.01*1.2*float(c)*mu40; ap = f*a
            pred = 1/np.trapezoid(fL*np.exp(muhat*Lg*(ap-a)**2/(a*(1-a))), Lg)
            results.append(('d', f'A int, headroom {c}, sigma x{f}: ESS/N vs 1/E[W^2]', f"{ess:.3f}", f"{pred:.3f}",
                            f"diff {ess-pred:+.3f}", abs(ess-pred) < 0.01))
    d = load('d_ext_B_rw'); N = len(d)
    for k, (f, ref) in enumerate(((0.8, 'd_ext_B_32'), (1.2, 'd_ext_B_48'))):
        w = d[f'w{k}']
        est = np.mean(w*(d['nint'] > 0)); err = np.std(w*(d['nint'] > 0))/math.sqrt(N)
        dr = load(ref); pr = np.mean(dr['nint'] > 0); er = math.sqrt(pr*(1-pr)/len(dr))
        check_z('d', f'B (eikonal, Gauss) reweighted P_int at sigma x{f} vs regenerated', est, err, pr, exp_err=er)
        ess = np.sum(w)**2/np.sum(w**2)/N
        results.append(('d', f'B sigma x{f}: ESS/N', f"{ess:.3f}", '', '', True))
    d = load('d_int_A_pvar'); N = len(d)
    for k, a in ((0, 0.3), (1, 0.9)):
        w = d[f'w{k}']
        est = np.mean(w*(d['nint'] == 0)); err = np.std(w*(d['nint'] == 0))/math.sqrt(N)
        check_z('d', f'A int, factor 0.6 -> {a} by reweighting: T', est, err, int_T(a*mu40))

# ---------------- (e) A vs B on QMC ----------------
REF = {  # python reference (design doc 2.E; s5_qmc.py, s6_qmc_optical.py)
    'e_ext_A_rho': 229.6, 'e_ext_A_fold': 266.4, 'e_ext_A_foldFA': 271.9,
    'e_ext_B_k1_dec0': 275.6, 'e_ext_B_k1_dec1': 271.9, 'e_ext_B_eik_dec0': 270.0, 'e_ext_B_eik_dec1': 266.3,
    'e_int_B_Atom_dec0': 0.5164, 'e_int_B_Atom_dec1': 0.5106, 'e_int_B_Box_dec0': 0.4882,
    'e_int_B_Box_dec1': 0.4759, 'e_int_B_Gauss_dec0': 0.4889, 'e_int_B_Gauss_dec1': 0.4777,
    'e_int_A_rho': 0.4792}
LABEL = {'e_ext_A_rho': 'A on rho (optical)', 'e_ext_A_fold': 'A on rho (x) p (= exp(-k1))',
         'e_ext_A_foldFA': 'A on rho (x) p + finite-A factors (= B uncorrelated IPM)',
         'e_ext_B_k1_dec0': 'B kappa=1, QMC', 'e_ext_B_k1_dec1': 'B kappa=1, decorrelated',
         'e_ext_B_eik_dec0': 'B eikonal, QMC', 'e_ext_B_eik_dec1': 'B eikonal, decorrelated',
         'e_int_A_rho': 'A on rho (zero range)', 'e_int_A_foldFA': 'A on rho (x) p + finite-A'}
def section_e():
    bmax = 8.0
    rows = []
    for name in ('e_ext_A_rho', 'e_ext_A_fold', 'e_ext_A_foldFA', 'e_ext_B_k1_dec0', 'e_ext_B_k1_dec1',
                 'e_ext_B_eik_dec0', 'e_ext_B_eik_dec1'):
        d = load(name); N = len(d); p = np.mean(d['nint'] > 0)
        s = 10*math.pi*bmax**2*p; e = 10*math.pi*bmax**2*math.sqrt(p*(1-p)/N)
        check_z('e', f'sigma_R [mb] {LABEL[name]}', s, e, REF[name], exp_err=0.5)
        rows.append((name, s, e))
    for name in ('e_int_A_rho', 'e_int_A_foldFA') + tuple(f'e_int_B_{p}_dec{k}' for p in ('Atom', 'Box', 'Gauss') for k in (0, 1)):
        d = load(name); N = len(d); T = np.mean(d['nint'] == 0); e = math.sqrt(T*(1-T)/N)
        lab = LABEL.get(name, name.replace('e_int_B_', 'B eikonal ').replace('_dec0', ', QMC').replace('_dec1', ', decorrelated'))
        if name in REF:
            check_z('e', f'T {lab}', T, e, REF[name], exp_err=0.0015)
        else:
            results.append(('e', f'T {lab}', f"{T:.4f} ± {e:.4f}", '(new)', '', True))

section_a(); m_eik = section_a_B(); section_b(); section_c(m_eik); section_d(m_eik); section_e()
npass = sum(r[5] for r in results)
with open(os.path.join(O, 'REPORT.md'), 'w') as f:
    f.write(f"# Sudakov validation report\n\n{npass}/{len(results)} checks pass.\n\n")
    f.write("| step | test | result | expected | statistic | pass |\n|---|---|---|---|---|---|\n")
    for r in results:
        f.write(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {'PASS' if r[5] else 'FAIL'} |\n")
print(open(os.path.join(O, 'REPORT.md')).read())
