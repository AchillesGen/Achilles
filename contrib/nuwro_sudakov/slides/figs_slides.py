"""Figures for the pedagogical slides.  Every input is either a NuWro run, a NuWro table
read from data/, or the NuWro formula re-coded (marked 'illustration' on the slide)."""
import numpy as np, re, json, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
sys.path.insert(0, '/home/claude/work/ana')
OUT = '/home/claude/slides/figs'; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 13, 'axes.grid': True, 'grid.alpha': 0.3, 'figure.dpi': 150,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'mathtext.fontset': 'cm'})
CA, CB, CX, CD = '#1f77b4', '#d62728', '#2ca02c', 'k'
GEN = '/home/claude/slides/gen'
M = 938.919; HC = 197.327
rng = np.random.default_rng(7)

def save(fig, name):
    fig.tight_layout(); fig.savefig(f'{OUT}/{name}.pdf'); plt.close(fig)

# ------------------------------------------------------------------ NuWro inputs
nuc = np.loadtxt(f'{GEN}/nucl.txt')          # r, rho, kF(p), Ef
r_, rho_, kf_, ef_ = nuc.T
corr = np.loadtxt(f'{GEN}/corr.txt')          # s, corr_ii, corr_ij, norm_ii, norm_ij

def table(fn, fields):
    rows = []; cur = {}
    for line in open(fn):
        line = line.strip()
        if line == '-':
            if cur: rows.append(cur); cur = {}
            continue
        m = re.match(r'(\w+):\s*([-\deE.+]+)', line)
        if m: cur[m.group(1)] = float(m.group(2))
    if cur: rows.append(cur)
    E = np.array([r['energy'] for r in rows])
    out = [E]
    for f in fields:
        y = np.array([r.get(f, np.nan) for r in rows]); m = ~np.isnan(y)
        out.append(np.interp(E, E[m], y[m]))
    return out
EX, SII, SIJ = table('/home/claude/src/nuwro/data/input/kaskada_NN_xsec_2.dat', ['xsec_ii', 'xsec_ij'])
EI, FII, FIJ = table('/home/claude/src/nuwro/data/input/kaskada_NN_inel_2.dat', ['inel_ii', 'inel_ij'])
sig_ii = lambda E: np.interp(E, EX, SII); sig_ij = lambda E: np.interp(E, EX, SIJ)

# ------------------------------------------------------------------ 1. density, kF, Ef
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot(r_, rho_, color=CA, lw=2.2)
ax[0].axhline(0.16, color='gray', ls=':', lw=1); ax[0].text(5.2, 0.163, r'$\rho_0=0.16$ fm$^{-3}$', color='gray', fontsize=11)
ax[0].set_xlabel('$r$ [fm]'); ax[0].set_ylabel(r'$\rho(r)$ [fm$^{-3}$]'); ax[0].set_xlim(0, 8.2)
ax[0].axvline(8.0004, color='k', ls='--', lw=1); ax[0].text(6.0, 0.12, 'cascade\nradius\n$R=8.0$ fm', fontsize=10)
ax[0].set_title(r'$^{12}$C point density (Fourier–Bessel)', fontsize=13)
ax[1].plot(r_, kf_, color=CB, lw=2.2, label=r'$k_F(r)$ [MeV]')
ax2 = ax[1].twinx(); ax2.plot(r_, ef_, color=CX, lw=2, ls='--', label=r'$E_F(r)$ [MeV]'); ax2.set_ylabel(r'$E_F$ [MeV]', color=CX)
ax2.grid(False)
ax[1].set_xlabel('$r$ [fm]'); ax[1].set_ylabel(r'$k_F^{\rm loc}(r)$ [MeV]', color=CB); ax[1].set_xlim(0, 8.2)
ax[1].set_title(r'local Fermi gas: $k_F=(3\pi^2\rho_p)^{1/3}$', fontsize=13)
save(fig, 'density')

# ------------------------------------------------------------------ 2. NN cross sections
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
E = np.logspace(-0.5, 3.5, 400)
ax[0].loglog(E, sig_ii(E), color=CA, lw=2, label=r'$pp,\ nn$  ($ii$)')
ax[0].loglog(E, sig_ij(E), color=CB, lw=2, label=r'$np$  ($ij$)')
ax[0].set_xlabel(r'$E_k^{\rm eff}$ [MeV]'); ax[0].set_ylabel(r'$\sigma_{\rm free}$ [mb]'); ax[0].legend()
ax[0].axvspan(0.1, 40, color='orange', alpha=0.12); ax[0].text(0.45, 25, 'slow secondaries\nlive here', fontsize=10)
ax[0].set_title('free NN cross section (PDG2016 table)', fontsize=12)
ax[1].plot(EI, FII, 'o-', color=CA, ms=3, label='$ii$'); ax[1].plot(EI, FIJ, 'o-', color=CB, ms=3, label='$ij$')
ax[1].set_xlim(0, 3000); ax[1].set_xlabel(r'$E_k^{\rm eff}$ [MeV]'); ax[1].set_ylabel(r'$f_{\rm inel}=\sigma_{\rm inel}/\sigma_{\rm tot}$'); ax[1].legend()
ax[1].axvline(280, color='gray', ls=':'); ax[1].text(300, 0.6, r'$f_{\rm inel}\equiv0$ below 280 MeV', fontsize=10)
ax[1].set_title('inelastic fraction (picks elastic / 1π / 2π)', fontsize=12)
save(fig, 'nnxs')

# ------------------------------------------------------------------ Ekeff machinery
def ekeff(Ek, q, c, m1=938.272, m2=939.565):
    p1 = np.sqrt(Ek**2 + 2 * m1 * Ek); E1 = Ek + m1; E2 = np.sqrt(q**2 + m2**2)
    return (E1 * E2 - p1 * q * c) / m2 - m1
def fermi_targets(kf, n):
    q = kf * rng.uniform(0, 1, n)**(1 / 3); c = rng.uniform(-1, 1, n); return q, c

# ------------------------------------------------------------------ 3. spread of sigma over the Fermi sea
fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), sharey=False)
kf0 = 263.0
for ax, Ek in zip(axs, [40, 200, 800]):
    q, c = fermi_targets(kf0, 200000); s = sig_ij(ekeff(Ek, q, c))
    ax.hist(s, bins=np.logspace(np.log10(max(s.min(), 1)), np.log10(s.max()), 60), color=CA, alpha=0.75, density=True)
    ax.set_xscale('log'); ax.axvline(s.mean(), color=CX, lw=2, label=f'mean {s.mean():.0f} mb')
    ax.axvline(s.max(), color=CB, lw=2, ls='--', label=f'max {s.max():.0f} mb')
    ax.set_title(f'$np$, projectile $E_k$ = {Ek} MeV', fontsize=12); ax.set_xlabel(r'$\sigma_{np}(E_k^{\rm eff}(T))$ [mb]')
    ax.legend(fontsize=10); ax.set_yticks([])
axs[0].set_ylabel('targets $T$ in the Fermi sea')
save(fig, 'sigma_spread')

# ------------------------------------------------------------------ 4. Pandharipande-Pieper factor and its bound
def mod(x, k1, k2):  # x = rho/rho0, k1,k2 vectors [MeV]
    beta = -116 * x; lam = (3.29 - 0.373 * x) * HC
    em = lambda p2: M / (1 - 2 * M * beta / lam**2 / (1 + p2 / lam**2)**2)
    p1s = (k1**2).sum(-1); p2s = (k2**2).sum(-1)
    m1, m2, m3 = em(p1s), em(p2s), em((p1s + p2s) / 2)
    return np.linalg.norm(k1 - k2, axis=-1) / M / np.linalg.norm(k1 / m1[..., None] - k2 / m2[..., None], axis=-1) * m3 / M
fig, ax = plt.subplots(2, 1, figsize=(6.2, 6.6))
kk = np.linspace(50, 2000, 120)
for x, col in [(0.5, '#9467bd'), (1.0, CA), (1.5, '#8c564b')]:
    mx = []; mn = []; av = []
    for k in kk:
        q, c = fermi_targets(kf0, 3000); s = np.sqrt(1 - c**2)
        k2 = np.stack([q * s, 0 * q, q * c], -1); k1 = np.broadcast_to([0, 0, k], k2.shape)
        m = mod(x, k1, k2); mx.append(m.max()); mn.append(m.min()); av.append(m.mean())
    ax[0].fill_between(kk, mn, mx, color=col, alpha=0.2); ax[0].plot(kk, av, color=col, lw=2, label=rf'$\rho={x}\rho_0$')
    lam = (3.29 - 0.373 * x) * HC; a = 2 * M * 116 * x / lam**2
    ax[0].axhline(1 / (1 - a / 4), color=col, ls='--', lw=1.3)
ax[0].set_xlabel(r'projectile $|\mathbf{k}_1|$ [MeV]'); ax[0].set_ylabel('mod'); ax[0].legend(fontsize=10)
ax[0].set_title('P&P factor over the Fermi sea (band); bound (dashed)', fontsize=11)
xs = np.linspace(0, 2, 100); lam = (3.29 - 0.373 * xs) * HC; a = 2 * M * 116 * xs / lam**2
ax[1].plot(xs, a, color=CA, lw=2, label=r'$a(\rho)=2M|\beta|/\lambda^2$'); ax[1].plot(xs, 1 / (1 - a / 4), color=CB, lw=2, label=r'bound $1/(1-a/4)$')
ax[1].axhline(4, color='gray', ls=':'); ax[1].set_xlabel(r'$\rho/\rho_0$'); ax[1].legend(fontsize=11); ax[1].set_ylim(0, 2.0)
ax[1].set_title('the bound needs $a<4$: always true here', fontsize=11)
save(fig, 'pp_bound')

# ------------------------------------------------------------------ 5. correlation tables
s, cii, cij, nii, nij = corr.T
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot(s, cii, color=CA, lw=2, label=r'corr$_{ii}(s)$'); ax[0].plot(s, cij, color=CB, lw=2, label=r'corr$_{ij}(s)$')
ax[0].set_xlabel('$s$ = distance travelled since last collision [fm]'); ax[0].legend(); ax[0].set_xlim(0, 8)
ax[0].set_title('correlation hole around the struck nucleon', fontsize=12)
ax[1].plot(s, nii, color=CA, lw=2, ls='--', label=r'norm$_{ii}(r)$'); ax[1].plot(s, nij, color=CB, lw=2, ls='--', label=r'norm$_{ij}(r)$')
ax[1].set_xlabel('$r$ = distance from the nuclear centre [fm]'); ax[1].legend(); ax[1].set_xlim(0, 8)
ax[1].set_title(r'normalisation (factor $=$ corr$(s)/$norm$(r)$)', fontsize=12)
save(fig, 'corr')

# ------------------------------------------------------------------ 6. survival + inverse transform
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
l = np.linspace(0, 6, 300); G = 0.6
ax[0].plot(l, np.exp(-G * l), color=CA, lw=2.5, label=r'$S(\ell)=e^{-\Gamma\ell}$')
ax[0].fill_between(l, 0, G * np.exp(-G * l), color=CB, alpha=0.25, label=r'$p(\ell)=\Gamma e^{-\Gamma\ell}$')
ax[0].set_xlabel(r'path length $\ell$ [fm]'); ax[0].legend(); ax[0].set_title(r'constant $\Gamma=0.6$ fm$^{-1}$', fontsize=12)
F = 1 - np.exp(-G * l); ax[1].plot(l, F, color=CA, lw=2.5, label=r'$F(\ell)=1-e^{-\Gamma\ell}$')
for u in [0.25, 0.7]:
    lam_ = -np.log(1 - u) / G
    ax[1].annotate('', xy=(lam_, u), xytext=(0, u), arrowprops=dict(arrowstyle='->', color=CB, lw=1.8))
    ax[1].annotate('', xy=(lam_, 0), xytext=(lam_, u), arrowprops=dict(arrowstyle='->', color=CB, lw=1.8))
    ax[1].text(0.1, u + 0.03, f'$u={u}$', color=CB)
ax[1].set_xlabel(r'$\ell$ [fm]'); ax[1].set_ylabel('CDF'); ax[1].set_title(r'inverse transform: $\lambda=-\ln(1-u)/\Gamma$', fontsize=12)
save(fig, 'survival')

# ------------------------------------------------------------------ 7. Jensen
fig, ax = plt.subplots(1, 2, figsize=(10, 3.7))
x = np.linspace(0, 2.2, 200); ax[0].plot(x, np.exp(-x), color=CA, lw=2.5, label=r'$e^{-x}$ (convex)')
x1, x2 = 0.2, 1.8; ax[0].plot([x1, x2], [np.exp(-x1), np.exp(-x2)], 'o-', color='gray', lw=1)
xm = (x1 + x2) / 2; ax[0].plot([xm], [(np.exp(-x1) + np.exp(-x2)) / 2], 's', color=CB, ms=9, label=r'$\mathbb{E}[e^{-\Gamma\Delta}]$')
ax[0].plot([xm], [np.exp(-xm)], 'D', color=CX, ms=9, label=r'$e^{-\langle\Gamma\rangle\Delta}$')
ax[0].set_xlabel(r'$x=\Gamma(T)\Delta$'); ax[0].legend(fontsize=11); ax[0].set_title('one random target per step survives more', fontsize=12)
# toy: Gamma = G(u1+u2) (two independent channels, uniform), per-step survival vs exact
R = 1.2 * 12**(1 / 3); rho = 1 / (4 / 3 * np.pi * 1.2**3); G = rho * 4.0
u = np.sqrt(rng.uniform(0, 1, 1000000)); b = 0.999 * R * u; z0 = -0.999 * R * np.sqrt(1 - u**2); L = np.sqrt(R**2 - b**2) - z0
Sx = lambda xx: ((1 - np.exp(-G * np.maximum(xx, 1e-12))) / (G * np.maximum(xx, 1e-12)))**2
exact = 1 - np.exp(-G * L).mean()
D = np.logspace(-2.2, 0.3, 40); pred = []
for d in D:
    n = np.floor(L / d); rr = L - n * d; pred.append(1 - (Sx(d)**n * Sx(rr)).mean() - exact)
ax[1].plot(D, np.array(pred), color=CA, lw=2, label='prediction for (a) [derived]')
uk = json.load(open('/home/claude/work/unit/ks_results.json'))
for k, d in [('a_s0.05_r1', 0.05), ('a_s0.2_r1', 0.2), ('a_s1.0_r1', 1.0)]:
    r = uk[k]; e = np.sqrt(r['Pint'] * (1 - r['Pint']) / r['n'])
    ax[1].errorbar([d], [r['Pint'] - r['Pint_exact']], [e], fmt='o', color=CB, ms=6, label='(a) NuWro runs' if d == 0.05 else None)
for k in ['b_seg1_r1', 'b_seg1_r1_h3']:
    r = uk[k]; e = np.sqrt(r['Pint'] * (1 - r['Pint']) / r['n'])
    ax[1].errorbar([0.012 if k.endswith('r1') else 0.015], [r['Pint'] - r['Pint_exact']], [e], fmt='D', color=CX, ms=6, label='(b) NuWro runs' if k.endswith('r1') else None)
ax[1].axhline(0, color='k', lw=0.8); ax[1].set_xscale('log'); ax[1].set_xlabel(r'step $\Delta$ [fm]'); ax[1].set_ylabel(r'$P_{\rm int}-P_{\rm int}^{\rm exact}$')
ax[1].legend(fontsize=10); ax[1].set_title(r'uniform sphere, $\sigma\sim U[0,2\bar\sigma]$ per draw', fontsize=12)
save(fig, 'jensen')

# ------------------------------------------------------------------ 8. the veto algorithm, pedagogical picture
def Gam(l): return 0.22 + 0.17 * np.sin(0.9 * l) * np.exp(-0.08 * l) + 0.2 * np.exp(-((l - 6.5) / 0.8)**2)
for seed in range(2000):
    rng2 = np.random.default_rng(seed); segs = np.arange(0, 10.01, 2.0); cands = []; first = None
    for a0, a1 in zip(segs[:-1], segs[1:]):
        gh = Gam(np.linspace(a0, a1, 200)).max() * 2.0; pos = a0
        while True:
            pos += rng2.exponential(1 / gh)
            if pos > a1: break
            y = rng2.uniform() * gh; acc = y < Gam(pos); cands.append((pos, y, acc))
            if acc and first is None: first = pos
    nrej = sum(1 for c in cands if first and c[0] < first)
    if first and 5.0 < first < 7.5 and nrej >= 3: break
fig, ax = plt.subplots(figsize=(11, 4.0))
l = np.linspace(0, 10, 500); ax.plot(l, Gam(l), color=CA, lw=2.5, label=r'true hazard $\Gamma(\ell)$ (never integrated)')
for a0, a1 in zip(segs[:-1], segs[1:]):
    gh = Gam(np.linspace(a0, a1, 200)).max() * 2.0
    ax.hlines(gh, a0, a1, color=CB, lw=2.2, label=r'majorant $\hat\Gamma$ (constant per segment)' if a0 == 0 else None)
    ax.vlines([a0, a1], 0, gh, color=CB, lw=0.6, ls=':')
k = 0
for pos, y, acc in cands:
    later = pos > first + 1e-9
    if not later: k += 1; ax.annotate(str(k), (pos, y), textcoords='offset points', xytext=(4, 6), fontsize=10)
    ax.plot(pos, y, 'o' if acc else 'x', color=(CX if acc else 'dimgray'), ms=9 if acc else 8, mew=2, alpha=0.25 if later else 1)
ax.plot([], [], 'x', color='dimgray', ms=8, mew=2, label=r'rejected: $u\hat\Gamma>\Gamma$ (go on from here)')
ax.plot([], [], 'o', color=CX, ms=9, label=r'accepted: $u\hat\Gamma<\Gamma$')
ax.axvline(first, color=CX, lw=1.5, ls='--'); ax.text(first + 0.08, 0.85, 'first accepted\n= the collision', color=CX, fontsize=11)
ax.text(first + 1.2, 0.05, 'later candidates are never drawn\n(shown faded)', fontsize=9, color='gray')
ax.set_xlabel(r'path length $\ell$ [fm]'); ax.set_ylabel(r'rate [fm$^{-1}$]'); ax.set_ylim(0, 1.25); ax.legend(loc="upper left", fontsize=10, ncol=2)
save(fig, 'veto_picture')

# ------------------------------------------------------------------ 9. hazard and majorant along a real chord through 12C (illustration with NuWro inputs)
Ek = 200.0; bimp = 1.0; Rc = 5.6; lmax = 2 * np.sqrt(Rc**2 - bimp**2); ll = np.linspace(0, lmax, 240)
rr = np.sqrt(bimp**2 + (ll - lmax / 2)**2)
rho_l = np.interp(rr, r_, rho_); kf_l = np.interp(rr, r_, kf_)
def avg_gamma(rh, kf, n=6000):
    if rh <= 0: return 0.0
    q, c = fermi_targets(max(kf, 1), n); E_ = ekeff(Ek, q, c); x = rh / 0.16
    s_ = np.sqrt(1 - c**2); k2 = np.stack([q * s_, 0 * q, q * c], -1)
    p1 = np.sqrt(Ek**2 + 2 * 938.272 * Ek); k1 = np.broadcast_to([0, 0, p1], k2.shape); md = mod(x, k1, k2)
    return rh / 2 * (np.mean(sig_ii(E_) * md) + np.mean(sig_ij(E_) * md)) * 0.1   # mb -> fm^2
Gl = np.array([avg_gamma(a, b) for a, b in zip(rho_l, kf_l)])
fig, ax = plt.subplots(figsize=(11, 3.9))
ax.plot(ll, Gl, color=CA, lw=2.5, label=r'$\Gamma_n+\Gamma_p=\langle\rho\,\sigma\,{\rm mod}\rangle_T$ along the chord')
for a0 in np.arange(0, lmax, 2.0):
    a1 = min(a0 + 2.0, lmax); m = (ll >= a0) & (ll <= a1)
    rhohi = rho_l[m].max(); kfhi = kf_l[m].max()
    p1 = np.sqrt(Ek**2 + 2 * 938.272 * Ek); E1 = Ek + 938.272; E2 = np.sqrt(kfhi**2 + 939.565**2)
    Emin = (E1 * E2 - p1 * kfhi) / 939.565 - 938.272; Emax = (E1 * E2 + p1 * kfhi) / 939.565 - 938.272
    grid = np.r_[Emin, Emax, EX[(EX > Emin) & (EX < Emax)]]; x = rhohi / 0.16; lam = (3.29 - 0.373 * x) * HC; a = 2 * M * 116 * x / lam**2
    gh = rhohi / 2 * (sig_ii(grid).max() + sig_ij(grid).max()) / (1 - a / 4) * 0.1
    ax.hlines(gh, a0, a1, color=CB, lw=2.2, label=r'$\hat\Gamma_n+\hat\Gamma_p$ on 2 fm segments' if a0 == 0 else None)
    ax.vlines([a0, a1], 0, gh, color=CB, lw=0.5, ls=':')
ax.set_xlabel(r'$\ell$ along a chord at $b=1$ fm [fm]'); ax.set_ylabel(r'fm$^{-1}$'); ax.legend(fontsize=11)
ax.set_title(r'200 MeV proton in $^{12}$C (illustration with NuWro inputs; SRC/mfp off as for the first collision)', fontsize=11)
save(fig, 'chord')

# ------------------------------------------------------------------ 10. probability integral transform (KS test)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
Rs = 1.2 * 12**(1 / 3); rho_s = 1 / (4 / 3 * np.pi * 1.2**3)
for axx, tag, col, ttl in [(ax[0], 'b_seg1_r1_h3', CX, r'(b): random $\sigma$, $\hat\Gamma\times3$'), (ax[1], 'a_s1.0_r1', CB, r'(a): random $\sigma$, $\Delta=1$ fm')]:
    a = np.loadtxt(f'/home/claude/work/unit/{tag}/sud_test.txt'); x0, y0, z0, c, xc, yc, zc = a.T
    Gs = rho_s * 4.0; L = np.sqrt(Rs**2 - x0**2 - y0**2) - z0; m = c == 1; d = zc[m] - z0[m]
    uu = (1 - np.exp(-Gs * d)) / (1 - np.exp(-Gs * L[m]))
    axx.hist(uu, bins=40, color=col, alpha=0.7, density=True); axx.axhline(1, color='k', lw=1.2)
    from scipy import stats; ks = stats.kstest(uu, 'uniform')
    axx.set_title(ttl + f'\nKS $D$={ks.statistic:.4f}, $p$={ks.pvalue:.2g}', fontsize=11); axx.set_xlabel(r'$u=F(d\,|\,L)$'); axx.set_ylim(0.85, 1.15)
save(fig, 'pit')

# ------------------------------------------------------------------ 11. why acceptance is low
fig, ax = plt.subplots(figsize=(7.5, 3.8))
Es = np.logspace(np.log10(15), np.log10(1500), 40)
for kf, ls in [(200, '--'), (263, '-')]:
    ratio_ij = []; ratio_ii = []
    for Ek_ in Es:
        q, c = fermi_targets(kf, 40000); E_ = ekeff(Ek_, q, c)
        p1 = np.sqrt(Ek_**2 + 2 * 938.272 * Ek_); E1 = Ek_ + 938.272; E2 = np.sqrt(kf**2 + 939.565**2)
        Emax = (E1 * E2 + p1 * kf) / 939.565 - 938.272; Emin = 0 if p1 <= kf else (E1 * E2 - p1 * kf) / 939.565 - 938.272
        g = np.r_[Emin, Emax, EX[(EX > Emin) & (EX < Emax)]]
        ratio_ij.append(sig_ij(E_).mean() / sig_ij(g).max()); ratio_ii.append(sig_ii(E_).mean() / sig_ii(g).max())
    ax.semilogx(Es, ratio_ij, color=CB, ls=ls, lw=2, label=f'$np$, $k_F$={kf} MeV'); ax.semilogx(Es, ratio_ii, color=CA, ls=ls, lw=2, label=f'$pp$, $k_F$={kf} MeV')
ax.set_xlabel('projectile $E_k$ [MeV]'); ax.set_ylabel(r'$\langle\sigma\rangle_T\,/\,\sup_T\sigma$'); ax.legend(fontsize=10); ax.set_ylim(0, 1)
ax.set_title('acceptance cost of bounding over the Fermi sea', fontsize=12)
save(fig, 'eff')
print('done')
