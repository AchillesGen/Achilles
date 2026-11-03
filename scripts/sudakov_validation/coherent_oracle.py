"""Usage: python3 coherent_oracle.py <density: r rho ...>
IPM eikonal oracle for the coherent correction, no Pauli.
Pure-phase (elastic-only, below-threshold) NN profile Gamma(s) = 1 - exp(i chi(s)),
chi = chi0 exp(-s^2/2w^2).  12C one-body thickness from the QMC density.
Per nucleon at fixed b (average over its transverse position, t normalised to 1):
  classical survival      c = 1 - <|G|^2>
  exact Glauber (IPM)     g = |1 - <G>|^2 = c + |<G>|^2
  F^2 rule (b-uniform)    u = 1 - <|G|^2> (1 - fbar),  fbar = int|<G>|^2 / int<|G|^2>
sigma_R = int d^2b [1 - s(b)^A]."""
import numpy as np, sys
d = np.loadtxt(sys.argv[1], comments="#"); r, rho = d[:, 0], d[:, 1]
# thickness per nucleon on a 2D grid
L, n = 12.0, 241
x = np.linspace(-L, L, n); dx = x[1] - x[0]
X, Y = np.meshgrid(x, x, indexing="ij"); B = np.hypot(X, Y)
z = np.linspace(-L, L, 481)
T = np.array([np.trapezoid(np.interp(np.hypot(bb, z), r, rho, right=0), z) for bb in x[x >= 0]])
bgrid = x[x >= 0]
t2 = np.interp(B, bgrid, T); t2 /= t2.sum() * dx * dx
A = 12
F = np.fft.fft2
def conv(a, b):  # (a * b)(b) on the grid
    return np.real(np.fft.fftshift(np.fft.ifft2(F(np.fft.ifftshift(a)) * F(np.fft.ifftshift(b))))) * dx * dx
def convc(a, b):
    return np.fft.fftshift(np.fft.ifft2(F(np.fft.ifftshift(a)) * F(np.fft.ifftshift(b)))) * dx * dx
def run(chi0, w):
    chi = chi0 * np.exp(-B**2 / (2 * w * w))
    G = 1 - np.exp(1j * chi)
    G2 = np.abs(G) ** 2
    sig = G2.sum() * dx * dx  # sigma_tot = sigma_el (fm^2)
    mG = convc(t2, G)        # <G>(b)
    mG2 = conv(t2, G2)       # <|G|^2>(b)
    coh = np.abs(mG) ** 2
    fbar = coh.sum() / mG2.sum()
    c = 1 - mG2; g = c + coh; u = 1 - mG2 * (1 - fbar)
    sR = lambda s: (1 - np.clip(s, 0, None) ** A).sum() * dx * dx * 10  # mb
    return sig * 10, fbar, sR(c), sR(g), sR(u), mG2.max()
print("chi0   w     sigNN[mb] fbar   sR_class  sR_Glauber  sR_F2rule  F2/Gl-1  Gl/cl-1  max<P>")
for chi0 in (np.pi / 3, 0.6):
    for w in (0.5, 0.8, 1.1, 1.5, 2.0, 2.6, 3.2):
        s, fb, c, g, u, pm = run(chi0, w)
        print(f"{chi0:.2f} {w:4.1f} {s:9.0f} {fb:6.3f} {c:9.1f} {g:10.1f} {u:10.1f} {u/g-1:+8.3f} {g/c-1:+8.3f} {pm:6.2f}")
