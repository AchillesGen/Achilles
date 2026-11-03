"""Exact toy-nucleus first-collision law (shared by analyse.py and plots.py)."""
import math
import numpy as np

sigf = 4.0

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

