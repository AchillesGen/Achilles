"""IPM coherent veto table for the Sudakov collision model (YAML key Sudakov/CoherentTable).

For momentum transfer q to the nucleus, an exact independent-particle (Slater determinant)
nucleus leaves the nucleus excited with probability (closure, static)
    S_inc(q) = 1 - Tr(rho e rho e^dag) / N,   e = exp(i q.r),
which contains both Pauli blocking and the removal of the coherent (ground-state) piece.
The cascade's local-Fermi-gas Pauli factor gives instead, on average over the struck nucleon,
    L(q) = < 3x/4 - x^3/16 >_rho,   x = q / kF(r)  (1 for x >= 2).
Vetoing Pauli-allowed elastic collisions with v(q) = max(0, 1 - S_inc(q)/L(q)) makes the
average match S_inc(q). The F^2 rule instead vetoes with F(q)^2, which double counts the
coherent piece that L(q) already removes.

Nucleus: 12C as (0s1/2)^2 (0p3/2)^4 per isospin, HO length b fitted to the rms of the density.
Usage: python3 ipm_veto.py <density file: r rho_p [rho_n]> [qmax MeV] [dq MeV]  (prints YAML)"""
import sys
import numpy as np

HBARC = 197.327


def p32_projector():
    I2 = np.eye(2)
    S = [np.array([[0, 1], [1, 0]]) / 2, np.array([[0, -1j], [1j, 0]]) / 2, np.diag([0.5, -0.5])]
    eps = np.zeros((3, 3, 3))
    for i, j, k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        eps[i, j, k], eps[i, k, j] = 1, -1
    LS = sum(np.kron(-1j * eps[k], S[k]) for k in range(3))
    P = (2 / 3) * (np.eye(6) + LS)
    assert np.allclose(P @ P, P)
    return I2, P


def s_inc(q, b):
    """q in fm^-1; returns S_inc and F for (0s1/2)^2 (0p3/2)^4."""
    I2, P = p32_projector()
    y = (q * b) ** 2 / 4
    ey = np.exp(-y)
    E = np.zeros((4, 4), complex)  # basis s, px, py, pz; e = exp(i q z)
    E[0, 0] = E[1, 1] = E[2, 2] = ey
    E[3, 3] = (1 - 2 * y) * ey
    E[0, 3] = E[3, 0] = 1j * (q * b / np.sqrt(2)) * ey
    Ef = np.kron(E, I2)
    rho = np.zeros((8, 8), complex)
    rho[:2, :2] = I2
    rho[2:, 2:] = P
    N = np.trace(rho).real
    return 1 - np.trace(rho @ Ef @ rho @ Ef.conj().T).real / N, np.trace(rho @ Ef).real / N


def lfg_allowed(q, r, rho):
    kf = np.cbrt(3 * np.pi**2 * rho)
    x = np.where(kf > 0, q / np.maximum(kf, 1e-12), 1e9)
    allow = np.where(x < 2, 0.75 * x - x**3 / 16, 1.0)
    w = rho * r**2
    return np.trapezoid(allow * w, r) / np.trapezoid(w, r)


def table(r, rho, qmax=800.0, dq=10.0):
    rms = np.sqrt(np.trapezoid(rho * r**4, r) / np.trapezoid(rho * r**2, r))
    b = rms / np.sqrt(26 / 12)
    qs = np.arange(0.0, qmax + 0.5 * dq, dq)
    vs = []
    for qm in qs:
        q = max(qm, 1e-3) / HBARC
        S, _ = s_inc(q, b)
        L = lfg_allowed(q, r, rho)
        vs.append(float(np.clip(1 - S / L, 0, 1)) if L > 0 else 1.0)
    return qs, np.array(vs), b, rms


if __name__ == "__main__":
    d = np.loadtxt(sys.argv[1], comments="#")
    qmax = float(sys.argv[2]) if len(sys.argv) > 2 else 800.0
    dq = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
    qs, vs, b, rms = table(d[:, 0], d[:, 1], qmax, dq)
    print(f"    # IPM coherent veto, 12C (0s1/2)^2(0p3/2)^4, HO b = {b:.3f} fm (rms {rms:.3f} fm)")
    print("    CoherentTable:")
    print("      q: [" + ", ".join(f"{q:g}" for q in qs) + "]")
    print("      v: [" + ", ".join(f"{v:.4f}" for v in vs) + "]")
