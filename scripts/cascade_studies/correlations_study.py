"""How much do n-body correlations in the nuclear configurations matter for the cascade?

Frozen spectators, straight lines, Achilles' linear per-pair survival 1 - P(b)
(no Pauli, a = 1). Observables:
  T(sigma)   survival of a nucleon leaving a vertex on a nucleon of the configuration
  sigma_R    external-beam reaction cross section, all A nucleons as targets
Ensembles: QMC (full correlations), Jastrow (rho1 + g(r) only), mixed (rho1 only).
Models built on rho1 alone: the optical / mean-free-path limit exp(-sigma int rho dz).

Run: python3 scripts/cascade_studies/correlations_study.py [outdir]
"""

import json
import pathlib
import sys

import numpy as np

import configs

RNG = np.random.default_rng(20261010)
SIGMAS_MB = [10, 20, 30, 40, 50]
PROFILES = {
    "Gaussian": lambda b2, s: np.exp(-np.pi * b2 / s),
    "Cylinder": lambda b2, s: (b2 < s / np.pi).astype(float),
}


def random_dirs(n):
    v = RNG.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def vertex_survival(pos, sigma_fm2, prof):
    """Per-event survival S, F1 = sum P_j, and sum P_j^2 for every nucleon as vertex."""
    nc, a, _ = pos.shape
    n = random_dirs(nc * a).reshape(nc, a, 3)
    d = pos[:, None, :, :] - pos[:, :, None, :]  # d[c, i, j] = r_j - r_i
    z = np.einsum("cijk,cik->cij", d, n)
    b2 = np.einsum("cijk,cijk->cij", d, d) - z * z
    p = PROFILES[prof](b2, sigma_fm2) * (z > 0)
    idx = np.arange(a)
    p[:, idx, idx] = 0.0
    s = np.prod(1.0 - p, axis=2).ravel()
    f1 = p.sum(axis=2).ravel()
    f2 = (p.sum(axis=2) ** 2 - (p * p).sum(axis=2)).ravel()
    return s, f1, f2


def sigma_r(pos, sigma_fm2, prof, bmax=7.0, nrep=40):
    """sigma_R [mb] for a beam along a random axis, impact parameter uniform in a disk."""
    pos = np.tile(pos, (nrep, 1, 1))
    nc = pos.shape[0]
    n = random_dirs(nc)
    # orthonormal transverse basis
    t1 = np.cross(n, np.where(np.abs(n[:, :1]) < 0.9, [[1, 0, 0]], [[0, 1, 0]]))
    t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    t2 = np.cross(n, t1)
    rb = bmax * np.sqrt(RNG.uniform(size=nc))
    phi = RNG.uniform(0, 2 * np.pi, size=nc)
    bvec = rb[:, None] * (np.cos(phi)[:, None] * t1 + np.sin(phi)[:, None] * t2)
    d = pos - bvec[:, None, :]
    z = np.einsum("cjk,ck->cj", d, n)
    b2 = np.einsum("cjk,cjk->cj", d, d) - z * z
    s = np.prod(1.0 - PROFILES[prof](b2, sigma_fm2), axis=1)
    area = np.pi * bmax**2 * 10.0  # mb
    return area * (1 - s.mean()), area * s.std() / np.sqrt(nc)


def optical_vertex(pos, sigma_fm2, a, nev=200000, smax=12.0, ns=240):
    """Mean-free-path transparency from rho1 alone: <exp(-sigma (A-1)/A int_0^inf rho dz)>.

    rho(r) is the angle-averaged 1-body density of the ensemble; the vertex is
    drawn from rho and the direction is isotropic. This is what an MFP cascade
    on a density (no configurations) computes.
    """
    rr = np.linalg.norm(pos.reshape(-1, 3), axis=1)
    edges = np.linspace(0, 8, 161)
    h, _ = np.histogram(rr, bins=edges)
    shell = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    rho = a * h / len(rr) / shell
    rc = 0.5 * (edges[1:] + edges[:-1])
    flat = pos.reshape(-1, 3)
    r0 = flat[RNG.integers(0, len(flat), size=nev)]
    n = random_dirs(nev)
    s = np.linspace(0, smax, ns)
    pts = r0[:, None, :] + s[None, :, None] * n[:, None, :]
    col = np.trapezoid(np.interp(np.linalg.norm(pts, axis=-1), rc, rho, right=0.0), s, axis=1)
    return np.exp(-sigma_fm2 * col * (a - 1) / a).mean()


def main(outdir):
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    results = {}

    sets = {
        "12C": ("QMC_configs.out.gz", "MF_configs.out.gz"),
        "16O": ("16O_AFDMC_configs.out.gz", "16O_MF_configs.out.gz"),
    }
    for nuc, (fq, fmf) in sets.items():
        pos, isp, w = configs.load(fq, max_configs=36000)
        pos, _ = configs.unweight(pos, isp, w, RNG)
        mf, mfp, mfw = configs.load(fmf, max_configs=36000)
        mf, _ = configs.unweight(mf, mfp, mfw, RNG)
        a = pos.shape[1]
        nev = pos.shape[0]
        mix = configs.mixed(pos, RNG, nev)
        r, g = configs.g_of_r(pos, RNG)
        _, g_mf = configs.g_of_r(mf, RNG)
        jas = configs.jastrow(pos, r, g, RNG, nev)
        _, g_jas = configs.g_of_r(jas, RNG)
        ens = {"QMC": pos, "Jastrow": jas, "Mixed": mix, "MF": mf}
        res = {
            "n_configs": nev,
            "g_r": {"r": r.tolist(), "QMC": g.tolist(), "MF": g_mf.tolist(), "Jastrow": g_jas.tolist()},
            "rms_radius": {k: float(np.sqrt((v**2).sum(-1).mean())) for k, v in ens.items()},
        }
        for prof in PROFILES:
            for sig in SIGMAS_MB:
                sf = sig / 10.0
                row = {}
                for name, x in ens.items():
                    s, f1, f2 = vertex_survival(x, sf, prof)
                    m1, m2 = f1.mean(), f2.mean()
                    row[name] = {
                        "T": float(s.mean()),
                        "T_err": float(s.std() / np.sqrt(len(s))),
                        # first factorial cumulant: optical with the exact conditional density
                        "T_F1": float(np.exp(-m1)),
                        # second order: adds spectator-spectator (and 3-body) correlations
                        "T_F2": float(np.exp(-m1 + 0.5 * (m2 - m1**2))),
                        "sigmaR": sigma_r(x, sf, prof),
                    }
                row["Optical_rho1"] = {"T": float(optical_vertex(pos, sf, a))}
                res[f"{prof}_{sig}"] = row
                print(nuc, prof, sig, {k: round(v["T"], 4) for k, v in row.items()},
                      {k: round(v["sigmaR"][0], 1) for k, v in row.items() if "sigmaR" in v}, flush=True)
        results[nuc] = res
    (outdir / "correlations.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "cascade_study_out")
