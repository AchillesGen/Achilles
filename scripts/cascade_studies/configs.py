"""Load Achilles nuclear configuration files and build comparison ensembles.

Ensembles used by the cascade studies:
  * ``load``     -- configurations exactly as Achilles uses them (weighted).
  * ``mixed``    -- every nucleon drawn from a different configuration. Same
                    1-body density, no correlations of any order (an
                    independent-particle sampler from rho(r)).
  * ``jastrow``  -- independent nucleons from the same rho(r) with a pair
                    weight prod g(r_ij) (Metropolis). Tests whether the
                    1-body density plus g(r) reproduces the full ensemble.
"""

import gzip
import pathlib

import numpy as np

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "configurations"


def load(name, max_configs=None):
    """Return positions (Nc, A, 3) [fm], is_proton (Nc, A) and weights (Nc,)."""
    path = DATA / name
    with gzip.open(path, "rt") as f:
        header = f.readline().split()
        nnuc, nconf = int(header[0]), int(header[1])
        if max_configs is not None:
            nconf = min(nconf, max_configs)
        pos = np.empty((nconf, nnuc, 3))
        isp = np.empty((nconf, nnuc), dtype=bool)
        wgt = np.empty(nconf)
        for c in range(nconf):
            for k in range(nnuc):
                tok = f.readline().split()
                isp[c, k] = tok[0] == "1"
                pos[c, k] = [float(x) for x in tok[1:4]]
            wgt[c] = float(f.readline())
            f.readline()
    return pos, isp, wgt


def unweight(pos, isp, wgt, rng):
    """Accept-reject on w/w_max, as DensityConfiguration::GetConfiguration does."""
    keep = rng.uniform(size=len(wgt)) < wgt / wgt.max()
    return pos[keep], isp[keep]


def mixed(pos, rng, nout=None, isp=None):
    """Independent-particle ensemble with the same 1-body density as ``pos``.

    With ``isp`` given, proton slots draw from the proton pool and neutron slots
    from the neutron pool, so rho_p and rho_n are kept separately.
    """
    nc, a, _ = pos.shape
    nout = nc if nout is None else nout
    if isp is None:
        flat = pos.reshape(-1, 3)
        return flat[rng.integers(0, flat.shape[0], size=(nout, a))]
    out = np.empty((nout, a, 3))
    for flag in (True, False):
        pool = pos[isp == flag]
        slots = np.flatnonzero(isp[0] == flag)
        out[:, slots] = pool[rng.integers(0, len(pool), size=(nout, len(slots)))]
    return out


def pair_distribution(pos, rmax=4.0, nbins=40):
    """Unnormalised histogram of all pair distances r_ij, i<j."""
    a = pos.shape[1]
    iu = np.triu_indices(a, 1)
    d = np.linalg.norm(pos[:, :, None, :] - pos[:, None, :, :], axis=-1)[:, iu[0], iu[1]]
    h, edges = np.histogram(d.ravel(), bins=nbins, range=(0, rmax))
    return h / pos.shape[0], edges


def g_of_r(pos, rng, rmax=4.0, nbins=40):
    """Pair correlation g(r) = rho2(r) / rho2_uncorrelated(r) (mixed-event method)."""
    h, edges = pair_distribution(pos, rmax, nbins)
    hm, _ = pair_distribution(mixed(pos, rng, 4 * pos.shape[0]), rmax, nbins)
    with np.errstate(invalid="ignore", divide="ignore"):
        return 0.5 * (edges[1:] + edges[:-1]), h / hm


def jastrow(pos_ref, g_r, g_val, rng, nout, nsweep=40, isp=None, niter=4):
    """Sample prod_i u(r_i) rho1(r_i) prod_{i<j} g(r_ij), with u tuned so the result keeps rho1.

    A move replaces one nucleon by a fresh draw from the 1-body pool of its
    isospin, picked with probability proportional to u(|r|) (an independence
    proposal), so the 1-body factor cancels in the Metropolis acceptance and
    only the pair weight g remains. The pair weight alone pushes nucleons
    outward; u is iterated (u <- u rho1_target / rho1_sampled, radial bins)
    until the sampled 1-body density matches the reference again.
    """
    if isp is None:
        isp = np.ones(pos_ref.shape[:2], dtype=bool)
    a = pos_ref.shape[1]
    slot_flag = isp[0]
    g_tab = np.clip(np.nan_to_num(g_val, nan=0.0), 1e-6, None)
    edges = np.linspace(0, 8, 41)

    def g(r):
        return np.interp(r, g_r, g_tab, left=g_tab[0], right=1.0)

    pools, pool_bin, u = {}, {}, {}
    for flag in (True, False):
        pools[flag] = pos_ref[isp == flag]
        pool_bin[flag] = np.clip(np.digitize(np.linalg.norm(pools[flag], axis=1), edges) - 1, 0, len(edges) - 2)
        u[flag] = np.ones(len(edges) - 1)
    target = {f: np.histogram(np.linalg.norm(pools[f], axis=1), edges)[0] + 1.0 for f in pools}

    def draw(flag, n):
        p = u[flag][pool_bin[flag]]
        return pools[flag][rng.choice(len(p), size=n, p=p / p.sum())]

    for _ in range(niter + 1):
        x = np.empty((nout, a, 3))
        for k in range(a):
            x[:, k] = draw(bool(slot_flag[k]), nout)
        for _ in range(nsweep):
            for k in range(a):
                new = draw(bool(slot_flag[k]), nout)
                others = np.delete(x, k, axis=1)
                r_old = np.linalg.norm(others - x[:, k, None, :], axis=-1)
                r_new = np.linalg.norm(others - new[:, None, :], axis=-1)
                ratio = np.prod(g(r_new) / g(r_old), axis=1)
                acc = rng.uniform(size=nout) < ratio
                x[acc, k] = new[acc]
        for flag in (True, False):
            got = np.histogram(np.linalg.norm(x[:, slot_flag == flag], axis=-1), edges)[0] + 1.0
            ratio = (target[flag] / target[flag].sum()) / (got / got.sum())
            u[flag] *= np.clip(ratio, 0.5, 2.0)
    return x


def write(path, pos, isp_slots):
    """Write an Achilles configuration file (unit weights)."""
    nc, a, _ = pos.shape
    with gzip.open(path, "wt") as f:
        f.write(f"{a} {nc} 1.0 1.0\n")
        for c in range(nc):
            for k in range(a):
                x, y, z = pos[c, k]
                f.write(f"{1 if isp_slots[k] else -1} {x:.6f} {y:.6f} {z:.6f}\n")
            f.write("1.0\n\n")
