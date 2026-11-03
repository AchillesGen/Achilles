// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef CURVED_PATH_HH
#define CURVED_PATH_HH

#include <limits>
#include <memory>
#include <vector>

#include "Achilles/Potential.hh"
#include "Achilles/SymplecticIntegrator.hh"
#include "Achilles/ThreeVector.hh"

namespace achilles {

class Potential;

/// Pair collision kernel of the Achilles profiles in two forms:
///  - the per-pair probability P(b^2, sigma) of the closest-approach (Veto) trial, and
///  - the local hazard density kappa(r; sigma) [1/fm] in the pair separation r, defined as the
///    inverse Abel transform of X(b) = -ln[1 - P(b^2, sigma)], so that the integral of kappa
///    along any straight line with impact parameter b equals X(b) and the pair survives with
///    exactly 1 - P(b^2, sigma).
/// Gaussian: kappa = sigma^{-1/2} Li_{1/2}(exp(-pi r^2 / sigma)) ~ 1/r as r -> 0.
/// Pion:     kappa = (c/pi) sum_n K_0(n c r), c = sqrt(2 pi / sigma) ~ 1/(2r) as r -> 0.
/// Cylinder: X is infinite for b^2 < sigma/pi and zero outside, so kappa is a hard sphere of
///           radius sqrt(sigma/pi): a deterministic contact event, not a rate.
class PairKernel {
  public:
    enum class Profile { Gaussian, Pion, Cylinder };
    PairKernel() = default;
    explicit PairKernel(Profile profile) : m_profile{profile} {}

    Profile Type() const { return m_profile; }
    bool HardSphere() const { return m_profile == Profile::Cylinder; }
    /// Radius of the hard sphere (Cylinder only), sigma in fm^2
    static double ContactRadius(double sigma);

    /// P(b^2, sigma), sigma in fm^2
    double Probability(double b2, double sigma) const;
    /// X(b) = -ln(1 - P)
    double Exponent(double b2, double sigma) const;
    /// kappa(r; sigma) in 1/fm; decreasing in r, infinite at r = 0
    double Kappa(double r, double sigma) const;

    /// Li_{1/2}(e^{-t}) for t >= 0
    static double LiHalfExp(double t);
    /// sum_{n>=1} K_0(n x) for x > 0
    static double SumK0(double x);

  private:
    Profile m_profile{Profile::Gaussian};
};

/// Bicubic (Catmull-Rom) table of the real parts of a potential V(|p|, r), for fast and smooth
/// evaluation along trajectories. Radial nodes sit at (i + 1/2) dr and the table is extended
/// evenly to r < 0, so r = 0 (where the Schroedinger-equivalent potential is 0/0) is never
/// evaluated. Momenta above the table are clamped.
class TabulatedPotential : public Potential {
  public:
    TabulatedPotential(const Potential &base, double pmax = 2500, double dp = 5, double rmax = 12,
                       double dr = 0.02);
    PotentialVals operator()(double p, double r) const override;
    static std::string Name() { return "Tabulated"; }
    /// Value and analytic derivatives of the interpolant: {V, dV/dp, dV/dr} for the vector
    /// and the scalar part
    void Evaluate(double p, double r, double *vec, double *scal) const;

  private:
    double Interp(const std::vector<double> &table, double p, double r) const;
    void Interp3(const std::vector<double> &table, double p, double r, double *out) const;
    size_t m_np{}, m_nr{};
    double m_dp{}, m_dr{};
    std::vector<double> m_vector{}, m_scalar{};
};

/// Trajectory of one particle between collisions in a static, real mean field, with dense
/// output. The nodes come from Josh's SymplecticIntegrator (Tao's extended phase space,
/// Step<2>, omega = 20) on a fixed time grid t0 + k h; between nodes the position and the
/// momentum are cubic Hermite interpolants of (q, dq/dt) and (p, dp/dt). The cascade process
/// for this particle is defined on this dense path, which is a deterministic function of the
/// starting state. Without a potential the path is the exact straight line.
class CurvedPath {
  public:
    struct Node {
        double t;
        ThreeVector q, p, v, f;
    };
    /// Rigorous bounds for the dense path on [ta, tb] inside one segment:
    /// chord (a, b), maximal distance of the path from the chord (delta), and an upper bound
    /// on the speed |dq/dt| (vmax)
    struct Bound {
        ThreeVector a, b;
        double delta, vmax;
        double vmin{0}, amax{0}, ta{0}, tb{0};
    };

    CurvedPath() = default;
    CurvedPath(double t0, const ThreeVector &q0, const ThreeVector &p0, double mass,
               std::shared_ptr<Potential> potential, double step, double rexit, bool external);

    bool Straight() const { return !m_potential; }
    double T0() const { return m_t0; }
    double Mass() const { return m_mass; }
    double Step() const { return m_step; }
    /// Conserved energy: H(q0, p0) with the potential, sqrt(p^2 + m^2) without
    double Energy() const { return m_energy; }
    /// First time after t0 at which the particle is outside the sphere of radius rexit (and,
    /// for an external test particle, also past the plane z = rexit), as for Cascade::Escaped
    double ExitTime();
    /// True if the orbit can never reach the exit radius: bound behind the centrifugal barrier
    /// of the central field (classically trapped although H > m); such particles are captured
    bool Trapped();

    ThreeVector Position(double t);
    ThreeVector Velocity(double t);
    /// d^2q/dt^2 of the dense output (linear on each segment)
    ThreeVector Acceleration(double t);
    ThreeVector SegmentAcceleration(size_t k, double t);
    /// Kinetic momentum from the dense output
    ThreeVector Momentum(double t);
    /// Dense-output momentum rescaled onto the energy shell H(q(t), p) = Energy()
    ThreeVector ShellMomentum(double t);
    /// Hamiltonian of this particle at (q, p)
    double Hamiltonian(const ThreeVector &q, const ThreeVector &p) const;

    /// Index of the segment [t_k, t_{k+1}) containing t (the path is grown as needed)
    size_t Segment(double t);
    double SegmentStart(size_t k) { return Nodes(k + 1)[k].t; }
    double SegmentEnd(size_t k) { return Nodes(k + 2)[k + 1].t; }
    Bound Bounds(double ta, double tb);

    /// Largest |dq/dt| over the whole path up to tmax (rigorous, from the segment bounds)
    double MaxSpeed(double tmax);

    /// Bound for the block of 2^level segments starting at segment index << level, built from
    /// its children: the chord of the block, delta >= the distance of the path from that chord,
    /// vmax >= |dq/dt| >= vmin, amax >= |d^2q/dt^2|. Cached; shared by all pairs of the path.
    const Bound &Block(size_t level, size_t index);
    /// Index of the last segment that has to exist to reach the exit
    size_t ExitSegment();
    static constexpr size_t cMaxLevel = 14;

  private:
    const std::vector<Node> &Nodes(size_t n);
    void Grow();
    void Hermite(size_t k, double t, ThreeVector *q, ThreeVector *v) const;
    Node MakeNode(double t, const ThreeVector &q, const ThreeVector &p) const;
    static ThreeVector TableGradient(const TabulatedPotential &, const ThreeVector &q,
                                     const ThreeVector &p, double mass, bool radial);

    double m_t0{}, m_mass{}, m_step{}, m_rexit{}, m_energy{};
    bool m_external{false};
    double m_exit{-1};
    std::shared_ptr<Potential> m_potential{};
    std::shared_ptr<TabulatedPotential> m_table{};
    SymplecticIntegrator m_integrator{};
    std::vector<Node> m_nodes{};
    std::vector<std::vector<Bound>> m_blocks{};
    std::vector<std::vector<char>> m_block_set{};
};

} // namespace achilles

#endif
