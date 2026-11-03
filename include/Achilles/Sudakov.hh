// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Sudakov-veto intranuclear cascade.
//
// One sampling engine (global time ordering, piecewise-constant majorants, a veto chain of
// independent accept factors, variation weights) with two interchangeable rate models:
//   Track A  DensityRate        Gamma = sigma rho(r) v, target created at the collision point
//   Track B  ConfigurationRate  one hazard per (particle, configuration nucleon) pair
// and two collision models:
//   ForwardCollision   constant sigma, no deflection (validation against closed forms)
//   AchillesCollision  Achilles InteractionHandler: channels, 2->2 kinematics, Pauli factor
//
// Design and derivations: sudakov_inc_design.md (Steps 1-4).

#ifndef ACHILLES_SUDAKOV_HH
#define ACHILLES_SUDAKOV_HH

#include "Achilles/Particle.hh"

#include <array>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace YAML {
class Node;
}

namespace achilles {
class Event;
class Nucleus;
class InteractionHandler;
} // namespace achilles

namespace achilles::sudakov {

inline constexpr double kInf = std::numeric_limits<double>::infinity();
inline constexpr double kMbToFm2 = 0.1;
inline constexpr double kQcap = 0.9; ///< cap of the finite-A per-nucleon probability Q

/// Uniform random numbers in (0, 1)
using Uniform01 = std::function<double()>;

/// Piecewise-constant majorant in absolute time, plus optional point trials ("atoms").
/// Continuous segments must not overlap; they are summed by the rate models before Add.
class Majorant {
  public:
    struct Atom {
        double time;
        size_t label;
    };

    void Add(double t0, double t1, double rate);
    void AddAtom(double time, size_t label) { m_atoms.push_back({time, label}); }
    void SortAtoms();
    /// Next continuous trial strictly after t for uniform u: solves Lambda(t, T) = -ln u.
    double NextContinuous(double t, double u) const;
    /// First atom with time > t, or nullptr.
    const Atom *NextAtom(double t) const;
    double RateAt(double t) const;
    double Integral(double t0, double t1) const;
    bool Empty() const { return m_rate.empty() && m_atoms.empty(); }
    double End() const;
    size_t NSegments() const { return m_rate.size(); }

  private:
    std::vector<double> m_t0, m_t1, m_rate;
    std::vector<Atom> m_atoms;
};

/// One factor of the veto chain. `vars` holds the value under each variation
/// (empty: same as the nominal value for every variation).
struct Factor {
    std::string name;
    double value{1};
    std::vector<double> vars{};
};

/// Everything known about one trial.
struct Trial {
    size_t proj{};
    double time{};
    double majorant_rate{}; ///< value of the majorant at `time` (continuous trials)
    bool atom{false};
    size_t atom_label{};
    size_t pair{SIZE_MAX}; ///< Track B: index of the proposed pair in the line cache
    size_t n_rate_factors{
        0}; ///< factors [0, n) come from the rate model, the rest from the collision model
    size_t target{SIZE_MAX};
    bool created_target{false};
    ThreeVector point{};
    std::vector<Factor> factors{};
    std::vector<PID> out_ids{};
    std::vector<Particle> outgoing{};
    bool deflects{true};
};

/// Record of an accepted collision.
struct Collision {
    double time;
    size_t proj, target;
    ThreeVector point;
    PID proj_id, target_id;
};

/// Run-level diagnostics: trial and violation counters.
struct Diagnostics {
    size_t trials{}, accepted{}, violations{}, events{}, capped{};
    double max_ratio{};
    void Reset() { *this = Diagnostics{}; }
};

struct EngineSettings {
    size_t max_trials_per_event{2000000};
    bool clip_violations{true};
    /// Collision-model vetoes (Pauli, channel, kinematics) act as a single attempt per pair:
    /// a pair whose trial passes the rate factors but fails them is not tried again on this
    /// line (probability thinning, as Achilles Veto/Continuous). false: hazard thinning.
    /// MakeEngine defaults to true (YAML PauliAttempt: Single).
    bool single_attempt{false};
    std::vector<std::string> variations{};
};

/// Per-event state shared between the engine and the models.
class Context {
  public:
    double time{};
    bool new_line{true}; ///< set by the engine: the majorant is built for a new straight line
    size_t n_struck{};
    std::vector<Collision> collisions{};
    std::vector<ThreeVector> hole_centres{};
    std::vector<double> weights{};
    Diagnostics *diag{nullptr};
    Uniform01 uniform;
    size_t NVar() const { return weights.size(); }
};

/// Track-specific rate model.
class RateModel {
  public:
    virtual ~RateModel() = default;
    virtual std::string Name() const = 0;
    virtual void BeginEvent(Event &, Nucleus *, Context &) {}
    /// Majorant along particle idx's current line, from t0 (absolute time).
    virtual Majorant BuildMajorant(Event &, size_t idx, double t0, Context &) = 0;
    /// Choose the target and push the track factors. Return false if the trial cannot be
    /// a collision at all (e.g. the pair's nucleon is gone); that is a rejection with a = 0.
    virtual bool Propose(Event &, Trial &, Context &) = 0;
    virtual void OnAccept(Event &, const Trial &, Context &) {}
    virtual void OnReject(Event &, const Trial &, Context &) {}
    /// Called when the rate factors accepted the trial but a collision-model factor vetoed it
    /// and single-attempt mode is on.
    virtual void OnCollisionVeto(Event &, const Trial &, Context &) {}
    /// True if an accepted collision can raise other particles' true rates above their
    /// current majorants (then all majorants are rebuilt).
    virtual bool CommitInvalidatesOthers() const { return false; }
};

/// Shared NN physics: cross sections, channels, final states, Pauli factor.
class CollisionModel {
  public:
    virtual ~CollisionModel() = default;
    virtual std::string Name() const = 0;
    virtual void BeginEvent(Event &, Nucleus *) {}
    virtual bool Interacts(const Particle &) const = 0;
    /// Total cross section of the pair (mb). Indices into event.Hadrons().
    virtual double SigmaTot(Event &, size_t proj, size_t target) const = 0;
    /// Upper bound of SigmaTot over the targets the rate model can produce (mb).
    virtual double SigmaBound(Event &, const Particle &proj, double kmax) const = 0;
    /// Isospin-averaged SigmaTot with a target at rest (mb), for the finite-A void factor.
    virtual double SigmaAtRest(Event &e, const Particle &proj) const {
        return SigmaBound(e, proj, 0);
    }
    /// Multiplicative in-medium modification of SigmaTot for a pair at the given positions
    /// (1: none). Kept separate from SigmaTot so rate models can evaluate it where the
    /// collision would happen (Track B: closest approach) rather than where it is scheduled.
    virtual double Medium(const Particle &, const Particle &) const { return 1.0; }
    /// Choose the channel, generate the final state, push the Pauli (and kinematic) factors.
    virtual void Propose(Event &, Trial &, Context &) = 0;
    /// Apply an accepted collision. Returns indices of new particles that propagate.
    virtual std::vector<size_t> Commit(Event &, Trial &, Context &) = 0;
};

/// Smooth Pauli occupation f_tau(r, k) = Phi(k / kF_tau(r)), kF from the species density.
class Occupancy {
  public:
    enum class Shape { Step, FermiDirac };
    Occupancy() = default;
    Occupancy(Shape shape, double delta);
    void SetDensities(std::function<double(double)> rho_p, std::function<double(double)> rho_n);
    bool HasDensities() const { return static_cast<bool>(m_rho_p); }
    double KFermi(double r, PID pid) const; // MeV
    double Phi(double x) const;
    double f(double r, PID pid, double k) const;
    /// Sample |k| from 2 f(r, k) d^3k / rho(r) by inverting the tabulated shape.
    double SampleK(double r, PID pid, double u) const;
    double XMax() const { return m_xmax; }
    Shape GetShape() const { return m_shape; }

  private:
    Shape m_shape{Shape::Step};
    double m_delta{0.05}, m_norm{1}, m_xmax{1};
    std::vector<double> m_xgrid{}, m_cdf{};
    std::function<double(double)> m_rho_p{}, m_rho_n{};
    void Tabulate();
};

/// The Sudakov veto engine.
class Engine {
  public:
    Engine(std::shared_ptr<RateModel> rate, std::shared_ptr<CollisionModel> coll,
           EngineSettings settings, Uniform01 uniform);
    ~Engine();
    Engine(const Engine &) = delete;
    Engine &operator=(const Engine &) = delete;
    void Evolve(Event &, Nucleus *);
    const Context &GetContext() const { return m_ctx; }
    const Diagnostics &GetDiagnostics() const { return m_diag; }
    void ResetDiagnostics() { m_diag.Reset(); }
    RateModel &Rate() { return *m_rate; }
    CollisionModel &Coll() { return *m_coll; }
    const EngineSettings &Settings() const { return m_settings; }

  private:
    std::shared_ptr<RateModel> m_rate;
    std::shared_ptr<CollisionModel> m_coll;
    EngineSettings m_settings;
    Context m_ctx;
    Diagnostics m_diag;

    struct Entry {
        double time;
        size_t idx;
        size_t gen;
        bool atom;
        size_t atom_label;
        bool operator>(const Entry &o) const { return time > o.time; }
    };
    std::vector<size_t> m_gen{};
    std::vector<Majorant> m_major{};
    std::vector<Entry> m_heap{};

    void Schedule(Event &, size_t idx, double t0, bool rebuild, bool new_line);
    void Push(const Entry &);
    Entry Pop();
    void PropagateTo(Event &, double t);
};

// ------------------------------------------------------------------------------------------
// Track A: density backbone
// ------------------------------------------------------------------------------------------
struct DensityRateSettings {
    size_t shells{16};
    double rmax{-1};         ///< fm; < 0: nucleus radius (or the supplied value)
    double headroom{1.0};    ///< c >= 1 multiplies the majorant
    double fold_width{0};    ///< fm, rms per axis of a 3D Gaussian folding of rho
    bool finite_A{false};    ///< (A-m)/A depletion and 1/(1-Q) void factor (IPM-exact)
    double kmax_factor{1.6}; ///< target |k| <= kmax_factor * kF(0) (for the sigma bound)
    std::vector<double> hole_d{}, hole_g{}; ///< tabulated g(d); empty: no hole
    double hole_gmax{1};
    bool vertex_hole{true};
    double rho_hat_override{-1}; ///< > 0: use this as the density bound in every shell (tests)
    std::vector<double> sigma_scale_vars{}; ///< per variation: sigma -> s * sigma
};

class DensityRate : public RateModel {
  public:
    DensityRate(DensityRateSettings, std::shared_ptr<CollisionModel>, std::shared_ptr<Occupancy>);
    std::string Name() const override { return "SudakovA"; }
    /// Supply densities directly (standalone use); otherwise taken from the Nucleus.
    void SetDensities(std::function<double(double)> rho_p, std::function<double(double)> rho_n,
                      double A, double rmax);
    void BeginEvent(Event &, Nucleus *, Context &) override;
    Majorant BuildMajorant(Event &, size_t idx, double t0, Context &) override;
    bool Propose(Event &, Trial &, Context &) override;
    void OnAccept(Event &, const Trial &, Context &) override;
    void OnReject(Event &, const Trial &, Context &) override;
    bool CommitInvalidatesOthers() const override {
        return !m_s.hole_d.empty() && m_s.hole_gmax > 1;
    }
    double RhoTotal(double r) const;
    double KFermiMax() const { return m_kfmax; }

  private:
    DensityRateSettings m_s;
    std::shared_ptr<CollisionModel> m_coll;
    std::shared_ptr<Occupancy> m_occ;
    std::function<double(double)> m_rho_p{}, m_rho_n{};
    std::vector<double> m_r{}, m_rho_tab{}, m_rhop_tab{}, m_rhon_tab{};
    std::vector<double> m_shell_edges{}, m_shell_max{};
    double m_A{12}, m_rmax{-1}, m_kfmax{0};
    bool m_ready{false};
    // per-line cache
    struct Line {
        size_t gen_key{SIZE_MAX};
        double t0{}, sigma_hat{}, vhat{}, sigma_ref{}, qexit{};
        ThreeVector r0{}, beta{};
    };
    std::vector<Line> m_lines{};
    Nucleus *m_nuc{nullptr};
    void Tabulate();
    double Q(const Line &, double t) const;
    double HoleG(double d) const;
};

// ------------------------------------------------------------------------------------------
// Track B: configuration backbone
// ------------------------------------------------------------------------------------------
struct ConfigurationRateSettings {
    enum class Law { Kappa1, Eikonal };
    enum class Profile { Atom, Box, Gauss };
    Law law{Law::Eikonal};
    Profile profile{Profile::Gauss};
    double x0{2.6};
    double headroom{1.0};
    double gauss_nsigma{5.0};
    size_t gauss_bins{20};
    std::vector<double> sigma_scale_vars{};
};

class ConfigurationRate : public RateModel {
  public:
    ConfigurationRate(ConfigurationRateSettings, std::shared_ptr<CollisionModel>);
    std::string Name() const override { return "SudakovB"; }
    Majorant BuildMajorant(Event &, size_t idx, double t0, Context &) override;
    bool Propose(Event &, Trial &, Context &) override;
    void OnCollisionVeto(Event &, const Trial &, Context &) override;
    /// Per-pair probability for a full line (public for tests).
    double PairProb(double b2, double sigma_fm2) const;
    double PairX(double b2, double sigma_fm2) const;
    double Width2(double sigma_fm2) const; ///< transverse s^2 (fm^2)

  private:
    ConfigurationRateSettings m_s;
    std::shared_ptr<CollisionModel> m_coll;
    struct Pair {
        size_t j;
        double tc, b2, v;
        std::vector<double> sigmas; ///< nominal, then each variation (fm^2)
        double hat_lo, hat_hi;      ///< support of the bound
        bool blocked{false};        ///< single-attempt mode: vetoed once on this line
    };
    struct LineCache {
        size_t key{SIZE_MAX};
        std::vector<Pair> pairs{};
    };
    std::vector<LineCache> m_cache{};
    double PairHazard(const Pair &, double sigma, double t) const;
    void PairBoundPieces(const Pair &, double t0, std::vector<std::array<double, 3>> &) const;
    double PairBound(const Pair &, double t) const;
};

// ------------------------------------------------------------------------------------------
// Collision models
// ------------------------------------------------------------------------------------------
/// Constant sigma, no deflection, optional constant Pauli-like factor (validation).
class ForwardCollision : public CollisionModel {
  public:
    ForwardCollision(double sigma_mb, double pauli_a = 1.0, std::vector<double> pauli_vars = {})
        : m_sigma(sigma_mb), m_a(pauli_a), m_avars(std::move(pauli_vars)) {}
    std::string Name() const override { return "Forward"; }
    bool Interacts(const Particle &p) const override { return p.Info().IsNucleon(); }
    double SigmaTot(Event &, size_t, size_t) const override { return m_sigma; }
    double SigmaBound(Event &, const Particle &, double) const override { return m_sigma; }
    void Propose(Event &, Trial &, Context &) override;
    std::vector<size_t> Commit(Event &, Trial &, Context &) override;

  private:
    double m_sigma, m_a;
    std::vector<double> m_avars;
};

/// Achilles InteractionHandler: channels, 2->2 kinematics, Pauli factor (1-f3)(1-f4).
struct AchillesCollisionSettings {
    /// sigma_cap (mb, <= 0: none) caps the NN cross section. The Achilles NN elastic
    /// parametrisations diverge as 1/p_cm^2 at zero relative momentum, which makes any
    /// majorant over Fermi-moving targets unbounded; the cap is a stated model choice.
    double sigma_cap{-1};
    /// Pauli blocking of outgoing nucleons (false: Fermi motion kept, blocking off)
    bool pauli{true};
    /// Coherent veto: an elastic NN collision whose final state passes the Pauli check is
    /// vetoed with probability F_t(q)^2 (q: momentum transfer to the nucleus, F_t: one-body
    /// form factor of the struck species), i.e. the nucleus stays in its ground state.
    bool coherent{false};
    /// If non-empty, the coherent veto uses this table instead of F_t(q)^2: veto probability
    /// v(q) (q in MeV, linear interpolation, constant beyond the ends) for Pauli-allowed
    /// elastic collisions. Used for the IPM correction v = 1 - S_inc(q) / L_LFG(q), which
    /// removes only what the local-Fermi-gas Pauli factor misses (scripts/.../ipm_veto.py).
    std::vector<double> veto_q{}, veto_v{};
    /// Bound on the in-medium factor used in the Track A majorant (violations are counted)
    double medium_max{1.0};
};

class AchillesCollision : public CollisionModel {
  public:
    AchillesCollision(std::shared_ptr<Occupancy> occ, double sigma_cap = -1)
        : m_occ(std::move(occ)) {
        m_s.sigma_cap = sigma_cap;
    }
    AchillesCollision(std::shared_ptr<Occupancy> occ, AchillesCollisionSettings s)
        : m_occ(std::move(occ)), m_s(s) {}
    void Bind(InteractionHandler *handler) { m_handler = handler; }
    void SetMedium(std::function<double(const Particle &, const Particle &)> f) {
        m_medium = std::move(f);
    }
    /// Tabulate F_p, F_n from the species densities (any normalisation).
    void SetFormFactors(const std::function<double(double)> &rho_p,
                        const std::function<double(double)> &rho_n, double rmax = 15.0);
    bool HasFormFactors() const { return !m_ff_p.empty(); }
    /// One-body form factor of the species, q in MeV; F(0) = 1.
    double FormFactor(PID pid, double q) const;
    /// Coherent veto probability for momentum transfer q (MeV) on a struck nucleon of pid.
    double CoherentVetoProb(PID pid, double q) const;
    const AchillesCollisionSettings &Settings() const { return m_s; }
    std::string Name() const override { return "Achilles"; }
    bool Interacts(const Particle &p) const override { return p.Info().IsNucleon(); }
    double SigmaTot(Event &, size_t, size_t) const override;
    double SigmaBound(Event &, const Particle &, double kmax) const override;
    double SigmaAtRest(Event &, const Particle &) const override;
    double Medium(const Particle &a, const Particle &b) const override {
        return m_medium ? m_medium(a, b) : 1.0;
    }
    void Propose(Event &, Trial &, Context &) override;
    std::vector<size_t> Commit(Event &, Trial &, Context &) override;
    void BeginEvent(Event &, Nucleus *) override;

  private:
    std::shared_ptr<Occupancy> m_occ;
    AchillesCollisionSettings m_s{};
    InteractionHandler *m_handler{nullptr};
    std::function<double(const Particle &, const Particle &)> m_medium{};
    // bounds are cached per (projectile PID, |p| in 1 MeV bins, kmax in 1 MeV bins)
    mutable std::map<std::tuple<long, long, long>, double> m_bound_cache{}, m_rest_cache{};
    double SigmaBoundUncached(Event &, const Particle &, double kmax) const;
    double SigmaAtRestUncached(Event &, const Particle &) const;
    static constexpr double kFFdq = 2.0; // MeV
    std::vector<double> m_ff_p{}, m_ff_n{};
};

// ------------------------------------------------------------------------------------------
// Utilities (public for tests)
// ------------------------------------------------------------------------------------------
double Ein(double x); ///< Ein(x) = int_0^x (1 - e^-t)/t dt

/// Build an engine from a YAML node (Cascade/Sudakov). `track` is "A" or "B".
std::shared_ptr<Engine> MakeEngine(const YAML::Node &, char track, Uniform01 uniform);

} // namespace achilles::sudakov

#endif
