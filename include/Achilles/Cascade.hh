// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef CASCADE_HH
#define CASCADE_HH

#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <optional>
#include <queue>
#include <vector>

#include "Achilles/CurvedPath.hh"
#include "Achilles/DecayHandler.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/InteractionHandler.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/SymplecticIntegrator.hh"
#include "Achilles/ThreeVector.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic pop

namespace achilles {

class Nucleus;
class Particle;
class Event;
class PID;

using Particles = std::vector<Particle>;
using InteractionDistances = std::vector<std::pair<std::size_t, double>>;

/// The Cascade class performs a cascade of the nucleons inside the nucleus. The nucleons that
/// are struck in the hard interaction propagate through the nuclear medium. To determine if an
/// interaction occurs, we calculate the interaction cross-section of Np and Nn, where N is the
/// propagating nucleon.
class Cascade {
    static constexpr int cMaxSteps = 100000;

  public:
    // Probability Enums
    enum ProbabilityType { Gaussian, Pion, Cylinder };
    std::string ToString(const ProbabilityType &type) {
        switch(type) {
        case Gaussian:
            return "Gaussian";
        case Pion:
            return "Pion";
        case Cylinder:
            return "Cylinder";
        }
        return "Unknown";
    }

    // In-Medium Effects Enums
    enum InMedium { None, NonRelativistic, Relativistic };
    std::string ToString(const InMedium &type) {
        switch(type) {
        case None:
            return "None";
        case NonRelativistic:
            return "NonRelativistic";
        case Relativistic:
            return "Relativistic";
        }
        return "Unknown";
    }

    // Algorithms Enums
    enum Algorithm { Base, MFP, Veto, Continuous, CurvedVeto, Hazard };
    std::string ToString(const Algorithm &type) {
        switch(type) {
        case Base:
            return "Base";
        case MFP:
            return "MFP";
        case Veto:
            return "Veto";
        case Continuous:
            return "Continuous";
        case CurvedVeto:
            return "CurvedVeto";
        case Hazard:
            return "Hazard";
        }
        return "Unknown";
    }

    /// Options of the curved-path algorithms (CurvedVeto = option a, Hazard = option b)
    struct CurvedOptions {
        /// Thinning: a Pauli-blocked candidate leaves the channel running (b) or retries the
        /// atom Poisson(X) times (a). SingleAttempt: one Pauli test per pair (b) or per
        /// closest approach (a), as in Algorithm: Veto. Default: Thinning for Hazard,
        /// SingleAttempt for CurvedVeto.
        enum class Pauli { Default, Thinning, SingleAttempt };
        /// Momentum at which the pair cross section is evaluated: the asymptotic momentum of
        /// the propagating particle (constant between collisions) or the local one
        enum class Sigma { Asymptotic, Local };
        /// Majorant of the Hazard algorithm: adaptive per step, one constant per pair over the
        /// whole path, or the per-step one inflated by 2 (for invariance tests)
        enum class Majorant { PerStep, Constant, Inflated };
        Pauli pauli{Pauli::Default};
        Sigma sigma{Sigma::Asymptotic};
        Majorant majorant{Majorant::PerStep};
        double path_step{0.04};         ///< integrator step [fm/c]; Tao (omega = 20) is unstable above ~0.05
        double max_step_exponent{1.0};  ///< subdivide while the majorant exponent exceeds this
        bool tabulate{true};            ///< propagate in a bicubic table of the potential
    };

    /// Counters filled by the event-driven algorithms (reset by every Evolve)
    struct Diagnostics {
        size_t candidates{}, accepted{}, pauli_blocked{}, collisions{}, violations{};
        size_t subdivisions{}, minima{}, prescreened{}, contacts{}, captured{}, unresolved{};
        double max_ratio{}, weight{1}, energy_mismatch{};
        /// First accepted (Pauli-allowed) collision of the initially propagating particle
        bool first_set{false};
        double first_time{}, first_separation{}, first_b{}, first_z{};
        size_t first_partner{};
        ThreeVector first_position{};
    };
    /// @name Constructor and Destructor
    ///@{

    /// Create the Cascade object
    ///@param interactions: The interaction model for pp, pn, and np interactions
    ///@param prob: The interaction probability function to be used
    ///@param dist: The maximum distance step to take when propagating
    /// TODO: Should the ProbabilityType be part of the interaction class or the cascade class?
    Cascade() = default;
    Cascade(InteractionHandler, const ProbabilityType &, Algorithm, const InMedium &,
            std::string_view decay_file = default_decay_file, bool potential_prob = false,
            double dist = 0.03);
    Cascade(Cascade &&) = default;
    Cascade &operator=(Cascade &&) = default;

    /// Default destructor
    ~Cascade() = default;
    ///@}

    /// @name Getters
    ///@{

    // TODO: Convert to InteractionHandler
    /// Get the name of the interaction model used
    ///@return std::string: Name of the interaction model
    // std::string InteractionModel() const { return m_interactions->Name(); }

    /// Get the probability model used
    ///@return std::string: Name of the probability model
    std::string ProbabilityModel() const { return m_probability_name; }

    /// Get InMedium setting used
    ///@return std::string: InMedium setting
    std::string InMediumSetting() const {
        std::string name;
        switch(m_medium) {
        case None:
            name = "None";
            break;
        case NonRelativistic:
            name = "NonRelativistic";
            break;
        case Relativistic:
            name = "Relativistic";
            break;
        }
        return name;
    }

    /// Get potential prop option
    ///@return bool: PotentialProp option
    bool UsePotentialProp() const { return m_potential_prop; }

    /// Get step size
    ///@return double: default step size
    double StepSize() const { return distance; }

    /// @name Functions
    ///@{

    /// Give a random particle in the nucleus a kick. The kick is defined by a input
    /// four vector. To determine whether a proton or a neutron is kicked, the cross-section
    /// for protons and neutrons needs to be supplied.
    ///@param particles: The list of particles inside the nucleus
    ///@param energyTransfer: The energy transfered to the nucleon during the kick
    ///@param sigma: An array representing the cross-section for different kicked nucleons
    void Kick(Event &, const FourVector &, const std::array<double, 2> &);

    /// Reset the cascade internal variables for the next cascade
    void Reset();

    /// Helper function to make a specific nucleon as the kicked nucleon
    ///@param idx: The index of the particle that has been kicked
    void SetKicked(const std::size_t &idx) { kickedIdxs.insert(idx); }

    /// Simulate the cascade on an event until all particles either escape,
    /// are recaptured, or are in the background.
    ///@param event: The event to run the cascade evolution on
    ///@param nucleus: The nucleus to use during evolution
    ///@param maxSteps: The maximum steps to take in the cascade
    void Evolve(Event &, Nucleus *, const std::size_t &maxSteps = cMaxSteps);

    /// Simulate evolution of a kicked particle until it interacts for the
    /// first time with another particle, accumulating the total distance
    /// traveled by the kicked particle before it interacts.
    ///@param nucleus: The nucleus to evolve according to the mean free path calculation
    ///@param maxSteps: The maximum steps to take in the particle evolution
    void MeanFreePath(Event &, Nucleus *, const std::size_t &maxSteps = cMaxSteps);

    /// Simulate the cascade until all particles either escape, are recaptured, or are in
    /// the background. This is done according to the NuWro algorithm.
    ///@param nucleus: The nucleus to evolve according to the NuWro method of cascade
    ///@param maxSteps: The maximum steps to take in the particle evolution
    void NuWro(Event &, Nucleus *, const std::size_t &maxSteps = cMaxSteps);

    /// Simulate evolution of a kicked particle until it interacts for the
    /// first time with another particle, accumulating the total distance
    /// traveled by the kicked particle before it interacts.
    ///@param nucleus: The nucleus to evolve according to the mean free path calculation
    ///@param maxSteps: The maximum steps to take in the particle evolution
    void MeanFreePath_NuWro(Event &, Nucleus *, const std::size_t &maxSteps = cMaxSteps);
    /// Symplectic integrator used by PotentialProp for one particle: Tao's extended phase
    /// space, omega = 20, with the cascade's Hamiltonian derivatives.
    static SymplecticIntegrator MakeIntegrator(const Particle &, std::shared_ptr<Potential>);

    void SetCurvedOptions(const CurvedOptions &opts) { m_curved = opts; }
    /// Switch Pauli blocking of cascade final states off (validation runs); default on
    void SetPauliBlocking(bool on) { m_pauli = on; }
    /// Largest integrator step [fm/c] for PotentialProp in the step-based algorithms: a cascade
    /// step is subdivided into equal integrator steps no longer than this (0: no subdivision,
    /// the pre-fix behaviour). Tao's integrator with omega = 20 is unstable above ~0.05 fm/c.
    void SetMaxIntegratorStep(double step) { m_max_integrator_step = step; }
    /// Capture nucleons whose orbit is trapped behind the centrifugal barrier (PotentialProp,
    /// step-based algorithms); default on
    void SetCaptureTrapped(bool on) { m_capture_trapped = on; }
    const CurvedOptions &GetCurvedOptions() const { return m_curved; }
    const Diagnostics &GetDiagnostics() const { return m_diag; }

    /// No-collision probability S(t) of propagating particle idx on its path through the
    /// frozen background of the event, for the configured algorithm (Veto/CurvedVeto: product
    /// over closest approaches; Continuous: box hazards; Hazard: exp of the integrated pair
    /// hazards), with the Pauli-allowed fraction estimated from npauli sampled final states at
    /// each pair (npauli = 0: no Pauli blocking). Ported from the vetocascade branch.
    std::vector<double> SurvivalProfile(Event &, Nucleus *, std::size_t,
                                        const std::vector<double> &times, std::size_t npauli = 0);
    ///@}
  private:
    // Functions
    std::size_t GetInter(Particles &, const Particle &, double &stepDistance);
    void AdaptiveStep(const Particles &, const double &) noexcept;
    bool BetweenPlanes(const ThreeVector &, const ThreeVector &, const ThreeVector &,
                       double) const noexcept;
    const ThreeVector Project(const ThreeVector &, const ThreeVector &,
                              const ThreeVector &) const noexcept;
    const InteractionDistances AllowedInteractions(Particles &, const std::size_t &) noexcept;
    double GetXSec(Event &, size_t, size_t) const;
    std::size_t Interacted(Event &, size_t, const InteractionDistances &) noexcept;
    void Escaped(Particles &);
    void FinalizeMomentum(Event &, Particles &, size_t, size_t) noexcept;
    bool PauliBlocking(const Particle &) const noexcept;
    bool Absorption(Event &, Particle &, Particle &) noexcept;
    void AddIntegrator(size_t, const Particle &);
    void Propagate(size_t, Particle *);
    void PropagateSpace(size_t, Particle *, double);
    std::set<size_t> InitializeIntegrator(Event &);
    bool TrappedInBarrier(const Particle &) const;
    void UpdateKicked(Particles &, std::set<size_t> &);
    void Validate(Event &);
    size_t BaseAlgorithm(size_t, Event &);
    size_t MFPAlgorithm(size_t, Event &);
    double InMediumCorrection(const Particle &, const Particle &) const;
    void PropagateAll(Particles &, double) const;
    bool HasInteraction(Event &, size_t, size_t) const;
    bool Decay(Event &, size_t) const;
    bool DecayNow(Event &, size_t) const;
    double Lifetime(const Particle &) const;

    // Event-driven algorithms (Veto: one trial per pair at closest approach,
    // Continuous: a smooth hazard per pair along the line)
    void EvolveVeto(Event &, Nucleus *);
    void ScheduleTrials(Event &, size_t, bool formed = false);
    std::optional<double> SmoothPairTime(double, double, double, double, double) const;
    void ScheduleDecay(const Particle &, size_t);
    double ExitTime(const Particle &) const;

    // Curved-path algorithms (CurvedVeto: one trial per local minimum of the pair distance on
    // the dense path; Hazard: continuous pair hazard kappa(r) with a textbook veto)
    struct Track {
        CurvedPath path;
        double t_formed{};
        bool potential{false};
    };
    void EvolveCurved(Event &, Nucleus *);
    bool StartTrack(Event &, size_t, bool with_potential);
    void ScheduleCurvedPairs(Event &, size_t);
    void SchedulePair(Event &, size_t, size_t, double);
    void SyncParticles(Particles &, double);
    void FinishTrack(Particles &, size_t);
    double PairSigma(Event &, size_t, size_t);
    std::optional<double> NextHazardTime(size_t, const ThreeVector &, double sigma, double start);
    std::optional<double> NextContactTime(size_t, const ThreeVector &, double sigma, double start);
    std::optional<double> NextMinimumTrial(size_t, size_t, const ThreeVector &, double sigma,
                                           double start);
    std::optional<Track> MakeTrack(Event &, size_t, bool with_potential, double start) const;
    double HazardMajorant(Track &, const ThreeVector &, double sigma, double ta, double *tb);
    double PairHazard(Track &, const ThreeVector &, double sigma, double t);
    double LocalSigma(double t);
    double KappaBound(double r, double sigma_lo, double sigma_hi) const;
    struct Bracket {
        double ta, tb, dlow;
    };
    std::vector<Bracket> CurvedMinima(Track &, const ThreeVector &, double start, double stop);
    std::vector<double> CumulativeHazard(Track &, const ThreeVector &, double sigma, double start,
                                         const std::vector<double> &times);
    CurvedOptions::Pauli PauliMode() const;
    std::shared_ptr<Potential> CurvedPotential() const;
    void RecordFirst(const Particles &, size_t, size_t);

    // Variables
    std::set<std::size_t> kickedIdxs;
    double distance{}, timeStep{}, currentTime{};
    InteractionHandler m_interactions{};
    DecayHandler m_decays;
    std::function<double(double, double)> probability;
    std::function<size_t(Cascade *, size_t, Event &)> algorithm;
    Algorithm m_algorithm{Algorithm::Base};
    Nucleus *m_nucleus;
    InMedium m_medium;
    bool m_potential_prop;
    std::map<size_t, SymplecticIntegrator> integrators;
    std::string m_probability_name;
    PairKernel m_kernel{};
    CurvedOptions m_curved{};
    Diagnostics m_diag{};
    std::map<size_t, Track> m_tracks;
    std::map<std::pair<size_t, size_t>, double> m_pair_sigma;
    std::map<std::pair<size_t, size_t>, size_t> m_atom_attempts;
    std::set<size_t> m_first_kicked;
    bool m_pauli{true};
    Event *m_event{nullptr};
    size_t m_chan_idx{}, m_chan_hit{};
    bool m_local_sigma{false};
    double m_max_integrator_step{0.04};
    bool m_capture_trapped{true};
    mutable std::shared_ptr<Potential> m_table{};
    mutable const Potential *m_table_source{nullptr};
    ThreeVector m_first_q0{}, m_first_v0{};

    /// Marks a queue entry as a decay trial rather than a pair trial
    static constexpr size_t cDecayTrial = SIZE_MAX;
    /// Marks the end of a formation zone in the Continuous algorithm
    static constexpr size_t cFormedTrial = SIZE_MAX - 1;
    /// Marks the exit of a particle from the nucleus in the curved-path algorithms
    static constexpr size_t cExitTrial = SIZE_MAX - 2;

    struct queue_entry {
        double time;
        std::pair<size_t, size_t> idxs;
        bool operator>(const queue_entry &rhs) const { return time > rhs.time; }
    };

    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<queue_entry>>
        m_time_steps;
};

} // namespace achilles

namespace YAML {
template <> struct convert<achilles::Cascade> {
    static bool decode(const Node &node, achilles::Cascade &cascade) {
        auto handler = node["Interactions"].as<achilles::InteractionHandler>();
        auto probType = node["Probability"].as<achilles::Cascade::ProbabilityType>();
        auto mediumType = node["InMedium"].as<achilles::Cascade::InMedium>();
        auto potentialProp = node["PotentialProp"].as<bool>();
        auto distance = node["Step"].as<double>();
        auto algorithm = node["Algorithm"].as<achilles::Cascade::Algorithm>();
        std::string_view decay_file = achilles::default_decay_file;
        if(node["DecayFile"]) decay_file = node["DecayFile"].as<std::string_view>();
        cascade = achilles::Cascade(std::move(handler), probType, algorithm, mediumType, decay_file,
                                    potentialProp, distance);
        if(node["PauliBlocking"]) cascade.SetPauliBlocking(node["PauliBlocking"].as<bool>());
        if(node["CaptureTrapped"]) cascade.SetCaptureTrapped(node["CaptureTrapped"].as<bool>());
        if(node["MaxIntegratorStep"])
            cascade.SetMaxIntegratorStep(node["MaxIntegratorStep"].as<double>());
        if(node["Curved"]) {
            const auto curved = node["Curved"];
            achilles::Cascade::CurvedOptions opts;
            using Opt = achilles::Cascade::CurvedOptions;
            if(curved["Pauli"]) {
                const auto mode = curved["Pauli"].as<std::string>();
                if(mode == "Thinning")
                    opts.pauli = Opt::Pauli::Thinning;
                else if(mode == "SingleAttempt")
                    opts.pauli = Opt::Pauli::SingleAttempt;
                else
                    return false;
            }
            if(curved["Sigma"]) {
                const auto mode = curved["Sigma"].as<std::string>();
                if(mode == "Asymptotic")
                    opts.sigma = Opt::Sigma::Asymptotic;
                else if(mode == "Local")
                    opts.sigma = Opt::Sigma::Local;
                else
                    return false;
            }
            if(curved["Majorant"]) {
                const auto mode = curved["Majorant"].as<std::string>();
                if(mode == "PerStep")
                    opts.majorant = Opt::Majorant::PerStep;
                else if(mode == "Constant")
                    opts.majorant = Opt::Majorant::Constant;
                else if(mode == "Inflated")
                    opts.majorant = Opt::Majorant::Inflated;
                else
                    return false;
            }
            if(curved["PathStep"]) opts.path_step = curved["PathStep"].as<double>();
            if(curved["Tabulate"]) opts.tabulate = curved["Tabulate"].as<bool>();
            if(curved["MaxStepExponent"])
                opts.max_step_exponent = curved["MaxStepExponent"].as<double>();
            cascade.SetCurvedOptions(opts);
        }
        return true;
    }
};

template <> struct convert<achilles::Cascade::ProbabilityType> {
    static bool decode(const Node &node, achilles::Cascade::ProbabilityType &type) {
        if(node.as<std::string>() == "Gaussian")
            type = achilles::Cascade::ProbabilityType::Gaussian;
        else if(node.as<std::string>() == "Pion")
            type = achilles::Cascade::ProbabilityType::Pion;
        else if(node.as<std::string>() == "Cylinder")
            type = achilles::Cascade::ProbabilityType::Cylinder;
        else
            return false;
        return true;
    }
};

template <> struct convert<achilles::Cascade::InMedium> {
    static bool decode(const Node &node, achilles::Cascade::InMedium &type) {
        if(node.as<std::string>() == "None")
            type = achilles::Cascade::InMedium::None;
        else if(node.as<std::string>() == "NonRelativistic")
            type = achilles::Cascade::InMedium::NonRelativistic;
        else if(node.as<std::string>() == "Relativistic")
            type = achilles::Cascade::InMedium::Relativistic;
        else
            return false;
        return true;
    }
};

template <> struct convert<achilles::Cascade::Algorithm> {
    static bool decode(const Node &node, achilles::Cascade::Algorithm &type) {
        if(node.as<std::string>() == "Base")
            type = achilles::Cascade::Algorithm::Base;
        else if(node.as<std::string>() == "MFP")
            type = achilles::Cascade::Algorithm::MFP;
        else if(node.as<std::string>() == "Veto")
            type = achilles::Cascade::Algorithm::Veto;
        else if(node.as<std::string>() == "Continuous")
            type = achilles::Cascade::Algorithm::Continuous;
        else if(node.as<std::string>() == "CurvedVeto")
            type = achilles::Cascade::Algorithm::CurvedVeto;
        else if(node.as<std::string>() == "Hazard")
            type = achilles::Cascade::Algorithm::Hazard;
        else
            return false;
        return true;
    }
};
} // namespace YAML

#endif // end of include guard: CASCADE_HH
