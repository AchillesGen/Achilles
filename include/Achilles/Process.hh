#pragma once

#include "Achilles/CombinedCuts.hh"
#include "Achilles/MultiChannel.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Statistics.hh"
#include "Achilles/Unweighter.hh"
#include "Achilles/XSecBackend.hh"
#include "Achilles/XSecSpline.hh"
#include "fmt/base.h"

#include <optional>

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "Could not find includes: <filesystem> or <experimental/filesystem>"
#endif

namespace achilles {

class Beam;
class FourVector;
class NuclearModel;
class Nucleus;
class SherpaInterface;
class Settings;

struct ProcessMetadata {
    int id;
    std::string name, description, inspireHEP;
};

class Process {
  public:
    Process(ProcessInfo info, std::unique_ptr<Unweighter> unweighter)
        : m_info{std::move(info)}, m_unweighter{std::move(unweighter)} {}
    Process(Process &&) = default;
    MOCK ~Process() = default;
    double TotalCrossSection() const { return m_xsec.Mean(); }
    ProcessInfo &Info() { return m_info; }
    MOCK const ProcessInfo &Info() const { return m_info; }
    void SetupHadrons(Event &) const;
    MOCK void AddWeight(double weight) {
        m_unweighter->AddEvent(weight);
        m_xsec += weight;
    }
    MOCK double Unweight(double weight) { return m_unweighter->AcceptEvent(weight); }
    double MaxWeight() { return m_unweighter->MaxValue(); }
    MOCK void ExtractMomentum(const Event &, FourVector &, std::vector<FourVector> &,
                              std::vector<FourVector> &, std::vector<FourVector> &,
                              std::vector<FourVector> &) const;
    MOCK void ExtractParticles(const Event &, Particle &, std::vector<Particle> &,
                               std::vector<Particle> &, std::vector<Particle> &,
                               std::vector<Particle> &) const;
    FourVector ExtractQ(const Event &) const;
    double UnweightEff() const { return std::abs(m_xsec.Mean()) / m_unweighter->MaxValue(); }
    bool operator==(const Process &other) const { return m_info == other.m_info; }

    // Metadata handlers
    void SetID(NuclearModel *);
    int ID() const { return m_id; }
    ProcessMetadata Metadata(XSecBackend *) const;

    // Cache results
    bool SaveState(std::ostream &) const;
    bool LoadState(std::istream &);

  private:
    // Helper functions
    refParticles SelectParticles(const refParticles &, const refParticles &,
                                 const std::vector<PID> &, const std::vector<FourVector> &,
                                 ParticleStatus) const;

    // Variables
    ProcessInfo m_info;
    StatsData m_xsec{};
    std::unique_ptr<Unweighter> m_unweighter;
    int m_id{};

    // Metadata handlers
    std::string Name(XSecBackend *) const;
    std::string Description(XSecBackend *) const;
    std::string InspireHEP(XSecBackend *) const;
};

class ProcessGroup {
  public:
    ProcessGroup() {}
    ProcessGroup(std::shared_ptr<Beam> beam, std::shared_ptr<Nucleus> nucleus)
        : m_beam{std::move(beam)}, m_nucleus{std::move(nucleus)} {}
    void CrossSection(Event &, std::optional<size_t>);
    // Select a process by weight; if required_nu is set, restrict to processes
    // whose incoming lepton matches it (used in geometry mode).
    size_t SelectProcess(std::optional<PID> required_nu = std::nullopt) const;

    // Handling individual processes
    const Process &GetProcess(size_t i) const { return m_processes[i]; }
    Process &GetProcess(size_t i) { return m_processes[i]; }
    void AddProcess(Process process) { m_processes.push_back(std::move(process)); }
    const std::vector<Process> &Processes() const { return m_processes; }
    // Absolute max-weight (same units as MaxWeight()) summed over only the
    // processes whose incoming lepton is nu; 0 if the group has no such process.
    // Used to weight group selection correctly in geometry mode.
    double NeutrinoMaxWeight(PID nu) const;

    // Handling physics objects
    Beam *GetBeam() { return m_beam.get(); }
    const Beam *GetBeam() const { return m_beam.get(); }
    Nucleus *GetNucleus() { return m_nucleus.get(); }
    const Nucleus *GetNucleus() const { return m_nucleus.get(); }
    void SetupBackend(const Settings &, std::unique_ptr<NuclearModel>, SherpaInterface *);
    void SetCuts(CutCollection cuts) { m_cuts = std::move(cuts); }
    void SetupLeptons(Event &, std::optional<size_t>) const;

    // Initialize processes and process groups
    static std::map<size_t, ProcessGroup> ConstructGroups(const Settings &, NuclearModel *,
                                                          std::shared_ptr<Beam>,
                                                          std::shared_ptr<Nucleus>);

    Integrand<FourVector> &GetIntegrand() { return m_integrand; }
    const Integrand<FourVector> &GetIntegrand() const { return m_integrand; }
    bool SetupIntegration(const Settings &);
    void Optimize();
    void Summary() const;
    Event GenerateEvent(std::optional<PID> required_nu = std::nullopt);
    Event SingleEvent(const std::vector<FourVector> &, double,
                      std::optional<PID> required_nu = std::nullopt);
    const double &MaxWeight() const { return m_maxweight; }
    double &MaxWeight() { return m_maxweight; }
    const std::vector<double> &ProcessWeights() const { return m_process_weights; }
    std::vector<double> &ProcessWeights() { return m_process_weights; }
    // Max-weight (w_max_p) of the process selected by the most recent
    // GenerateEvent/SingleEvent call. Used to apply the TotalXSecRetry weight
    // convention (max(w_max_p, raw_w)) in geometry mode.
    double LastProcessMaxWeight() const { return m_last_max_weight; }
    void SetOptimize(bool optimize) { b_optimize = optimize; }
    size_t Multiplicity() const { return m_processes[0].Info().Multiplicity(); }

    // Handling of XSecSplines for Geometry. One sigma(E) spline per incoming
    // neutrino species (the target is the group's nucleus); summed over the
    // group's processes for that species.
    void SetupSplines();
    // sigma(E) for the given incoming neutrino on this group's nucleus, summed
    // over the matching processes; 0 if the group has no such process or splines
    // are disabled (mono-energetic beam).
    double TotalXSec(double energy, PID nu) const;

    // Metadata handlers
    std::vector<ProcessMetadata> Metadata() const;
    std::vector<int> ProcessIds() const;
    // std::string UniqueID() const;

    // Cache results
    bool Save(const fs::path &) const;
    bool Load(const fs::path &);

    friend std::hash<ProcessGroup>;

  private:
    // Physics components
    std::vector<Process> m_processes;
    std::shared_ptr<Beam> m_beam = nullptr;
    std::shared_ptr<Nucleus> m_nucleus = nullptr;
    std::unique_ptr<XSecBackend> m_backend = nullptr;
    CutCollection m_cuts;

    // Numerical components
    bool NeedsOptimization() const;
    MultiChannel m_integrator;
    Integrand<FourVector> m_integrand;
    StatsData m_xsec;
    GeometryXSecSpline m_splines{};
    bool m_use_spline{false};

    // Parameters
    std::vector<double> m_process_weights;
    double m_maxweight{};
    double m_last_max_weight{};
    bool b_optimize{true}, b_calc_weights{};
};

std::vector<int> AllProcessIDs(const std::vector<ProcessGroup> &);
std::vector<ProcessMetadata> AllProcessMetadata(const std::vector<ProcessGroup> &);

// Per-group selection weights for an injected (nu, target) ray in geometry mode:
// each group's matching-species max-weight, or 0 when the group's nucleus is not
// the target. Selecting proportionally to these (rather than a group's total
// weight) keeps an injected ray's species from being biased by a group's
// unrelated processes.
std::vector<double> GeometryGroupWeights(const std::vector<ProcessGroup> &groups, PID nu,
                                         PID target);

} // namespace achilles

template <> struct std::hash<achilles::Process> {
    std::size_t operator()(const achilles::Process &p) const {
        return std::hash<achilles::ProcessInfo>{}(p.Info());
    }
};

template <> struct std::hash<achilles::ProcessGroup> {
    std::size_t operator()(const achilles::ProcessGroup &p) const;
};

namespace fmt {

template <> struct formatter<achilles::Process> {
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        return ctx.begin();
    }

    auto format(const achilles::Process &process, format_context &ctx) const
        -> format_context::iterator {
        return format_to(ctx.out(), "Process[{}]", process.Info());
    }
};

template <> struct formatter<achilles::ProcessGroup> {
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        return ctx.begin();
    }

    auto format(const achilles::ProcessGroup &group, format_context &ctx) const
        -> format_context::iterator {
        return format_to(ctx.out(), "ProcessGroup[{}]", fmt::join(group.Processes(), ", "));
    }
};

} // namespace fmt
