#pragma once

#include "Achilles/CombinedCuts.hh"
#include "Achilles/MultiChannel.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Statistics.hh"
#include "Achilles/Unweighter.hh"
#include "Achilles/XSecBackend.hh"

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
    size_t SelectProcess() const;

    // Handling individual processes
    const Process &GetProcess(size_t i) const { return m_processes[i]; }
    Process &GetProcess(size_t i) { return m_processes[i]; }
    void AddProcess(Process process) { m_processes.push_back(std::move(process)); }
    const std::vector<Process> &Processes() const { return m_processes; }

    // Handling physics objects
    Beam *GetBeam() { return m_beam.get(); }
    Nucleus *GetNucleus() { return m_nucleus.get(); }
    void SetupBackend(const YAML::Node &, std::unique_ptr<NuclearModel>, SherpaInterface *);
    void SetCuts(CutCollection cuts) { m_cuts = std::move(cuts); }
    void SetupLeptons(Event &, std::optional<size_t>) const;

    // Initialize processes and process groups
    static std::map<size_t, ProcessGroup> ConstructGroups(const YAML::Node &, NuclearModel *,
                                                          std::shared_ptr<Beam>,
                                                          std::shared_ptr<Nucleus>);

    Integrand<FourVector> &GetIntegrand() { return m_integrand; }
    const Integrand<FourVector> &GetIntegrand() const { return m_integrand; }
    bool SetupIntegration(const YAML::Node &);
    void Optimize();
    void Summary() const;
    Event GenerateEvent();
    Event SingleEvent(const std::vector<FourVector> &, double);
    const double &MaxWeight() const { return m_maxweight; }
    double &MaxWeight() { return m_maxweight; }
    void SetOptimize(bool optimize) { b_optimize = optimize; }
    size_t Multiplicity() const { return m_processes[0].Info().Multiplicity(); }

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

    // Parameters
    std::vector<double> m_process_weights;
    double m_maxweight{};
    bool b_optimize{true}, b_calc_weights{};
};

std::vector<int> AllProcessIDs(const std::vector<ProcessGroup> &);
std::vector<ProcessMetadata> AllProcessMetadata(const std::vector<ProcessGroup> &);

} // namespace achilles

template <> struct std::hash<achilles::Process> {
    std::size_t operator()(const achilles::Process &p) const {
        return std::hash<achilles::ProcessInfo>{}(p.Info());
    }
};

template <> struct std::hash<achilles::ProcessGroup> {
    std::size_t operator()(const achilles::ProcessGroup &p) const;
};
