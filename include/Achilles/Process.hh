#pragma once

#include "Achilles/CombinedCuts.hh"
#include "Achilles/MultiChannel.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Statistics.hh"
#include "Achilles/Unweighter.hh"
#include "Achilles/XSecBackend.hh"

#include <optional>

namespace achilles {

class Beam;
class FourVector;
class NuclearModel;
class Nucleus;
class SherpaInterface;

class Process {
  public:
    Process(ProcessInfo info, std::unique_ptr<Unweighter> unweighter)
        : m_info{std::move(info)}, m_unweighter{std::move(unweighter)} {}
    double TotalCrossSection() const { return m_xsec.Mean(); }
    ProcessInfo &Info() { return m_info; }
    const ProcessInfo &Info() const { return m_info; }
    size_t NInitialStates(Nucleus *) const;
    void SetupHadrons(Event &) const;
    void AddWeight(double weight) {
        m_unweighter->AddEvent(weight);
        m_xsec += weight;
    }
    double Unweight(double weight) { return m_unweighter->AcceptEvent(weight); }
    double MaxWeight() { return m_unweighter->MaxValue(); }
    void ExtractMomentum(const Event &, FourVector &, std::vector<FourVector> &,
                         std::vector<FourVector> &, std::vector<FourVector> &) const;
    double UnweightEff() const { return m_xsec.Mean() / m_unweighter->MaxValue(); }

  private:
    ProcessInfo m_info;
    StatsData m_xsec{};
    std::unique_ptr<Unweighter> m_unweighter;
};

class ProcessGroup {
  public:
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
    void SetupLeptons(Event &) const;

    // Initialize processes and process groups
    static std::map<size_t, ProcessGroup> ConstructProcessGroups(const YAML::Node &, NuclearModel *,
                                                                 std::shared_ptr<Beam>,
                                                                 std::shared_ptr<Nucleus>);

    Integrand<FourVector> &GetIntegrand() { return m_integrand; }
    const Integrand<FourVector> &GetIntegrand() const { return m_integrand; }
    void SetupIntegration(const YAML::Node &);
    void Optimize();
    void Summary() const;
    Event GenerateEvent();
    Event SingleEvent(const std::vector<FourVector> &, double);
    double MaxWeight() const { return m_maxweight; }

  private:
    // Physics components
    std::vector<Process> m_processes;
    std::shared_ptr<Beam> m_beam;
    std::shared_ptr<Nucleus> m_nucleus;
    std::unique_ptr<XSecBackend> m_backend;
    CutCollection m_cuts;

    // Numerical components
    MultiChannel m_integrator;
    Integrand<FourVector> m_integrand;

    // Parameters
    std::vector<double> m_process_weights;
    double m_maxweight{};
    bool b_optimize{};
};

} // namespace achilles
