#pragma once

#include "Achilles/MultiChannel.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Statistics.hh"
#include "Achilles/Unweighter.hh"

namespace achilles {

class Beam;
class FourVector;
class NuclearModel;
class Nucleus;

class Process {
    public:
        Process(Process_Info info) : m_info{std::move(info)} {}
        double TotalCrossSection() const { return m_xsec.Mean(); }
        Process_Info& Info() { return m_info; }
        const Process_Info& Info() const { return m_info; }
        double CrossSection(const Event&, NuclearModel*);
        void SelectInitialState(Event&) const;

    private:
        Process_Info m_info;
        StatsData m_xsec{};
};

class ProcessGroup {
    public:
        ProcessGroup(std::shared_ptr<NuclearModel> nuc_model) : m_nuc_model{std::move(nuc_model)} {}
        double CrossSection(Event&);
        Process SelectProcess(double ran);

        // Handling individual processes
        const Process& GetProcess(size_t i) const { return m_processes[i]; }
        Process& GetProcess(size_t i) { return m_processes[i]; }
        void AddProcess(Process process) { m_processes.push_back(std::move(process)); }

        // Handling physics objects
        NuclearModel* GetNuclearModel() { return m_nuc_model.get(); }
        Beam* GetBeam() { return m_beam.get(); }
        Nucleus* GetNucleus() { return m_nucleus.get(); }

        // Initialize processes and process groups
        static std::vector<ProcessGroup> ConstructProcessGroups(const YAML::Node&);

        Integrand<FourVector>& GetIntegrand() { return m_integrand; }
        const Integrand<FourVector>& GetIntegrand() const { return m_integrand; }
        void Optimize();
        void Summary();
        Event GenerateEvent();

    private:
        // Physics components
        std::vector<Process> m_processes;
        std::shared_ptr<NuclearModel> m_nuc_model;
        std::shared_ptr<Beam> m_beam;
        std::shared_ptr<Nucleus> m_nucleus;

        // Numerical components
        MultiChannel m_integrator;
        Integrand<FourVector> m_integrand;
        std::vector<std::unique_ptr<Unweighter>> m_unweighters;

        // Parameters
        double m_maxweight;
        bool b_optimize;
};

}
