#include "Achilles/Process.hh"
#include "Achilles/Beams.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"

using achilles::Process;
using achilles::ProcessGroup;

double Process::CrossSection(const Event &event, NuclearModel *nuclear) {
    double xsec = 0;
    return xsec;
}

void Process::SelectInitialState(Event &event) const {

}

double ProcessGroup::CrossSection(Event &event) {
    if(b_optimize) {
        return std::accumulate(m_processes.begin(), m_processes.end(), 0,
                               [&](const auto &process) {
                                    return process.CrossSection(event, GetNuclearModel(), b_optimize);
                               });
    } else {
        auto process = SelectProcess(Random::Instance().Uniform<double>(0.0, 1.0));
        return process.CrossSection(event, GetNuclearModel());
    }
}
