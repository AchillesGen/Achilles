#include "Achilles/XSecBackend.hh"
#include "Achilles/FourVector.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/SherpaInterface.hh"
#endif

#include <stdexcept>

using achilles::XSecBackend;
using achilles::XSecBackendFactory;
using achilles::XSecBuilder;

double achilles::DefaultBackend::CrossSection(const Event &event, const Process &process) const {
    const auto& process_info = process.Info();
    static constexpr size_t lepton_in_end_idx = 1;
    const size_t hadron_end_idx = process_info.m_hadronic.first.size()+lepton_in_end_idx;
    const size_t lepton_end_idx = process_info.m_leptonic.second.size()+hadron_end_idx;
    auto momentum = event.Momentum();
    FourVector lepton_momentum_in = momentum[0];
    std::vector<FourVector> hadron_momentum_in(momentum.begin()+lepton_in_end_idx,
                                               momentum.begin()+static_cast<int>(nhadrons_in));
    std::vector<FourVector> lepton_momentum_out(momentum.begin()+nleptons_in,
                                               momentum.begin()+nleptons_in+nhadrons_in);
    
    return 0;
}

void achilles::DefaultBackend::AddProcess(const Process &process) {
    const auto& process_info = process.Info();
    if(process_info.m_leptonic.second.size() != 1)
        throw std::runtime_error("Achilles::DefaultBackend: Can only handle single lepton final states");

    std::pair<PID,PID> leptons{process_info.m_leptonic.first, process_info.m_leptonic.second[0]};
    if(m_currents.find(leptons) != m_currents.end()) {
        spdlog::debug("Achilles::DefaultBackend: Leptonic current for {} -> {} already created",
                      leptons.first, leptons.second);
        return;
    }
    
    LeptonicCurrent current;
    current.Initialize(process_info);

    // TODO: Clean up how the form factors are loaded. Make part of the backend?
    if(form_factors.size() == 0)
        form_factors = current.GetFormFactor(); 
}

#ifdef ACHILLES_SHERPA_INTERFACE
double achilles::BSMBackend::CrossSection(const Event &) const {
    return 0;
}

double achilles::SherpaBackend::CrossSection(const Event &) const {
    return 0;
}
#endif

XSecBuilder::XSecBuilder(const std::string &name) : m_backend{XSecBackendFactory::Initialize(name)} {}

XSecBuilder& XSecBuilder::AddOptions(const YAML::Node &node) {
    m_backend -> SetOptions(node);
    return *this;
}

XSecBuilder& XSecBuilder::AddSherpa(achilles::SherpaInterface *sherpa) {
    m_backend -> SetSherpa(sherpa);
    return *this;
}

XSecBuilder& XSecBuilder::AddNuclearModel(std::shared_ptr<NuclearModel> model) {
    m_backend -> SetNuclearModel(model);
    return *this;
}

std::unique_ptr<XSecBackend> XSecBuilder::build() {
    if(!m_backend -> Validate())
        throw std::runtime_error("Backend is not properly configured!");
    return std::move(m_backend); 
}
