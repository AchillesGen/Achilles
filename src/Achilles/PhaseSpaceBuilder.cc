#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/PhaseSpaceFactory.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/FinalStateMapper.hh"

#ifdef ENABLE_BSM
#include "plugins/Sherpa/Channels.hh"
#endif

using achilles::PSBuilder;

PSBuilder& PSBuilder::Beam(std::shared_ptr<achilles::Beam> beam, size_t idx) {
    phase_space->lbeam = std::make_shared<BeamMapper>(idx, beam); 
    return *this;
}

PSBuilder& PSBuilder::Hadron(const std::string &mode, const std::vector<double> &masses, size_t idx) {
    phase_space->hbeam = PSFactory<HadronicBeamMapper, size_t>::Build(mode, idx);
    phase_space->hbeam->SetMasses(masses);
    return *this;
}

PSBuilder& PSBuilder::FinalState(const std::string &channel, const std::vector<double> &masses2) {
    phase_space->main = PSFactory<FinalStateMapper, std::vector<double>>::Build(channel, masses2);
    return *this;
}

#ifdef ENABLE_BSM
PSBuilder& PSBuilder::SherpaFinalState(const std::string &channel, const std::vector<double> &masses2) {
    auto sherpaMap = PSFactory<PHASIC::Channels, std::vector<double>>::Build(channel, masses2);
    phase_space->main = std::make_unique<SherpaMapper>(m_nlep+m_nhad-2, std::move(sherpaMap));
    return *this;
}

PSBuilder& PSBuilder::GenFinalState(std::unique_ptr<PHASIC::Channels> channel) {
    phase_space->main = std::make_unique<SherpaMapper>(m_nlep+m_nhad-2, std::move(channel));
    return *this;
}
#endif
