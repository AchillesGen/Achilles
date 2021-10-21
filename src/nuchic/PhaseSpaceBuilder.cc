#include "nuchic/PhaseSpaceBuilder.hh"
#include "nuchic/BeamMapper.hh"
#include "nuchic/Beams.hh"
#include "nuchic/PhaseSpaceFactory.hh"
#include "nuchic/HadronicMapper.hh"
#include "nuchic/FinalStateMapper.hh"
#include "plugins/Sherpa/Channels.hh"

using nuchic::PSBuilder;

PSBuilder& PSBuilder::Beam(std::shared_ptr<nuchic::Beam> beam, size_t idx) {
    phase_space->lbeam = std::make_shared<BeamMapper>(idx, beam); 
    return *this;
}

PSBuilder& PSBuilder::Hadron(const std::string &mode, size_t idx) {
    phase_space->hbeam = PSFactory<HadronicBeamMapper, size_t>::Build(mode, idx);
    return *this;
}

PSBuilder& PSBuilder::FinalState(const std::string &channel, const std::vector<double> &masses2) {
    phase_space->main = PSFactory<FinalStateMapper, std::vector<double>>::Build(channel, masses2);
    return *this;
}

PSBuilder& PSBuilder::SherpaFinalState(const std::string &channel, const std::vector<double> &masses2) {
    auto sherpaMap = PSFactory<PHASIC::Channels, std::vector<double>>::Build(channel, masses2);
    phase_space->main = std::make_unique<SherpaMapper>(m_nlep+m_nhad-2, std::move(sherpaMap));
    return *this;
}
