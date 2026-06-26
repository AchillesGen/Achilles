#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/HadronicMapper.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "Plugins/Sherpa/Channels.hh"
#endif

using achilles::PSBuilder;

PSBuilder::PSBuilder(const ProcessInfo &info) : m_info(info) {
    m_nlep = 1 + m_info.m_leptonic.second.size();
    m_nhad = m_info.m_hadronic.first.size() + m_info.m_hadronic.second.size();
    m_nspec = m_info.m_spectator.size();
    phase_space = std::make_unique<PSMapper>(m_nlep, m_nhad, m_nspec);
}

PSBuilder &PSBuilder::Beam(std::shared_ptr<achilles::Beam> beam, size_t idx) {
    phase_space->lbeam = std::make_shared<BeamMapper>(idx, beam);
    phase_space->lbeam->SetMasses(m_info.Masses());
    return *this;
}

PSBuilder &PSBuilder::Hadron(const std::string &mode) {
	try {
		// If beam has been initialized, HadronIdx can be computed
		return Hadron(mode,phase_space->HadronIdx());
	} catch(...) {
		// otherwise, HadronIdx is assumed to eventually be 1
		return Hadron(mode,1);
	}
}

PSBuilder &PSBuilder::Hadron(const std::string &mode, size_t idx) {
    // BUG: The idx parameter needs to be cast to size_t otherwise it is an error about rvalue and
    // lvalues
    phase_space->hbeam = Factory<HadronicBeamMapper, const ProcessInfo &, size_t>::Initialize(
        mode, m_info, std::move(idx));
    return *this;
}

PSBuilder &PSBuilder::FinalState(const std::string &channel,
                                 std::optional<double> gauge_boson_mass) {
	size_t idx;
	try {
		// If beam and hadrons have been initialized, FinalStateIdx can be computed
		idx=phase_space->FinalStateIdx();
	} catch(...) {
		// otherwise, FinalStateIdx is assumed to eventually be 2
		idx=2;
	}
    phase_space->finalstate =
        Factory<FinalStateMapper, std::vector<double>, double>::Initialize(channel, m_info.Masses(),idx);
    if(gauge_boson_mass.has_value()) phase_space->finalstate->SetGaugeBosonMass(gauge_boson_mass.value());
    return *this;
}

#ifdef ACHILLES_SHERPA_INTERFACE
PSBuilder &PSBuilder::SherpaFinalState(const std::string &channel) {
    auto sherpaMap =
        Factory<PHASIC::Channels, std::vector<double>>::Initialize(channel, m_info.Masses());
    phase_space->finalstate = std::make_unique<SherpaMapper>(m_nlep + m_nhad - 2, std::move(sherpaMap));
    return *this;
}

PSBuilder &PSBuilder::GenFinalState(std::unique_ptr<PHASIC::Channels> channel) {
    phase_space->finalstate = std::make_unique<SherpaMapper>(m_nlep + m_nhad - 2, std::move(channel));
    return *this;
}
#endif
