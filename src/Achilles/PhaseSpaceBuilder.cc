#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/HadronicMapper.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/Channels.hh"
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

PSBuilder &PSBuilder::Hadron(const std::string &mode, size_t idx) {
    // BUG: The idx parameter needs to be cast to size_t otherwise it is an error about rvalue and
    // lvalues
    phase_space->hbeam = Factory<HadronicBeamMapper, const ProcessInfo &, size_t>::Initialize(
        mode, m_info, static_cast<size_t>(idx));
    return *this;
}

PSBuilder &PSBuilder::FinalState(const std::string &channel,
                                 std::optional<double> gauge_boson_mass) {
    phase_space->main =
        Factory<FinalStateMapper, std::vector<double>>::Initialize(channel, m_info.Masses());
    if(gauge_boson_mass.has_value()) phase_space->main->SetGaugeBosonMass(gauge_boson_mass.value());
    return *this;
}

#ifdef ACHILLES_SHERPA_INTERFACE
PSBuilder &PSBuilder::SherpaFinalState(const std::string &channel) {
    auto sherpaMap =
        Factory<PHASIC::Channels, std::vector<double>>::Initialize(channel, m_info.Masses());
    phase_space->main = std::make_unique<SherpaMapper>(m_nlep + m_nhad - 2, std::move(sherpaMap));
    return *this;
}

PSBuilder &PSBuilder::GenFinalState(std::unique_ptr<PHASIC::Channels> channel) {
    phase_space->main = std::make_unique<SherpaMapper>(m_nlep + m_nhad - 2, std::move(channel));
    return *this;
}
#endif
