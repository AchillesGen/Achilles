#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/PhaseSpaceFactory.hh"

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
    phase_space->hbeam = PSFactory<HadronicBeamMapper, size_t>::Build(mode, idx);
    // TODO: Right now we only need the intial hadronic mass for Coherent, but we may need to change
    // this
    // TODO: Also, assumes that the zeroth index is the beam particle
    auto masses = {0.0, ParticleInfo(m_info.m_hadronic.first[0]).Mass()};
    phase_space->hbeam->SetMasses(masses);
    return *this;
}

PSBuilder &PSBuilder::FinalState(const std::string &channel,
                                 std::optional<double> gauge_boson_mass) {
    phase_space->main =
        PSFactory<FinalStateMapper, std::vector<double>>::Build(channel, m_info.Masses());
    if(gauge_boson_mass.has_value()) phase_space->main->SetGaugeBosonMass(gauge_boson_mass.value());
    return *this;
}

#ifdef ACHILLES_SHERPA_INTERFACE
PSBuilder &PSBuilder::SherpaFinalState(const std::string &channel) {
    auto sherpaMap =
        PSFactory<PHASIC::Channels, std::vector<double>>::Build(channel, m_info.Masses());
    phase_space->main = std::make_unique<SherpaMapper>(m_nlep + m_nhad - 2, std::move(sherpaMap));
    return *this;
}

PSBuilder &PSBuilder::GenFinalState(std::unique_ptr<PHASIC::Channels> channel) {
    phase_space->main = std::make_unique<SherpaMapper>(m_nlep + m_nhad - 2, std::move(channel));
    return *this;
}
#endif
