#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Event.hh"

using achilles::BeamMapper;

void BeamMapper::GeneratePoint(Event& event, const std::vector<double> &rans) {
    auto beam_id = *m_beam->BeamIDs().begin();
    // TODO: Should Masses().back() be the mass of the final state hadronic system or the initial?
    event.Momentum()[m_idx] = m_beam->Flux(beam_id, rans, (Smin() - Masses()[1]) / (2 * sqrt(Masses()[1])));
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
#endif
}

double BeamMapper::GenerateWeight(const Event& event, std::vector<double> &rans) {
    auto beam_id = *m_beam->BeamIDs().begin();
    // TODO: Should Masses().back() be the mass of the final state hadronic system or the initial?
    double wgt = m_beam->GenerateWeight(beam_id, event.Momentum()[m_idx], rans,
                                      (Smin() - Masses()[1]) / (2 * sqrt(Masses()[1])));
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
    spdlog::trace("  Beam weight = {}", wgt);
#endif
    return 1.0 / wgt;
}

size_t BeamMapper::NDims() const {
    return static_cast<size_t>(m_beam->NVariables());
}