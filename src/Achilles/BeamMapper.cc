#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/FourVector.hh"

using achilles::BeamMapper;

void BeamMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) {
    auto beam_id = *m_beam->BeamIDs().begin();
    // TODO: Should Masses().back() be the mass of the final state hadronic system or the initial?
    point[m_idx] = m_beam->Flux(beam_id, rans, (Smin() - Masses()[1]) / (2 * sqrt(Masses()[1])));
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
#endif
}

double BeamMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) {
    auto beam_id = *m_beam->BeamIDs().begin();
    // TODO: Should Masses().back() be the mass of the final state hadronic system or the initial?
    auto wgt = m_beam->GenerateWeight(beam_id, point[m_idx], rans,
                                      (Smin() - Masses()[1]) / (2 * sqrt(Masses()[1])));
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  Beam weight = {}", wgt);
#endif
    m_density = 1.0 / wgt;
    return m_density;
}

size_t BeamMapper::NDims() const {
    return static_cast<size_t>(m_beam->NVariables());
}

std::optional<size_t> BeamMapper::FixedDim() const {
    // Geometry mode serves the injected ray, so this mapper's variable is not
    // drawn from the density it reports; the caller has to divide it back out.
    if(NDims() == 0 || !m_beam->IsFixed(*m_beam->BeamIDs().begin())) return std::nullopt;
    return 0;
}
