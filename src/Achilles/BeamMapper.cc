#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/FourVector.hh"

using achilles::BeamMapper;

void BeamMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) {
    auto beam_id = *m_beam->BeamIDs().begin();
    // TODO: 1. Resolve this arbitrary eps with a cut
    //       2. Replace the sqrt(Smin())-sqrt(Masses()[0]) with final state hadronic mass
    // static constexpr double eps=5;
    point[m_idx] =
        m_beam->Flux(beam_id, rans, (Smin() - Masses().back()) / (2 * sqrt(Masses().back())));
    // Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
}

double BeamMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) {
    auto beam_id = *m_beam->BeamIDs().begin();
    // TODO: 1. Resolve this arbitrary eps with a cut
    //       2. Replace the sqrt(Smin())-sqrt(Masses()[0]) with final state hadronic mass
    // static constexpr double eps=5;
    auto wgt = m_beam->GenerateWeight(beam_id, point[m_idx], rans,
                                      (Smin() - Masses().back()) / (2 * sqrt(Masses().back())));
    // Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    // spdlog::trace("  Beam weight = {}", wgt);
    return 1.0 / wgt;
}

size_t BeamMapper::NDims() const {
    return static_cast<size_t>(m_beam->NVariables());
}
