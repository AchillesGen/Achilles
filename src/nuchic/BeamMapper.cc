#include "nuchic/BeamMapper.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Beams.hh"

using nuchic::BeamMapper;

void BeamMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    auto beam_id = *m_beam -> BeamIDs().begin();
    point[m_idx] = m_beam -> Flux(beam_id, rans);
}

double BeamMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    auto beam_id = *m_beam -> BeamIDs().begin();
    return 1.0/m_beam -> GenerateWeight(beam_id, point[m_idx], rans);
}


size_t BeamMapper::NDims() const {
    return static_cast<size_t>(m_beam -> NVariables());
}
