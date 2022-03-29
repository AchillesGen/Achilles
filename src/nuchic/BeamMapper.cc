#include "nuchic/BeamMapper.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Beams.hh"

using nuchic::BeamMapper;

void BeamMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) {
    auto beam_id = *m_beam -> BeamIDs().begin();
    point[m_idx] = m_beam -> Flux(beam_id, rans);
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
}

double BeamMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) {
    auto beam_id = *m_beam -> BeamIDs().begin();
    auto wgt = m_beam -> GenerateWeight(beam_id, point[m_idx], rans);
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  Beam weight = {}", wgt);
    return 1.0/wgt;
}


size_t BeamMapper::NDims() const {
    return static_cast<size_t>(m_beam -> NVariables());
}
