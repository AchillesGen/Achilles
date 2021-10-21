#include "nuchic/BeamMapper.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Beams.hh"

using nuchic::BeamMapper;

void BeamMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    point[m_idx] = m_beam -> Flux(PID::nu_muon(), rans);
}

double BeamMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    return 1.0/m_beam -> GenerateWeight(PID::nu_muon(), point[m_idx], rans);
}


size_t BeamMapper::NDims() const {
    return static_cast<size_t>(m_beam -> NVariables());
}
