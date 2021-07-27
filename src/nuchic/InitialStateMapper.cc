#include "nuchic/InitialStateMapper.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Beams.hh"

using nuchic::InitialMapper;

void InitialMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    // Separate random numbers
    const std::vector<double> beamRans(rans.begin(), rans.begin() + m_beam -> NVariables());  
    const std::vector<double> qRans(rans.begin() + m_beam -> NVariables(), rans.end());

    // Initialize lepton from beam
    point[1] = m_beam -> Flux(PID::electron(), beamRans);

    // Initialize dummy lepton to calculate Q from
    const double dE = point[1].E();
    constexpr double dcosT = 2;
    constexpr double dphi = 2*M_PI;
    const double E = dE*qRans[0];
    const double cosT = dcosT*qRans[1] - 1;
    const double sinT = sqrt(1-cosT*cosT);
    const double phi = dphi*qRans[2];
    const FourVector dummyLepton{E, E*sinT*cos(phi), E*sinT*sin(phi), E*cosT};
    point[0] = point[1] - dummyLepton;
}

double InitialMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    // Initialize random numbers
    const std::vector<double> beamRans(rans.begin(), rans.begin() + m_beam -> NVariables());  
    const std::vector<double> qRans(rans.begin() + m_beam -> NVariables(), rans.end());

    // Get weights
    constexpr double dcosT = 2;
    constexpr double dphi = 2*M_PI;
    return dcosT*dphi*point[0].E();
}
