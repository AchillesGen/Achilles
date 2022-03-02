#include <utility>

#include "nuchic/Constants.hh"
#include "nuchic/QuasielasticTestMapper.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/Units.hh"

using nuchic::QuasielasticTestMapper;

QuasielasticTestMapper::QuasielasticTestMapper(const YAML::Node &node, std::shared_ptr<Beam> beam)
      : m_beam{std::move(beam)} {
    mode = node["Run Mode"].as<nuchic::RunMode>();
    switch(mode) {
        case RunMode::FixedAngleEnergy:
            nvars = 4;
            m_angle = node["Angle"].as<double>()*1.0_deg;
            m_lepton_energy = node["Energy"].as<double>()*1.0_deg;
            break;
        case RunMode::FixedAngle:
            nvars = 5;
            m_angle = node["Angle"].as<double>()*1.0_deg;
            m_lepton_energy = -1;
            break;
        case RunMode::FullPhaseSpace:
            nvars = 6;
            m_angle = -1;
            m_lepton_energy = -1;
            break;
    }
    nvars += static_cast<size_t>(m_beam->NVariables());
}

void QuasielasticTestMapper::GeneratePoint(std::vector<FourVector> &mom, const std::vector<double> &rans) {
    const std::vector<double> beamRans(rans.begin(), rans.begin() + m_beam -> NVariables());
    const std::vector<double> momRans(rans.begin() + m_beam -> NVariables(), rans.end());
    mom.resize(4);

    auto beam_id = *m_beam -> BeamIDs().begin();
    mom[1] = m_beam -> Flux(beam_id, beamRans); 
   
    size_t iRan = 0;
    double phi_l = dPhi*momRans[iRan++];
    double cosT{}, sinT{}, Elepton{};
    switch(mode) {
        case RunMode::FixedAngleEnergy:
            Elepton = m_lepton_energy;
            cosT = std::cos(m_angle);
            sinT = std::sin(m_angle);
            break;
        case RunMode::FixedAngle:
            Elepton = mom[1].E()*momRans[iRan++];
            cosT = std::cos(m_angle);
            sinT = std::sin(m_angle);
            break;
        case RunMode::FullPhaseSpace:
            Elepton = mom[1].E()*momRans[iRan++];
            cosT = dCos*momRans[iRan++] - 1;
            sinT = sqrt(1-cosT*cosT);
            break;
    }
    mom[3] = {Elepton, Elepton*sinT*cos(phi_l), Elepton*sinT*sin(phi_l), Elepton*cosT};

    auto Q = mom[1] - mom[3];
    double cosT_h = dCos*momRans[iRan++] - 1;
    double sinT_h = sqrt(1-cosT_h*cosT_h);
    double phi_h = dPhi*momRans[iRan++];
    double p_h = dp*momRans[iRan];

    ThreeVector tmp = {p_h*sinT_h*cos(phi_h), p_h*sinT_h*sin(phi_h), p_h*cosT_h};
    auto tmp2 = tmp + Q.Vec3();

    double Epp = sqrt(pow(Constant::mN, 2) + tmp2.P2());
    double Ep = Constant::mN+Q.E()-Epp;
    mom[0] = FourVector(tmp, Constant::mN-Ep);
    mom[2] = FourVector(tmp2, Epp);

    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
}

double QuasielasticTestMapper::GenerateWeight(const std::vector<FourVector> &mom, std::vector<double> &rans) {
    std::vector<double> beamRans(static_cast<size_t>(m_beam -> NVariables()));
    std::vector<double> momRans(nvars - beamRans.size());

    // Calculate the weights
    double wgt = 1.0;
    auto beam_id = *m_beam -> BeamIDs().begin();
    wgt /= m_beam -> GenerateWeight(beam_id, mom[1], beamRans);  
    size_t iRan = 0;
    momRans[iRan++] = mom[3].Phi()/dPhi;
    wgt /= dPhi;
    switch(mode) {
        case RunMode::FixedAngleEnergy:
            break;
        case RunMode::FixedAngle:
            momRans[iRan++] = mom[3].E()/mom[1].E();
            wgt /= mom[1].E();
            break;
        case RunMode::FullPhaseSpace:
            momRans[iRan++] = mom[3].E()/mom[1].E();
            momRans[iRan++] = (mom[3].CosTheta()+1)/dCos;
            wgt /= (mom[1].E()*dCos);
            break;
    }

    momRans[iRan++] = (mom[0].CosTheta()+1)/dCos;
    momRans[iRan++] = mom[0].Phi()/dPhi;
    momRans[iRan++] = mom[0].P()/dp;
    wgt /= (dCos*dPhi*dp*mom[0].P2()*mom[3].E()/(mom[2].E()));
    wgt *= 16*M_PI*M_PI;

    if(mom[0].E() > Constant::mN) wgt = std::numeric_limits<double>::infinity();

    beamRans.insert(
            beamRans.end(),
            std::make_move_iterator(momRans.begin()),
            std::make_move_iterator(momRans.end()));
    swap(rans, beamRans);

    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
    spdlog::trace("  Weight: {}", wgt);

    return wgt;
}
