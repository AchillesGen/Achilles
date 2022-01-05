#include "plugins/Sherpa/Channels1.hh"
#include "ATOOLS/Phys/Flavour.H"
#include "plugins/Sherpa/PrintVec.hh"

using namespace PHASIC;
using namespace ATOOLS;

void C1_0::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    CE.Isotropic2Momenta(p[0]+p[1],s2,s3,p[2],p[3],ran[0],ran[1],-1,1);
    Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C1_0::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &ran) {
    double wt = CE.Isotropic2Weight(p[2],p[3],ran[0],ran[1],-1,1);
    if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,2*3.-4.);
    Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
    return 1.0/wt;
}

void C1_1::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    CE.TChannelMomenta(p[0],p[1],p[2],p[3],s2,s3,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[0],ran[1]);
    Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C1_1::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &ran) {
    double wt = CE.TChannelWeight(p[0],p[1],p[2],p[3],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[0],ran[1]);
    if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,2*3.-4.);
    Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
    return 1.0/wt;
}

void C1_2::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    double tmass201 = Flavour((kf_code)(23)).Mass();
    CE.TChannelMomenta(p[0],p[1],p[2],p[3],s2,s3,tmass201,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[0],ran[1]);
    Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C1_2::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &ran) {
    double tmass201 = Flavour((kf_code)(23)).Mass();
    double wt = CE.TChannelWeight(p[0],p[1],p[2],p[3],tmass201,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[0],ran[1]);
    if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,2*3.-4.);
    Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
    return 1.0/wt;
}
