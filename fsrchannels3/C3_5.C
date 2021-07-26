// Channel_Generator_UniV
#include "PHASIC++/Channels/Channel_Elements.H"
#include "nuchic/Mapper.hh"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C3_5 : public nuchic::Mapper<Vec4D> {
    constexpr static double m_salpha=0.9,m_thexp=0.9;
    constexpr static double m_amct=1.,m_alpha=0.9,m_ctmax=1.,m_ctmin=-1.;
  public:
    double   GenerateWeight(std::vector<Vec4D> *, std::vector<double> &) const override;
    void   GeneratePoint(std::vector<Vec4D> *,const std::vector<double> &) const override;
  };
}

void C3_5::GeneratePoint(std::vector<Vec4D> * p,const std::vector<double> & ran) const
{
  double s34_max = sqr((p[0]+p[1]).Mass()-sqrt(s2));
  double s34_min = cuts->GetscutAmegic(std::string("34"));
  Vec4D  p34;
  double s34 = CE.MasslessPropMomenta(m_salpha,s34_min,s34_max,ran[0]);
  double tmass201 = Flavour((kf_code)(23)).Mass();
  CE.TChannelMomenta(p[0],p[1],p[2],p34,s2,s34,tmass201,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[1],ran[2]);
  CE.TChannelMomenta(p[0]-p[2],p[1],p[3],p[4],s3,s4,0.,m_alpha,1.,-1.,m_amct,0,ran[3],ran[4]);
}

double C3_5::GenerateWeight(std::vector<Vec4D> * p, std::vector<double> & rans) const
{
  double wt = 1.;
  double s34_max = sqr((p[0]+p[1]).Mass()-sqrt(s2));
  double s34_min = cuts->GetscutAmegic(std::string("34"));
  wt *= CE.MasslessPropWeight(m_salpha,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  double tmass201 = Flavour((kf_code)(23)).Mass();
  wt *= CE.TChannelWeight(p[0],p[1],p[2],p[3]+p[4],tmass201,m_alpha,m_ctmax,m_ctmin,m_amct,0,rans[1],rans[2]);
  wt *= CE.TChannelWeight(p[0]-p[2],p[1],p[3],p[4],0.,m_alpha,1.,-1.,m_amct,0,rans[3],rans[4]);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*3.-4.);

  weight = wt;
  return weight
}
