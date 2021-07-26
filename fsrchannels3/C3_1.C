// Channel_Generator_UniV
#include "PHASIC++/Channels/Channel_Elements.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C3_1 : public Single_Channel {
    constexpr static double m_salpha=0.9,m_thexp=0.9;
  public:
    void   GenerateWeight(Vec4D *,Cut_Data *, std::vector<double> &);
    void   GeneratePoint(Vec4D *,Cut_Data *,const std::vector<double> &);
  };
}

void C3_1::GeneratePoint(Vec4D * p,Cut_Data * cuts,const std::vector<double> & ran)
{
  double s34_max = sqr((p[0]+p[1]).Mass()-sqrt(s2));
  double s34_min = cuts->GetscutAmegic(std::string("34"));
  Flavour fl34 = Flavour((kf_code)(23));
  Vec4D  p34;
  double s34 = CE.MassivePropMomenta(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,ran[0]);
  CE.Isotropic2Momenta(p[0]+p[1],s2,s34,p[2],p34,ran[1],ran[2]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[3],ran[4]);
}

void C3_1::GenerateWeight(Vec4D* p,Cut_Data * cuts, std::vector<double> & rans)
{
  double wt = 1.;
  double s34_max = sqr((p[0]+p[1]).Mass()-sqrt(s2));
  double s34_min = cuts->GetscutAmegic(std::string("34"));
  Flavour fl34 = Flavour((kf_code)(23));
  Vec4D  p34 = p[3]+p[4];
  wt *= CE.MassivePropWeight(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  wt *= CE.Isotropic2Weight(p[2],p[3]+p[4],rans[1],rans[2]);
  wt *= CE.Isotropic2Weight(p[3],p[4],rans[3],rans[4]);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*3.-4.);

  weight = wt;
}
