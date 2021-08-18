#include "plugins/Channels3.hh"
#include "ATOOLS/Phys/Flavour.H"
#include <cmath>
#include "plugins/PrintVec.hh"

using namespace PHASIC;
using namespace ATOOLS;

void Channels3::GenerateNuclearPoint(std::vector<Vec4D> &p, const std::vector<double> &ran) const {
    double s234_max = sqr((p[0]+p[1]).Mass()-sqrt(s5));
    double s234_min = pow(sqrt(s3) + sqrt(s4) + sqrt(s2), 2);
    double s234 = CE.MasslessPropMomenta(0.5,s234_min,s234_max,ran[5]);
    Vec4D p234;
    CE.TChannelMomenta(p[0],p[1],p[5],p234,s5,s234,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[6],ran[7]);
}

double Channels3::GenerateNuclearWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) const {
    double wt = 1.;

    double s234_max = sqr((p[0]+p[1]).Mass()-sqrt(s5));
    double s234_min = pow(sqrt(s3) + sqrt(s4) + sqrt(s2), 2);
    wt *= CE.MasslessPropWeight(0.5,s234_min,s234_max,dabs((p[2]+p[3]+p[4]).Abs2()), rans[5]);
    wt *= CE.TChannelWeight(p[0],p[1],p[5],p[2]+p[3]+p[4],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,rans[6],rans[7]);

    return wt;
}

void C3_0::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> & ran) const
{
  GenerateNuclearPoint(p, ran);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Vec4D p34;
  double s34 = CE.MasslessPropMomenta(m_salpha,s34_min,s34_max,ran[0]);
  CE.Isotropic2Momenta(p[0]+p[1]-p[5],s2,s34,p[2],p34,ran[1],ran[2]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[3],ran[4]);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_0::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> & rans) const
{
  double wt = GenerateNuclearWeight(p, rans);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  wt *= CE.MasslessPropWeight(m_salpha,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  wt *= CE.Isotropic2Weight(p[2],p[3]+p[4],rans[1],rans[2],-1,1);
  wt *= CE.Isotropic2Weight(p[3],p[4],rans[3],rans[4],-1,1);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*4.-4.);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
  return 1/wt;
}

void C3_1::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> & ran) const
{ 
  GenerateNuclearPoint(p, ran);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Flavour fl34 = Flavour((kf_code)(23));
  Vec4D  p34;
  double s34 = CE.MassivePropMomenta(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,ran[0]);
  CE.Isotropic2Momenta(p[0]+p[1]-p[5],s2,s34,p[2],p34,ran[1],ran[2]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[3],ran[4]);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_1::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> & rans) const
{
  double wt = GenerateNuclearWeight(p, rans);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Flavour fl34 = Flavour((kf_code)(23));
  wt *= CE.MassivePropWeight(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  wt *= CE.Isotropic2Weight(p[2],p[3]+p[4],rans[1],rans[2],-1,1);
  wt *= CE.Isotropic2Weight(p[3],p[4],rans[3],rans[4],-1,1);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*4.-4.);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
  return 1/wt;
}

void C3_2::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> & ran) const
{
  GenerateNuclearPoint(p, ran);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Vec4D  p34;
  double s34 = CE.MasslessPropMomenta(m_salpha,s34_min,s34_max,ran[0]);
  CE.TChannelMomenta(p[0]-p[5],p[1],p[2],p34,s2,s34,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[1],ran[2]);
  CE.TChannelMomenta(p[0]-p[5],p[1]-p[2],p[4],p[3],s4,s3,0.,m_alpha,1.,-1.,m_amct,0,ran[3],ran[4]);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_2::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> & rans) const
{
  double wt = GenerateNuclearWeight(p, rans);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  wt *= CE.MasslessPropWeight(m_salpha,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  wt *= CE.TChannelWeight(p[0]-p[5],p[1],p[2],p[3]+p[4],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,rans[1],rans[2]);
  wt *= CE.TChannelWeight(p[0]-p[5],p[1]-p[2],p[4],p[3],0.,m_alpha,1.,-1.,m_amct,0,rans[3],rans[4]);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*4.-4.);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
  return 1/wt;
}

void C3_3::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> & ran) const
{
  GenerateNuclearPoint(p, ran);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Vec4D  p34;
  double s34 = CE.MasslessPropMomenta(m_salpha,s34_min,s34_max,ran[0]);
  double tmass201 = Flavour((kf_code)(23)).Mass();
  CE.TChannelMomenta(p[0]-p[5],p[1],p[2],p34,s2,s34,tmass201,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[1],ran[2]);
  CE.TChannelMomenta(p[0]-p[5],p[1]-p[2],p[4],p[3],s4,s3,0.,m_alpha,1.,-1.,m_amct,0,ran[3],ran[4]);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_3::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> & rans) const
{
  double wt = GenerateNuclearWeight(p, rans);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  wt *= CE.MasslessPropWeight(m_salpha,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  double tmass201 = Flavour((kf_code)(23)).Mass();
  wt *= CE.TChannelWeight(p[0]-p[5],p[1],p[2],p[3]+p[4],tmass201,m_alpha,m_ctmax,m_ctmin,m_amct,0,rans[1],rans[2]);
  wt *= CE.TChannelWeight(p[0]-p[5],p[1]-p[2],p[4],p[3],0.,m_alpha,1.,-1.,m_amct,0,rans[3],rans[4]);
  if (wt!=0.) wt = 1./wt/pow(2.*M_PI,3*4.-4.);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
  return 1/wt;
}

void C3_4::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> & ran) const
{
  GenerateNuclearPoint(p, ran);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Vec4D  p34;
  double s34 = CE.MasslessPropMomenta(m_salpha,s34_min,s34_max,ran[0]);
  CE.TChannelMomenta(p[0]-p[5],p[1],p[2],p34,s2,s34,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[1],ran[2]);
  CE.TChannelMomenta(p[0]-p[5],p[1]-p[2],p[3],p[4],s3,s4,0.,m_alpha,1.,-1.,m_amct,0,ran[3],ran[4]);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_4::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> & rans) const
{
  double wt = GenerateNuclearWeight(p, rans);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  wt *= CE.MasslessPropWeight(m_salpha,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  wt *= CE.TChannelWeight(p[0]-p[5],p[1],p[2],p[3]+p[4],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,rans[1],rans[2]);
  wt *= CE.TChannelWeight(p[0]-p[5],p[1]-p[2],p[3],p[4],0.,m_alpha,1.,-1.,m_amct,0,rans[3],rans[4]);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*4.-4.);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
  return 1/wt;
}

void C3_5::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> & ran) const
{
  GenerateNuclearPoint(p, ran);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Vec4D  p34;
  double s34 = CE.MasslessPropMomenta(m_salpha,s34_min,s34_max,ran[0]);
  double tmass201 = Flavour((kf_code)(23)).Mass();
  CE.TChannelMomenta(p[0]-p[5],p[1],p[2],p34,s2,s34,tmass201,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[1],ran[2]);
  CE.TChannelMomenta(p[0]-p[5],p[1]-p[2],p[3],p[4],s3,s4,0.,m_alpha,1.,-1.,m_amct,0,ran[3],ran[4]);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_5::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> & rans) const
{
  double wt = GenerateNuclearWeight(p, rans);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  wt *= CE.MasslessPropWeight(m_salpha,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  double tmass201 = Flavour((kf_code)(23)).Mass();
  wt *= CE.TChannelWeight(p[0]-p[5],p[1],p[2],p[3]+p[4],tmass201,m_alpha,m_ctmax,m_ctmin,m_amct,0,rans[1],rans[2]);
  wt *= CE.TChannelWeight(p[0]-p[5],p[1]-p[2],p[3],p[4],0.,m_alpha,1.,-1.,m_amct,0,rans[3],rans[4]);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*4.-4.);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
  return 1/wt;
}

void C3_6::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> & ran) const
{
  GenerateNuclearPoint(p, ran);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Vec4D  p34;
  double s34 = CE.MasslessPropMomenta(m_salpha,s34_min,s34_max,ran[0]);
  CE.TChannelMomenta(p[0]-p[5],p[1],p34,p[2],s34,s2,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[1],ran[2]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[3],ran[4]);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_6::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> & rans) const
{
  double wt = GenerateNuclearWeight(p, rans);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  wt *= CE.MasslessPropWeight(m_salpha,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  wt *= CE.TChannelWeight(p[0]-p[5],p[1],p[3]+p[4],p[2],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,rans[1],rans[2]);
  wt *= CE.Isotropic2Weight(p[3],p[4],rans[3],rans[4],-1,1);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*4.-4.);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
  return 1/wt;
}

void C3_7::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> & ran) const
{
  GenerateNuclearPoint(p, ran);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Flavour fl34 = Flavour((kf_code)(23));
  Vec4D  p34;
  double s34 = CE.MassivePropMomenta(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,ran[0]);
  CE.TChannelMomenta(p[0]-p[5],p[1],p34,p[2],s34,s2,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[1],ran[2]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[3],ran[4]);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_7::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) const
{
  double wt = GenerateNuclearWeight(p, rans);
  double s34_max = sqr((p[0]+p[1]-p[5]).Mass()-sqrt(s2));
  // TODO: set minimum from cuts
  // double s34_min = cuts->GetscutAmegic(std::string("34"));
  double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
  Flavour fl34 = Flavour((kf_code)(23));
  wt *= CE.MassivePropWeight(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,dabs((p[3]+p[4]).Abs2()),rans[0]);
  wt *= CE.TChannelWeight(p[0]-p[5],p[1],p[3]+p[4],p[2],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,rans[1],rans[2]);
  wt *= CE.Isotropic2Weight(p[3],p[4],rans[3],rans[4],-1,1);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,3*4.-4.);
  Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
  return 1/wt;
}
