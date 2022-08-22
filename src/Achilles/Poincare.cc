#include "Achilles/Poincare.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

Poincare::Poincare(const FourVector &v,const double &rsq):
  m_type(1), m_l(v), m_rsq(rsq>0.?rsq:v.M()) {}

Poincare::Poincare(const FourVector &v1,const FourVector &v2,int mode):
  m_type(mode?3:2), m_l(1.,0.,0.,0.), m_rsq(1.)
{
  if (m_type==3) {
    m_l=v1;
    m_t=v2;
    return;
  }
  FourVector b(v2.Vec3()/v2.P(),0.);
  m_l=FourVector(v1.Vec3()/v1.P(),0.);
  m_t=b+m_l*(m_l*b);
  double mt(m_t.P2());
  if (mt!=0.) m_t*=1./sqrt(mt);
  size_t l[3]{1,2,3};
  double ml[4]={0.,std::abs(m_l[1]),
        std::abs(m_l[2]),std::abs(m_l[3])};
  if (ml[l[2]]>ml[l[1]]) std::swap<size_t>(l[1],l[2]);
  if (ml[l[1]]>ml[l[0]]) std::swap<size_t>(l[0],l[1]);
  if (ml[l[2]]>ml[l[1]]) std::swap<size_t>(l[1],l[2]);
  double tdp(m_t[l[1]]*m_l[l[1]]+m_t[l[2]]*m_l[l[2]]);
  if (tdp!=0.) m_t[l[0]]=-tdp/m_l[l[0]];
  if (m_t.P2()==0.) m_t[l[1]]=1.;
  m_omct=m_l.SmallOMCT(b);
  m_st=-m_t*b;
}

void Poincare::Boost(FourVector &v) const
{
  double lv(m_l[1]*v[1]+m_l[2]*v[2]+m_l[3]*v[3]);
  double v0((m_l[0]*v[0]-lv)/m_rsq);
  double c1((v[0]+v0)/(m_rsq+m_l[0]));
  v=FourVector(v.Vec3()-c1*m_l.Vec3(),v0);
}

void Poincare::BoostBack(FourVector &v) const
{
  double lv(m_l[1]*v[1]+m_l[2]*v[2]+m_l[3]*v[3]);
  double v0((m_l[0]*v[0]+lv)/m_rsq);
  double c1((v[0]+v0)/(m_rsq+m_l[0]));
  v=FourVector(v.Vec3()+c1*m_l.Vec3(),v0);
}

void Poincare::Rotate(FourVector &v) const
{
  double vx(-m_l*v), vy(-m_t*v);
  v-=(m_omct*vx+m_st*vy)*m_l;
  v-=(-m_st*vx+m_omct*vy)*m_t;
}

void Poincare::RotateBack(FourVector &v) const
{
  double vx(-m_l*v), vy(-m_t*v);
  v-=(m_omct*vx-m_st*vy)*m_l;
  v-=(m_st*vx+m_omct*vy)*m_t;
}

void Poincare::Lambda(FourVector &v) const
{
  double m2=v.M2();
  v=v-2.0*(v*(m_l+m_t))/(m_l+m_t).M2()*(m_l+m_t)
    +2.0*(v*m_l)/m_l.M2()*m_t;
  v[0]=Sign(v[0])*sqrt(v.P2()+m2);
}

void Poincare::LambdaBack(FourVector &v) const
{
  double m2=v.M2();
  v=v-2.0*(v*(m_l+m_t))/(m_l+m_t).M2()*(m_l+m_t)
    +2.0*(v*m_t)/m_t.M2()*m_l;
  v[0]=Sign(v[0])*sqrt(v.P2()+m2);
}

void Poincare::Invert() 
{
  if (m_type==3) { std::swap<FourVector>(m_l,m_t); return; }
  if (m_type==2) { m_st=-m_st; return; }
  for (size_t i(1);i<4;++i) m_l[i]=-m_l[i];
}

FourVector Poincare_Sequence::operator*(const FourVector &p) const
{
  FourVector np(p);
  for(const auto &transform : *this) np=transform*np;
  return np;
}

void Poincare_Sequence::Invert()
{
  Poincare_Sequence copy(*this);
  reverse_iterator cit(copy.rbegin());
  for (iterator pit(begin());pit!=end();++pit,++cit) {
    cit->Invert();
    *pit=*cit;
  }
}
