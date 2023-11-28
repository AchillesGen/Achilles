#pragma once

#include "Achilles/FourVector.hh"

#include <vector>

namespace achilles {

  class Poincare {
  private:

    int m_type;
    FourVector m_l, m_t;
    double m_rsq, m_omct, m_st;

  public:

    Poincare(const FourVector &v=FourVector(0.,0.,0.,0.),
         const double &rsq=-1.);
    Poincare(const FourVector &v1,const FourVector &v2,int mode=0);
 
    void Boost(FourVector &v) const;
    void BoostBack(FourVector &v) const;

    void Rotate(FourVector &v) const;
    void RotateBack(FourVector &v) const;

    void Lambda(FourVector &v) const;
    void LambdaBack(FourVector &v) const;

    void Invert();

    inline FourVector operator*(const FourVector &vin) const
    {
      FourVector v(vin);
      if (m_type==1) Boost(v);
      if (m_type==2) Rotate(v);
      if (m_type==3) Lambda(v);
      return v;
    }

    inline const FourVector &PL() const { return m_l; }
    inline const FourVector &PT() const { return m_t; }

    inline double OMCTheta() const { return m_omct; }
    inline double SinTheta() const { return m_st;   }

  };// end of class Poincare

  class Poincare_Sequence: public std::vector<Poincare> {
  public:

    FourVector operator*(const FourVector &p) const;
    void Invert();

  };// end of class Poincare_Sequence

}// end of namespace apes
