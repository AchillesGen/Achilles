#ifndef PLUGINS_CHANNELS3_HH
#define PLUGINS_CHANNELS3_HH

// Channel_Generator_UniV
#include "PHASIC++/Channels/Channel_Elements.H"
#include "nuchic/Mapper.hh"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {

class Channels3 : public nuchic::Mapper<Vec4D> {
protected:
  constexpr static double m_salpha=0.9,m_thexp=0.9;
  constexpr static double m_amct=1.,m_alpha=0.9,m_ctmax=1.,m_ctmin=-1.;
  double s2, s3, s4, s5;

  void GenerateNuclearPoint(std::vector<Vec4D>&, const std::vector<double>&) const;
  double GenerateNuclearWeight(const std::vector<Vec4D>&, std::vector<double>&) const;
public:
  Channels3(double _s2, double _s3, double _s4, double _s5) 
      : s2{std::move(_s2)}, s3{std::move(_s3)}, s4{std::move(_s4)}, s5{std::move(_s5)} {}

  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double>&) const override = 0;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double>&) const override = 0;
  size_t NDims() const override { return 5; }
};

class C3_0 : public Channels3 {
public:
  C3_0(double _s2, double _s3, double _s4, double _s5) : Channels3(_s2, _s3, _s4, _s5) {}
  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C3_1 : public Channels3 {
public:
  C3_1(double _s2, double _s3, double _s4, double _s5) : Channels3(_s2, _s3, _s4, _s5) {}
  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C3_2 : public Channels3 {
public:
  C3_2(double _s2, double _s3, double _s4, double _s5) : Channels3(_s2, _s3, _s4, _s5) {}
  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C3_3 : public Channels3 {
public:
  C3_3(double _s2, double _s3, double _s4, double _s5) : Channels3(_s2, _s3, _s4, _s5) {}
  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C3_4 : public Channels3 {
public:
  C3_4(double _s2, double _s3, double _s4, double _s5) : Channels3(_s2, _s3, _s4, _s5) {}
  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C3_5 : public Channels3 {
public:
  C3_5(double _s2, double _s3, double _s4, double _s5) : Channels3(_s2, _s3, _s4, _s5) {}
  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C3_6 : public Channels3 {
public:
  C3_6(double _s2, double _s3, double _s4, double _s5) : Channels3(_s2, _s3, _s4, _s5) {}
  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C3_7 : public Channels3 {
public:
  C3_7(double _s2, double _s3, double _s4, double _s5) : Channels3(_s2, _s3, _s4, _s5) {}
  double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
  void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

}

#endif
