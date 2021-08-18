#ifndef PLUGINS_CHANNELS1_HH
#define PLUGINS_CHANNELS1_HH

// Channel_Generator_UniV
#include "PHASIC++/Channels/Channel_Elements.H"
#include "nuchic/Mapper.hh"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {

class Channels1 : public nuchic::Mapper<Vec4D> {
    protected:
        constexpr static double m_amct=1.,m_alpha=0.7,m_ctmax=1.0,m_ctmin=-1.;
        double s2, s3;
    public:
        Channels1(double _s2, double _s3) 
              : s2{std::move(_s2)}, s3{std::move(_s3)} {}

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double>&) const override = 0;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double>&) const override = 0;
        size_t NDims() const override { return 2; }
};

class C1_0 : public Channels1 {
    public:
        C1_0(double _s2, double _s3) : Channels1(_s2, _s3) {}
        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C1_1 : public Channels1 {
    public:
        C1_1(double _s2, double _s3) : Channels1(_s2, _s3) {}
        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

class C1_2 : public Channels1 {
    public:
        C1_2(double _s2, double _s3) : Channels1(_s2, _s3) {}
        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) const override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) const override;
};

}

#endif
