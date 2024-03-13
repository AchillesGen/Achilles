#ifndef LEPTONIC_MAPPER_HH
#define LEPTONIC_MAPPER_HH

#include "Achilles/Mapper.hh"
#include "Achilles/PhaseSpaceFactory.hh"

#ifdef ENABLE_BSM
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wimplicit-float-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#include "ATOOLS/Math/Vector.H"
#pragma GCC diagnostic pop
#endif // ENABLE_BSM

namespace achilles {

class FourVector;

class FinalStateMapper : public Mapper<FourVector> {
    public:
        FinalStateMapper(size_t _nout) : nout{std::move(_nout)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) override = 0;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) override = 0;
        size_t NDims() const override { return 3*nout - 4; }
        YAML::Node ToYAML() const override = 0;
        static std::string Name() { return "Final State"; }

    private:
        size_t nout;
};

class TwoBodyMapper : public FinalStateMapper,
                             RegistrablePS<FinalStateMapper, TwoBodyMapper, std::vector<double>> {
    public:
        TwoBodyMapper(const std::vector<double> &m) : FinalStateMapper(2), s2{m[0]}, s3{m[1]} {}
        static std::string Name() { return "TwoBody"; }
        static std::unique_ptr<FinalStateMapper> Construct(const std::vector<double> &m) {
            if(m.size() != 2) {
                auto msg = fmt::format("Incorrect number of masses. Expected 2. Got {}", m.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<TwoBodyMapper>(m);
        }

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) override;
        YAML::Node ToYAML() const override {
            YAML::Node result;
            result["Name"] = Name();
            result["Masses"] = std::vector<double>{s2, s3};
            return result;
        }

    private:
        const double s2, s3;
        static constexpr double dCos = 2;
        static constexpr double dPhi = 2*M_PI;
};

class ThreeBodyMapper : public FinalStateMapper,
                             RegistrablePS<FinalStateMapper, ThreeBodyMapper, std::vector<double>> {
    public:
        ThreeBodyMapper(const std::vector<double> &m) : FinalStateMapper(3), s2{m[0]}, s3{m[1]}, s4{m[2]} {}
        static std::string Name() { return "ThreeBody"; }
        static std::unique_ptr<FinalStateMapper> Construct(const std::vector<double> &m) {
            if(m.size() != 3) {
                auto msg = fmt::format("Incorrect number of masses. Expected 3. Got {}", m.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<ThreeBodyMapper>(m);
        }

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) override;
        double MassivePropagatorMomenta(double mass, double width, double smin, double smax, double ran);
        void Isotropic2Momenta(FourVector p, double s1_, double s2_, FourVector &p1, FourVector &p2, double ran1, double ran2, double ctmin, double ctmax); 
        int TChannelMomenta(FourVector p1in, FourVector p2in, FourVector &p1out, FourVector &p2out, double s1out,double s2out, double t_mass, double ctexp, double ctmax, double ctmin, double aminct, int, double ran1, double ran2);
        double SqLam(double s_,double s1_, double s2_);
        double Tj1(double cn,double amcxm,double amcxp,double ran);
        double Hj1(double cn, double amcxm, double amcxp);
        double MassivePropWeight(double mass, double width, int lim, double smin, double smax, double s, double&);
        double TChannelWeight(const FourVector &p1in, const FourVector &p2in, const FourVector &p1out, const FourVector &p2out, double t_mass, double ctexp, double ctmax, double ctmin, double aminct, int, double &ran1, double &ran2);
        double Isotropic2Weight(const FourVector& p1,const FourVector& p2,
                      double& ran1,double& ran2,double ctmin,double ctmax);
        void Boost(int lflag, const FourVector &q, const FourVector &ph, FourVector &p);
        YAML::Node ToYAML() const override {
            YAML::Node result;
            result["Name"] = Name();
            result["Masses"] = std::vector<double>{s2, s3, s4};
            return result;
        }
        void SetGaugeBosonMass(double mass) override { m_mass = mass; }

    private:
        const double s2, s3, s4;
        double m_mass; 

        static constexpr double dCos = 2;
        static constexpr double dPhi = 2*M_PI;
        constexpr static double m_salpha=0.9,m_thexp=0.9;
        constexpr static double m_amct=1.,m_alpha=0.9,m_ctmax=1.,m_ctmin=-1.;
};

#ifdef ENABLE_BSM
class SherpaMapper : public FinalStateMapper {
    public:
        SherpaMapper(size_t _nout, Mapper_ptr<ATOOLS::Vec4D> _mapper) 
            : FinalStateMapper(_nout), sherpa_mapper{std::move(_mapper)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) override;
        // size_t NDims() const override { return sherpa_mapper -> NDims(); }
        YAML::Node ToYAML() const override {
            YAML::Node result;
            result["Name"] = "SherpaMapper";
            result["Sherpa"] = sherpa_mapper -> ToYAML();
            result["Dim"] = NDims();
            return result;
        }

    private:
        Mapper_ptr<ATOOLS::Vec4D> sherpa_mapper;
};
#endif // ENABLE_BSM

}

#endif
