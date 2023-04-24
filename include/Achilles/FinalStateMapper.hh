#ifndef LEPTONIC_MAPPER_HH
#define LEPTONIC_MAPPER_HH

#include "Achilles/Mapper.hh"
#include "Achilles/PhaseSpaceFactory.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "ATOOLS/Math/Vector.H"
#pragma GCC diagnostic pop
#endif // ACHILLES_SHERPA_INTERFACE

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

#ifdef ACHILLES_SHERPA_INTERFACE
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
#endif // ACHILLES_SHERPA_INTERFACE

}

#endif
