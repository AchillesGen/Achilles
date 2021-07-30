#ifndef LEPTONIC_MAPPER_HH
#define LEPTONIC_MAPPER_HH

#include "nuchic/Mapper.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "ATOOLS/Math/Vector.H"
#pragma GCC diagnostic pop

namespace nuchic {

class FourVector;

class FinalStateMapper : public Mapper<FourVector> {
    public:
        FinalStateMapper(size_t _nout) : nout{std::move(_nout)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override = 0;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override = 0;
        size_t NDims() const override { return 3*nout - 4; }

    private:
        size_t nout;
};

class TwoBodyMapper : public FinalStateMapper {
    public:
        TwoBodyMapper(double _s2, double _s3) : FinalStateMapper(2), s2{std::move(_s2)}, s3{std::move(_s3)} {}
        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;

    private:
        const double s2, s3;
        static constexpr double dCos = 2;
        static constexpr double dPhi = 2*M_PI;
};

class SherpaMapper : public FinalStateMapper {
    public:
        SherpaMapper(size_t _nout, Mapper_ptr<ATOOLS::Vec4D> _mapper) 
            : FinalStateMapper(_nout), sherpa_mapper{std::move(_mapper)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        // size_t NDims() const override { return sherpa_mapper -> NDims(); }

    private:
        Mapper_ptr<ATOOLS::Vec4D> sherpa_mapper;
};

}

#endif
