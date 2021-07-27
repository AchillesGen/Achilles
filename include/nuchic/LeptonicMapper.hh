#ifndef LEPTONIC_MAPPER_HH
#define LEPTONIC_MAPPER_HH

#include "nuchic/Mapper.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "ATOOLS/Math/Vector.H"
#pragma GCC diagnostic pop

namespace nuchic {

class FourVector;

class LeptonicMapper : public Mapper<FourVector> {
    public:
        LeptonicMapper() = default;

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override { return 0; }
};

class SherpaMapper : public Mapper<FourVector> {
    public:
        SherpaMapper(Mapper_ptr<ATOOLS::Vec4D> _mapper) 
            : sherpa_mapper{std::move(_mapper)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override { return sherpa_mapper -> NDims(); }

    private:
        Mapper_ptr<ATOOLS::Vec4D> sherpa_mapper;
};

}

#endif
