#ifndef HADRONIC_MAPPER_HH
#define HADRONIC_MAPPER_HH

#include <cmath>

#include "nuchic/Mapper.hh"

namespace nuchic {

class FourVector;

class HadronicMapper : public Mapper<FourVector> {
    public:
        HadronicMapper() = default;

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override = 0;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override = 0;
        size_t NDims() const override = 0;
};

class QESpectralMapper : public HadronicMapper {
    public:
        QESpectralMapper() = default;

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override;

    private:
        static constexpr double dCos = 2;
        static constexpr double dPhi = 2*M_PI;
        static constexpr double dp = 1000;
};

}


#endif
