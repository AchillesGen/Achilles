#ifndef HADRONIC_MAPPER_HH
#define HADRONIC_MAPPER_HH

#include <cmath>

#include "nuchic/Mapper.hh"

namespace nuchic {

class FourVector;

class HadronicBeamMapper : public Mapper<FourVector> {
    public:
        HadronicBeamMapper(size_t idx) : m_idx{std::move(idx)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override = 0;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override = 0;
        size_t NDims() const override = 0;

    protected:
        size_t HadronIdx() const { return m_idx; }

    private:
        size_t m_idx;
};

class QESpectralMapper : public HadronicBeamMapper {
    public:
        QESpectralMapper(size_t idx) : HadronicBeamMapper(idx) {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override { return 4; }

    private:
        static constexpr double dCos = 2;
        static constexpr double dPhi = 2*M_PI;
        // static constexpr double dp = 800;
};

}


#endif
