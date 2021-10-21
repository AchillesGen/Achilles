#ifndef HADRONIC_MAPPER_HH
#define HADRONIC_MAPPER_HH

#include <cmath>

#include "nuchic/Mapper.hh"
#include "nuchic/PhaseSpaceFactory.hh"

namespace nuchic {

class FourVector;

class HadronicBeamMapper : public Mapper<FourVector> {
    public:
        HadronicBeamMapper(size_t idx, std::string name) 
            : m_idx{std::move(idx)}, m_name{std::move(name)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override = 0;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override = 0;
        size_t NDims() const override = 0;
        YAML::Node ToYAML() const override {
            YAML::Node result;
            result["Name"] = m_name;
            result["Idx"] = m_idx;
            return result;
        }
        static std::string Name() { return "Hadronic Initial State"; }

    protected:
        size_t HadronIdx() const { return m_idx; }

    private:
        size_t m_idx;
        std::string m_name;
};

class QESpectralMapper : public HadronicBeamMapper, RegistrablePS<HadronicBeamMapper, QESpectralMapper, size_t> {
    public:
        QESpectralMapper(size_t idx) : HadronicBeamMapper(idx, Name()) {}
        static std::string Name() { return "QESpectral"; }
        static std::unique_ptr<HadronicBeamMapper> Construct(const size_t &idx) {
            return std::make_unique<QESpectralMapper>(idx);
        }

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override { return 4; }

    private:
        static constexpr double dCos = 2;
        static constexpr double dPhi = 2*M_PI;
        static constexpr double dp = 800;
};

}


#endif
