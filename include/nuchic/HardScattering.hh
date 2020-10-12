#ifndef HARD_SCATTERING_HH
#define HARD_SCATTERING_HH

#include <vector>

namespace nuchic {

class FourVector;
class Particle;

using Particles = std::vector<Particle>;

enum class HardScatteringType {
    Quasielastic,
    MesonExchangeCurrent,
    Interference_QE_MEC,
    Resonance,
    ShallowInelastic,
    DeepInelastic
};

class HardScattering {
    public:
        HardScattering() = default;
        HardScattering(const HardScattering&) = default;
        HardScattering(HardScattering&&) = default;
        HardScattering& operator=(const HardScattering&) = default;
        HardScattering& operator=(HardScattering&&) = default;
        virtual ~HardScattering() = default;

        // Validation information
        virtual HardScatteringType ScatteringType() const = 0;

        // Calculation details
        virtual size_t PhaseSpaceSize() const = 0;
        virtual Particles GeneratePhaseSpace(const std::vector<double>&) const = 0;
        virtual double CrossSection(const Particles&) const = 0;
        virtual bool HasInitialState() const { return false; }
        virtual FourVector GenInitialState(const std::vector<FourVector>&) const;

};

class Quasielastic : HardScattering {
    public:
        HardScatteringType ScatteringType() const override { 
            return HardScatteringType::Quasielastic;
        }
};

class QESpectral : Quasielastic {
    public:
       size_t PhaseSpaceSize() const override { return 10; }
       Particles GeneratePhaseSpace(const std::vector<double>&) const override;
       double CrossSection(const Particles&) const override;

    private:
       Tensor<2> HadronicTensor(const Particles&) const;
};

class DIS : HardScattering {
    public:
        HardScatteringType ScatteringType() const override {
            return HardScatteringType::DeepInelastic;
        }

        size_t PhaseSpaceSize() const override { return 11; }
};

}

#endif
