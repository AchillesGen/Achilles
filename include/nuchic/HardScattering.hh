#ifndef HARD_SCATTERING_HH
#define HARD_SCATTERING_HH

#include <utility>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "nuchic/Beams.hh"
#include "nuchic/Nucleus.hh"

#include "nuchic/Histogram.hh"

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
        HardScattering(Beam leptonBeam, std::shared_ptr<Nucleus> nuc) 
            : m_leptonBeam{std::move(leptonBeam)}, m_nuc{std::move(nuc)} {};
 //       HardScattering(const HardScattering&) = default;
 //       HardScattering(HardScattering&&) = default;
 //       HardScattering& operator=(const HardScattering&) = default;
 //       HardScattering& operator=(HardScattering&&) = default;
        virtual ~HardScattering() = default;

        // Validation information
        virtual HardScatteringType ScatteringType() const = 0;

        // Calculation details
        virtual double CrossSection(const Particles&) const = 0;

        // Number of phase space variables
        int NVariables() const { return LeptonVariables() + HadronVariables(); }
        virtual int LeptonVariables() const = 0;
        virtual int HadronVariables() const = 0;

        // Generate phase space momentums
        virtual Particles GeneratePhaseSpace(const std::vector<double>&) const;
        virtual Particles GenerateLeptons(const std::vector<double>&) const;
        virtual Particles GenerateHadrons(const std::vector<double>&, const FourVector&) const = 0;

        // Phase space weight
        virtual double PhaseSpaceWeight(const Particles&) const;
        virtual double LeptonWeight(const Particles&) const;
        virtual double HadronWeight(const Particles&) const = 0;

        // Nucleus information
        const std::shared_ptr<Nucleus>& GetNucleus() const { return m_nuc; }
        std::shared_ptr<Nucleus>& GetNucleus() { return m_nuc; }

        // Beam information
        const Beam& GetBeam() const { return m_leptonBeam; }

        // Test Phasespace
        void SetHist(bool fill) { m_fill = fill; }
        const Histogram& GetHist() const { return hist; }
        double Test(const std::vector<double>&, const double&);

    private:
        bool m_fill{false};
        Histogram hist{500, 0, 1000, "Test"};
        Beam m_leptonBeam;
        std::shared_ptr<Nucleus> m_nuc;

};

class Quasielastic : public HardScattering {
    public:
        Quasielastic(Beam beam, std::shared_ptr<Nucleus> nuc) : HardScattering(beam, nuc) {}
        HardScatteringType ScatteringType() const override { 
            return HardScatteringType::Quasielastic;
        }

        int LeptonVariables() const override { return 3 + GetBeam().NVariables(); }

        double CrossSection(const Particles&) const override = 0;
};

class QESpectral : public Quasielastic {
    public:
        QESpectral(Beam beam, std::shared_ptr<Nucleus> nuc) : Quasielastic(beam, nuc) {}

        int HadronVariables() const override { return 4; }
        Particles GenerateHadrons(const std::vector<double>&, const FourVector&) const override;
        double HadronWeight(const Particles&) const override;

        double CrossSection(const Particles&) const override {
            return 10;
        }
};

class DIS : HardScattering {
    public:
        HardScatteringType ScatteringType() const override {
            return HardScatteringType::DeepInelastic;
        }
};

}

#endif
