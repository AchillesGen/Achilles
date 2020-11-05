#ifndef HARD_SCATTERING_HH
#define HARD_SCATTERING_HH

#include <utility>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "nuchic/Beams.hh"
#include "nuchic/RunModes.hh"
#include "nuchic/Histogram.hh"

namespace nuchic {

class FourVector;
class Particle;
class Nucleus;

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
        HardScattering(std::shared_ptr<Beam>, std::shared_ptr<Nucleus>, RunMode);
        HardScattering(const HardScattering&) = default;
        HardScattering(HardScattering&&) = default;
        HardScattering& operator=(const HardScattering&) = default;
        HardScattering& operator=(HardScattering&&) = default;
        virtual ~HardScattering() = default;

        // Validation information
        virtual HardScatteringType ScatteringType() const = 0;

        // Calculation details
        virtual double CrossSection(const Particles&) const = 0;

        // Number of phase space variables
        int NVariables() const { return LeptonVariables() + HadronVariables(); }
        virtual int LeptonVariables() const;
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
        const std::shared_ptr<Beam>& GetBeam() const { return m_leptonBeam; }

        // Special phase space routines
        void SetScatteringAngle(double angle) { m_angle = angle; }

        // Test Phasespace
        void SetHist(bool fill) { m_fill = fill; }
        const Histogram& GetHist() const { return hist; }
        double Test(const std::vector<double>&, const double&);

    protected:
        // Phase space factors
        double dp;
        static constexpr int nNucleonTypes = 2;
        static constexpr double dCos = 2;
        static constexpr double dPhi = 2*M_PI;

    private:
        bool m_fill{false};
        Histogram hist{500, 0, 1000, "Test"};
        std::shared_ptr<Beam> m_leptonBeam;
        std::shared_ptr<Nucleus> m_nuc;

        // TODO: Move to the driver class
        RunMode m_mode;
        double m_angle{};
};

class Quasielastic : public HardScattering {
    public:
        Quasielastic(std::shared_ptr<Beam> beam, std::shared_ptr<Nucleus> nuc,
                RunMode mode) : HardScattering(beam, nuc, mode) {}
        HardScatteringType ScatteringType() const override { 
            return HardScatteringType::Quasielastic;
        }

        double CrossSection(const Particles&) const override = 0;
};

class QESpectral : public Quasielastic {
    public:
        QESpectral(std::shared_ptr<Beam> beam, std::shared_ptr<Nucleus> nuc,
                RunMode mode) : Quasielastic(beam, nuc, mode) {}
        QESpectral(const YAML::Node&, std::shared_ptr<Beam> beam, std::shared_ptr<Nucleus> nuc);

        int HadronVariables() const override { return 4; }
        Particles GenerateHadrons(const std::vector<double>&, const FourVector&) const override;
        double HadronWeight(const Particles&) const override;

        double CrossSection(const Particles&) const override {
            return 10;
        }
};

class FQESpectral : public QESpectral {
    public:
        FQESpectral(const YAML::Node&, std::shared_ptr<Beam>,
                    std::shared_ptr<Nucleus>, RunMode);

        double CrossSection(const Particles&) const override;
        static std::string GetName() { return "QESpectral"; }
        static std::unique_ptr<HardScattering> Create(const YAML::Node &node,
                std::shared_ptr<Beam> beam, std::shared_ptr<Nucleus> nucleus,
                RunMode mode) {
            return std::make_unique<FQESpectral>(node, beam, nucleus, mode);
        }

    private:
        static bool registered;
};

class QEGlobalFermiGas : public Quasielastic {
    public:
        QEGlobalFermiGas(std::shared_ptr<Beam> beam, std::shared_ptr<Nucleus> nuc,
                RunMode mode) : Quasielastic(beam, nuc, mode) {}
        QEGlobalFermiGas(const YAML::Node&, std::shared_ptr<Beam> beam,
                         std::shared_ptr<Nucleus> nuc);

        int HadronVariables() const override { return 3; }
        Particles GenerateHadrons(const std::vector<double>&, const FourVector&) const override;
        double HadronWeight(const Particles&) const override;

        double CrossSection(const Particles&) const override {
            return 10;
        }
};

class FQEGlobalFermiGas : public QEGlobalFermiGas {
    public:
        FQEGlobalFermiGas(const YAML::Node&, std::shared_ptr<Beam>,
                          std::shared_ptr<Nucleus>, RunMode);

        double CrossSection(const Particles&) const override;
        static std::string GetName() { return "QEGlobalFermiGas"; }
        static std::unique_ptr<HardScattering> Create(const YAML::Node &node,
                std::shared_ptr<Beam> beam, std::shared_ptr<Nucleus> nucleus,
                RunMode mode) {
            return std::make_unique<FQEGlobalFermiGas>(node, beam, nucleus, mode);
        }

    private:
        static bool registered;
};

class DIS : HardScattering {
    public:
        HardScatteringType ScatteringType() const override {
            return HardScatteringType::DeepInelastic;
        }
};

}

#endif
