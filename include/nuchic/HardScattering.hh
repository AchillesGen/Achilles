#ifndef HARD_SCATTERING_HH
#define HARD_SCATTERING_HH

#include "plugins/SherpaMEs.hh"
#include <utility>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "nuchic/HardScatteringEnum.hh"
#include "nuchic/Beams.hh"
#include "nuchic/RunModes.hh"
#include "nuchic/Histogram.hh"

extern "C" {
    void Delete();
}

namespace nuchic {

class FourVector;
class Particle;
class Nucleus;
class Event;

struct InitialState;

using Particles = std::vector<Particle>;

class HardScattering {
    public:
        HardScattering(RunMode);
        HardScattering(const HardScattering&) = default;
        HardScattering(HardScattering&&) = default;
        HardScattering& operator=(const HardScattering&) = default;
        HardScattering& operator=(HardScattering&&) = default;
        virtual ~HardScattering() { if(p_sherpa) delete p_sherpa; }

        // Validation information
        virtual HardScatteringType ScatteringType() const = 0;

        // Calculation details
        virtual void CrossSection(Event&) const = 0;

        // Number of phase space variables
        int NVariables() const { return LeptonVariables() + HadronVariables(); }
        virtual int LeptonVariables() const;
        virtual int HadronVariables() const = 0;

        // Generate phase space momentums and return weights
        virtual void GeneratePhaseSpace(const std::vector<double>&, Event&) const;
        virtual void GenerateLeptons(const std::vector<double>&, Event&) const;
        virtual void GenerateHadrons(const std::vector<double>&,
                                     const FourVector&, Event&) const = 0;

        // Select initial state
        virtual bool InitializeEvent(Event&) = 0;
        size_t SelectMatrixElement(Event&) const;

        // Special phase space routines
        void SetScatteringAngle(double angle) { m_angle = angle; }
        void SetFinalLeptonEnergy(double energy) { m_lepton_energy = energy; }

        // Sherpa pointer operations
        void SetSherpa(SherpaMEs *const _sherpa, const std::vector<PID> &pids) { 
            p_sherpa = _sherpa;
            m_pids = pids;
        }
        std::vector<double> LeptonicTensor(const std::vector<FourVector>&,
                                           const double&) const;

        // Test Phasespace
        void SetHist(bool fill) { m_fill = fill; }
        const Histogram& GetHist() const { return hist; }

    protected:
        // Phase space factors
        static constexpr int nNucleonTypes = 2;
        static constexpr double dCos = 2;
        static constexpr double dPhi = 2*M_PI;

    private:
        SherpaMEs *p_sherpa{nullptr};
        std::vector<PID> m_pids;
        bool m_fill{false};
        Histogram hist{500, 0, 1000, "Test"};

        // TODO: Move to the driver class?
        RunMode m_mode;
        double m_angle{}, m_lepton_energy{};
};

class Quasielastic : public HardScattering {
    public:
        Quasielastic(RunMode mode) : HardScattering(mode) {}
        HardScatteringType ScatteringType() const override { 
            return HardScatteringType::Quasielastic;
        }

        void CrossSection(Event&) const override = 0;

        // Select initial state
        bool InitializeEvent(Event&) override;
};

class QESpectral : public Quasielastic {
    public:
        QESpectral(RunMode mode) : Quasielastic(mode) {}
        QESpectral(const YAML::Node&);

        int HadronVariables() const override { return 4; }
        void GenerateHadrons(const std::vector<double>&, const FourVector&, Event&) const override;

        void CrossSection(Event&) const override {}
};

class FQESpectral : public QESpectral {
    public:
        FQESpectral(const YAML::Node&, RunMode);
        ~FQESpectral() override { Delete(); }

        void CrossSection(Event&) const override;
        static std::string GetName() { return "QESpectral"; }
        static std::unique_ptr<HardScattering> Create(const YAML::Node &node,
                RunMode mode) {
            return std::make_unique<FQESpectral>(node, mode);
        }

    private:
        static bool registered;
};

class QEGlobalFermiGas : public Quasielastic {
    public:
        QEGlobalFermiGas(RunMode mode) : Quasielastic(mode) {}
        QEGlobalFermiGas(const YAML::Node&);

        int HadronVariables() const override { return 3; }
        void GenerateHadrons(const std::vector<double>&, const FourVector&, Event&) const override;

        void CrossSection(Event&) const override {}
};

class FQEGlobalFermiGas : public QEGlobalFermiGas {
    public:
        FQEGlobalFermiGas(const YAML::Node&, RunMode);

        void CrossSection(Event&) const override;
        static std::string GetName() { return "QEGlobalFermiGas"; }
        static std::unique_ptr<HardScattering> Create(const YAML::Node &node,
                RunMode mode) {
            return std::make_unique<FQEGlobalFermiGas>(node, mode);
        }

    private:
        static bool registered;
};

class DIS : HardScattering {
    public:
        HardScatteringType ScatteringType() const override {
            return HardScatteringType::DeepInelastic;
        }

        // Select initial state
        bool InitializeEvent(Event&) override;
};

}

#endif
