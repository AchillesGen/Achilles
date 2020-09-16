#ifndef METROPOLIS_HH
#define METROPOLIS_HH

#include <array>
#include <vector>
#include <limits>

#include "nuchic/Interactions.hh"
#include "nuchic/ParticleInfo.hh"
#include "nuchic/Rambo.hh"

namespace nuchic {

class Particle;

using Particles = std::vector<Particle>;
using RNG = randutils::mt19937_rng;

enum class ProbeType {
    nucleon, pion
};

//TODO: come up with better names for these
enum class ChargeMode {
    ii, ij, zero
};

enum class ProductionMethod {
    Probe, // produce pions at the probe location
    Target, // produce pions at the target location
    Random // produce pions at a random point between probe and target
};

using InteractionType = std::pair<ProbeType, ChargeMode>;

class MetropolisData {
    public:
        MetropolisData() = default;
        MetropolisData(const MetropolisData&) = default;
        MetropolisData(MetropolisData&&) = default;
        MetropolisData& operator=(const MetropolisData&) = default;
        MetropolisData& operator=(MetropolisData&&) = default;
        virtual ~MetropolisData() = default;

        virtual std::vector<double> XSec(ChargeMode, const Particle&, const Particle&) const = 0;
        virtual Particles Elastic(ChargeMode, const double&,
                                  const std::array<double, 2>&,
                                  const Particle&,
                                  const Particle&) const = 0;
        virtual Particles CEX(ChargeMode, const double&, const std::array<double, 2>&, 
                              const Particle&, const Particle&) const { return {}; }
        virtual Particles SinglePion(ChargeMode, const std::vector<double>&,
                                     const Particle&, const Particle&) = 0;
        virtual Particles DoublePion(ChargeMode, const std::vector<double>&,
                                     const Particle&, const Particle&) = 0;
        virtual bool IsAbsorption(const std::vector<double>&, const double&) const { return false; }
        virtual bool IsInelastic(ChargeMode, const double&, const double&) const { return false; }
        virtual bool IsCEX(ChargeMode, const double&, const double&) const { return false; }
        virtual bool IsSinglePion(ChargeMode, const double&, const double&) const { return false; }

    protected:
        Rambo rambo{2};
        ThreeVector GetPosition(ProductionMethod, const double&,
                                const ThreeVector&, const ThreeVector&) const;
};

class NucleonData : public MetropolisData {
    private:
        static constexpr size_t size = 8;
        using DataTable = std::array<double, size>;

    public:
        NucleonData(ProductionMethod mode) {
            production = mode;
        } 
        NucleonData(const NucleonData&) = default;
        NucleonData(NucleonData&&) = default;
        NucleonData& operator=(const NucleonData&) = default;
        NucleonData& operator=(NucleonData&&) = default;
        ~NucleonData() override = default;

        std::vector<double> XSec(ChargeMode, const Particle&, const Particle&) const override;
        Particles Elastic(ChargeMode, const double&,
                          const std::array<double, 2>&,
                          const Particle&,
                          const Particle&) const override;
        Particles SinglePion(ChargeMode, const std::vector<double>&,
                             const Particle&, const Particle&) override;
        Particles DoublePion(ChargeMode, const std::vector<double>&,
                             const Particle&, const Particle&) override;
        bool IsInelastic(ChargeMode, const double&, const double&) const override;
        bool IsSinglePion(ChargeMode, const double&, const double&) const override;

    private:
        ProductionMethod production;
        std::vector<PID> PionCharge(const double&, const PID&, const PID&) const;
        std::vector<PID> DoublePionCharge(const double&, const PID&, const PID&) const;

        // Data
        static constexpr DataTable energies = {335, 410, 510, 660, 840, 1160, 1780, 3900};
        static constexpr DataTable sigmaii = {24.5, 26.4, 30.4, 41.2, 47.2, 48.0, 44.2, 41.0};
        static constexpr DataTable sigmaij = {33.0, 34.0, 35.1, 36.5, 37.9, 40.2, 42.7, 42.0};
        static constexpr DataTable fiiInel = {0.07, 0.20, 0.31, 0.43, 0.58, 0.65, 0.69, 0.69};
        static constexpr DataTable fijInel = {0.04, 0.07, 0.15, 0.27, 0.37, 0.36, 0.35, 0.35};
        static constexpr DataTable fpi = {1.0, 1.0, 1.0, 1.0, 0.97, 0.80, 0.44, 0.44};
        static constexpr DataTable Aii = {0.1, 0.9, 2.7, 9.0, 14.3, 19.2, std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        static constexpr DataTable Bii = {0, 0, 0, 0, 0, 0, 0, 0};
        static constexpr DataTable Aij = {2.2, 1.8, 2.3, 8.8, 15.0, 29.4, std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        static constexpr DataTable Bij = {-1.0, -1.1, -0.7, -0.2, 0, 0, 0, 0};

        // Pion production parameters
        static constexpr std::array<double, 2> fpi0 = {0.11, 0.43};
        static constexpr std::array<double, 4> f2pi0 = {0.6, 0.2, 0.16, 0.10};
};

class PionData : public MetropolisData {
    private:
        static constexpr size_t size = 8;
        using DataTable = std::array<double, size>;

    public:
        PionData(ProductionMethod mode) {
            production = mode;
        }
        PionData(const PionData&) = default;
        PionData(PionData&&) = default;
        PionData& operator=(const PionData&) = default;
        PionData& operator=(PionData&&) = default;
        ~PionData() override = default;

        std::vector<double> XSec(ChargeMode, const Particle&, const Particle&) const override;
        Particles Elastic(ChargeMode, const double&,
                          const std::array<double, 2>&,
                          const Particle&,
                          const Particle&) const override;
        Particles CEX(ChargeMode, const double&, const std::array<double, 2>&,
                      const Particle&, const Particle&) const override;
        Particles SinglePion(ChargeMode, const std::vector<double>&,
                             const Particle&, const Particle&) override;
        Particles DoublePion(ChargeMode, const std::vector<double>&,
                             const Particle&, const Particle&) override;
        bool IsAbsorption(const std::vector<double>&, const double&) const override; 
        bool IsInelastic(ChargeMode, const double&, const double&) const override;
        bool IsCEX(ChargeMode, const double&, const double&) const override;
        bool IsSinglePion(ChargeMode, const double&, const double&) const override;


    private:
        ProductionMethod production;
        std::vector<PID> PionCharge(ChargeMode, const double&, const PID&, const PID&) const;
        std::vector<PID> DoublePionCharge(const double&, const PID&, const PID&) const;

        double XSecii(const Particle&) const;
        double XSecij(const Particle&) const;
        double XSecAbs(const Particle&) const;

        // Data
        static constexpr DataTable energies = {49, 85, 128, 184, 250, 350, 540, 1300};
        static constexpr DataTable sigmaii = {16, 50, 114, 200, 110, 51, 20, 30};
        static constexpr DataTable sigmaij = {15, 21, 43, 66, 44, 23, 22, 30};
        static constexpr DataTable sigmaabs = {20, 32, 45, 36, 18, 0, 0, 0};
        static constexpr DataTable fiiInel = {0, 0, 0, 0.03, 0.06, 0.16, 0.30, 0.88};
        static constexpr DataTable fiiCEX = {0, 0, 0, 0, 0, 0, 0, 0};
        static constexpr DataTable fijInel = {0.45, 0.57, 0.62, 0.64, 0.62, 0.56, 0.58, 0.94};
        static constexpr DataTable fijCEX = {1.0, 1.0, 1.0, 0.95, 0.89, 0.72, 0.51, 0.06};
        static constexpr DataTable f0Inel = {0.42, 0.36, 0.36, 0.37, 0.40, 0.50, 0.59, 0.94};
        static constexpr DataTable f0CEX = {1.0, 1.0, 1.0, 0.90, 0.84, 0.67, 0.50, 0.05};
        static constexpr DataTable fpi = {1.0, 1.0, 1.0, 1.0, 1.0, 0.98, 0.91, 0.24};
        static constexpr DataTable Aii = {3.2, 2.2, 1.9, 2.2, 2.6, 3.0, 3.0, 3.0};
        static constexpr DataTable Bii = {-1.8, -2.1, -1.5, -0.3, 2.0, 4.0, 4.0, 4.0};
        static constexpr DataTable Aij = {1.1, 1.9, 2.2, 2.2, 2.0, 2.7, 3.0, 3.0};
        static constexpr DataTable Bij = {0.8, 0.7, 0.8, 1.0, 1.4, 2.6, 3.6, 4.0};
        static constexpr DataTable A0 = {3.4, 2.1, 1.9, 2.1, 2.5, 3.0, 3.0, 3.0};
        static constexpr DataTable B0 = {-1.8, -2.0, -1.4, 0, 1.7, 4.0, 4.0, 4.0};

        // Pion production parameters
        static constexpr std::array<double, 2> fpi0 = {0.55, 0.44};
        static constexpr std::array<double, 2> f2pi0 = {0.1875, 0.0625};
};

class Metropolis : public Interactions {
    public:
        Metropolis(const std::string&);

        static std::unique_ptr<Interactions> Create(const std::string &data) {
            return std::make_unique<Metropolis>(data);
        }

        static std::string GetName() { return "Metropolis"; }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        std::vector<double> CrossSections(const Particle&, const Particle&) const override;
        std::vector<Particle> GenerateFinalState(RNG&, const Particle&,
                                                 const Particle&) const override;

        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2> &) const override {
            return {};
        }
    private:
        // Functions
        InteractionType  GetMode(const Particle&, const Particle&) const;
        std::vector<double> XSec(InteractionType, const Particle&, const Particle&) const;
        Particles Absorption(InteractionType, RNG&, const Particle&, const Particle&) const;
        Particles Elastic(InteractionType, RNG&, const double&, const Particle&, const Particle&) const;
        Particles CEX(InteractionType, RNG&, const double&, const Particle&, const Particle&) const;
        Particles SinglePion(InteractionType, RNG&, const Particle&, const Particle&) const;
        Particles DoublePion(InteractionType, RNG&, const Particle&, const Particle&) const;

        bool IsAbsorption(InteractionType, const std::vector<double>&, const double&) const;
        bool IsInelastic(InteractionType, const double&, const double&) const;
        bool IsCEX(InteractionType, const double&, const double&) const;
        bool IsSinglePion(InteractionType, const double&, const double&) const;

        std::shared_ptr<MetropolisData> GetProbe(const ProbeType&) const;

        // Data
        static bool registered;
        std::shared_ptr<MetropolisData> nucleon_data, pion_data;
};

}

#endif
