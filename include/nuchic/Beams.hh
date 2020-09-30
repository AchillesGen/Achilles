#ifndef BEAMS_HH
#define BEAMS_HH

#include <memory>

namespace nuchic {

class FluxType {
    public:
        FluxType() = default;
        FluxType(const FluxType&) = default;
        FluxType(FluxType&&) = default;
        FluxType& operator=(const FluxType&) = default;
        FluxType& operator=(FluxType&&) = default;
        virtual ~FluxType() = default;

        virtual double Flux(const double&) = 0;
        virtual double Weight(const double&) = 0;
};

class Monochromatic : public FluxType {
    public:
        Monochromatic(const double &energy) : m_energy(energy) {}

        double Flux(const double&) override {
            return m_energy;
        }

        double Weight(const double&) override {
            return 1;
        }

    private:
        double m_energy;
};


// TODO: Figure out how to handle spectrums
// How to handle this correctly?
// 1. Generate the energy according to some distribution
// 2. A. Multiply maximum energy by a fraction to get neutrino energy
//    B. Calculate the probability for this energy
// The first requires being able to appropriately sample from the distribution
// (i.e. having a form for the inverse of the CDF). The second would take advantage
// of the importance sampling of Vegas.
class Spectrum : public FluxType {
    public:
        Spectrum(const std::string&) {
            throw std::runtime_error("Spectrum Fluxes are not implemented");
        }

        double Flux(const double&) override {
            throw std::runtime_error("Spectrum Fluxes are not implemented");
        }

        double Weight(const double&) override {
            throw std::runtime_error("Spectrum Fluxes are not implemented");
        }
};

class Beam {
    public:
        Beam(std::unique_ptr<FluxType> &&flux) : m_flux(std::move(flux)) {}

    private:
        std::unique_ptr<FluxType> m_flux;
};


}

#endif
