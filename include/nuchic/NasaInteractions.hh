#ifndef NASAINTERACTIONS_HH
#define NASAINTERACTIONS_HH

#include "nuchic/Interactions.hh"
#include "nuchic/Interpolation.hh"

namespace nuchic {

class NasaInteractions : public Interactions {
    public:
        NasaInteractions(const std::string&) {};

        static std::unique_ptr<Interactions> Create(const std::string& data) {
            return std::make_unique<NasaInteractions>(data);
        }

        static std::string GetName() { return "NasaInteractions"; }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2>&) const override;
    private:
        static bool registered;
};

}

#endif
