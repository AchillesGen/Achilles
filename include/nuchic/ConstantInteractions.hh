#ifndef CONSTANTINTERACTIONS_HH
#define CONSTANTINTERACTIONS_HH

#include "nuchic/Interactions.hh"

namespace nuchic {

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class ConstantInteractions : public Interactions {
    public:
        ///@name Constructors and Destructors
        ///@{

        /// Initialize ConstantInteractions class. This loads data from an input file
        ///@param filename: The location of the Geant4 hdf5 data file
        ConstantInteractions(const std::string& xsec) : m_xsec(std::stod(xsec)) {}
        ConstantInteractions(const ConstantInteractions&) = default;
        ConstantInteractions(ConstantInteractions&&) = default;
        ConstantInteractions& operator=(const ConstantInteractions&) = default;
        ConstantInteractions& operator=(ConstantInteractions&&) = default;

        /// Generate a ConstantInteractions object. This is used in the InteractionFactory.
        ///@param data: The location of the data file to load containing the Geant4 cross-sections
        static std::unique_ptr<Interactions> Create(const std::string& data) {
            return std::make_unique<ConstantInteractions>(data);
        }

        /// Default Destructor
        ~ConstantInteractions() override = default;
        ///@}

        /// Returns the name of the class, used in the InteractionFactory
        ///@return std::string: The name of the class
        static std::string GetName() { return "ConstantInteractions"; }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override { return m_xsec; }
        ThreeVector MakeMomentum(bool, const double& pcm,
                                 const std::array<double, 2>& rans) const override {
            double ctheta = 2*rans[0]-1;
            double stheta = sqrt(1-ctheta*ctheta);
            double phi = 2*M_PI*rans[1];
            return pcm*ThreeVector(stheta*cos(phi), stheta*sin(phi), ctheta);
        }

    private:
        // Variables
        double m_xsec;
        static bool registered;
};

}

#endif
