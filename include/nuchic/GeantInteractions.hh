#ifndef GEANTINTERACTIONS_HH
#define GEANTINTERACTIONS_HH

#include "nuchic/Interactions.hh"
#include "nuchic/Interpolation.hh"

#include "H5Cpp.h"

namespace nuchic {

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class GeantInteractions : public Interactions {
    public:
        ///@name Constructors and Destructors
        ///@{

        /// Initialize GeantInteractions class. This loads data from an input file
        ///@param filename: The location of the Geant4 hdf5 data file
        GeantInteractions(const std::string&);
        GeantInteractions(const GeantInteractions&) = default;
        GeantInteractions(GeantInteractions&&) = default;
        GeantInteractions& operator=(const GeantInteractions&) = default;
        GeantInteractions& operator=(GeantInteractions&&) = default;

        /// Generate a GeantInteractions object. This is used in the InteractionFactory.
        ///@param data: The location of the data file to load containing the Geant4 cross-sections
        static std::unique_ptr<Interactions> Create(const std::string& data) {
            return std::make_unique<GeantInteractions>(data);
        }

        /// Default Destructor
        ~GeantInteractions() override = default;
        ///@}

        /// Returns the name of the class, used in the InteractionFactory
        ///@return std::string: The name of the class
        static std::string GetName() { return "GeantInteractions"; }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2>&) const override;
    private:
        // Functions
        double CrossSectionAngle(bool, const double&, const double&) const;
        void LoadData(bool, const H5::Group&);

        // Variables
        std::vector<double> m_theta, m_cdf;
        std::vector<double> m_pcmPP, m_xsecPP;
        std::vector<double> m_pcmNP, m_xsecNP;
        Interp1D m_crossSectionPP, m_crossSectionNP;
        Interp2D m_thetaDistPP, m_thetaDistNP;
        static bool registered;
};

}

#endif
