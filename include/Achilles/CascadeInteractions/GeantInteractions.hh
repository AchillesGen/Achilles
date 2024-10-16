#ifndef ACHILLES_CASCADEINTERACTIONS_GEANTINTERACTIONS
#define ACHILLES_CASCADEINTERACTIONS_GEANTINTERACTIONS

#include "Achilles/Interactions.hh"
#include "Achilles/Interpolation.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wuseless-cast"
#endif
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include "highfive/H5Group.hpp"
#pragma GCC diagnostic pop

namespace achilles {

class ThreeVector;

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class GeantInteraction : public Interaction, RegistrableInteraction<GeantInteraction> {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize GeantInteractions class. This loads data from an input file
    ///@param filename: The location of the Geant4 hdf5 data file
    GeantInteraction(const YAML::Node &);
    GeantInteraction(const GeantInteraction &) = default;
    GeantInteraction(GeantInteraction &&) = default;
    GeantInteraction &operator=(const GeantInteraction &) = default;
    GeantInteraction &operator=(GeantInteraction &&) = default;

    /// Generate a GeantInteractions object. This is used in the InteractionFactory.
    ///@param data: The location of the data file to load containing the Geant4 cross-sections
    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<GeantInteraction>(data);
    }

    /// Default Destructor
    ~GeantInteraction() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return GeantInteraction::Name(); }
    static std::string Name() { return "GeantInteraction"; }

    // These functions are defined in the base class
    std::vector<std::pair<PID, PID>> InitialStates() const override;
    InteractionResults CrossSection(Event &, size_t, size_t) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

  private:
    // Functions
    double CrossSectionAngle(bool, const double &, const double &) const;
    void LoadData(bool, const HighFive::Group &);
    ThreeVector MakeMomentum(bool, double, const std::vector<double> &) const;

    // Variables
    std::vector<double> m_theta, m_cdf;
    std::vector<double> m_pcmPP, m_xsecPP;
    std::vector<double> m_pcmNP, m_xsecNP;
    Interp1D m_crossSectionPP, m_crossSectionNP;
    Interp2D m_thetaDistPP, m_thetaDistNP;
};

} // namespace achilles

#endif
