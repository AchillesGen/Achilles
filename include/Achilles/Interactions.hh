#pragma once

#include <array>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Achilles/Factory.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/ThreeVector.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wuseless-cast"
#endif
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include "highfive/H5Group.hpp"
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

class Particle;
class FourVector;
class Potential;
class Random;

struct InteractionResult {
    std::vector<PID> particles;
    double cross_section;
};

using InteractionResults = std::vector<InteractionResult>;

/// Base class for implementing interaction models. These interaction models focus on the
/// interactions that occur between hadrons during the intranuclear cascade. This base class
/// is **not** to be used to inherit from for modeling the hard interactions between leptons
/// and the nucleus.
class Interaction {
  public:
    /// @name Constructors and Destructors
    ///@{

    /// Base class constructor
    Interaction() = default;
    Interaction(const Interaction &) = default;
    Interaction(Interaction &&) = default;
    Interaction &operator=(const Interaction &) = default;
    Interaction &operator=(Interaction &&) = default;

    /// Base class constructor for classes that need data files
    // Interactions(const std::string& data) {}

    /// Default destructor
    virtual ~Interaction() = default;
    ///@}

    /// Function to list all implemented interactions between two particles
    ///@return std::vector<std::pair<PID, PID>>: List of all interactions
    virtual std::vector<std::pair<PID, PID>> InitialStates() const = 0;

    /// Function that returns the total cross-section between two particles
    ///@param part1: The first particle involved with the interaction
    ///@param part2: The second particle involved with the interaction
    ///@return double: The cross-section
    virtual double TotalCrossSection(const Particle &, const Particle &) const;

    /// Function to determine the cross-section between two particles
    /// broken down into all possible final state channels
    ///@param part1: The first particle involved with the interaction
    ///@param part2: The second particle involved with the interaction
    ///@return double: The cross-section
    virtual InteractionResults CrossSection(const Particle &, const Particle &) const = 0;

    /// Function to generate the final state particles
    ///@param part1: The first particle involved with the interaction
    ///@param part2: The second particle involved with the interaction
    ///@param out_pids: The final state particles to generate
    ///@param rands: The random numbers to use in the generation
    ///@return std::vector<Particle>: The final state particles
    virtual std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                                   const std::vector<PID> &out_pids,
                                                   Random &) const = 0;

    virtual std::string GetName() const = 0;
    static std::string Name() { return "Interaction"; }

  protected:
    double CrossSectionLab(bool, const double &) const noexcept;
};

template <typename Derived>
using RegistrableInteraction = Registrable<Interaction, Derived, const YAML::Node &>;
using InteractionFactory = Factory<Interaction, const YAML::Node &>;

class NasaInteraction : public Interaction, RegistrableInteraction<NasaInteraction> {
  public:
    NasaInteraction() = default;
    NasaInteraction(const YAML::Node &){};

    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<NasaInteraction>(data);
    }

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return NasaInteraction::Name(); }
    static std::string Name() { return "NasaInteraction"; }

    // These functions are defined in the base class
    std::vector<std::pair<PID, PID>> InitialStates() const override {
        return {{PID::proton(), PID::proton()},
                {PID::neutron(), PID::proton()},
                {PID::neutron(), PID::neutron()}};
    }
    InteractionResults CrossSection(const Particle &, const Particle &) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

  private:
    ThreeVector MakeMomentum(bool, double, const std::vector<double> &) const;
};

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
    std::vector<std::pair<PID, PID>> InitialStates() const override {
        return {{PID::proton(), PID::proton()},
                {PID::neutron(), PID::proton()},
                {PID::neutron(), PID::neutron()}};
    }
    InteractionResults CrossSection(const Particle &, const Particle &) const override;
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

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class ConstantInteraction : public Interaction, RegistrableInteraction<ConstantInteraction> {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize ConstantInteractions class for testing purposes.
    ConstantInteraction() = default;
    /// Initialize ConstantInteractions class. This loads data from an input file
    ///@param node: The constant interactions to load
    ConstantInteraction(const YAML::Node &);
    ConstantInteraction(const ConstantInteraction &) = default;
    ConstantInteraction(ConstantInteraction &&) = default;
    ConstantInteraction &operator=(const ConstantInteraction &) = default;
    ConstantInteraction &operator=(ConstantInteraction &&) = default;

    /// Generate a ConstantInteractions object. This is used in the InteractionFactory.
    ///@param data: The location of the data file to load containing the Geant4 cross-sections
    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<ConstantInteraction>(data);
    }

    /// Default Destructor
    ~ConstantInteraction() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    std::string GetName() const override { return ConstantInteraction::Name(); }
    static std::string Name() { return "ConstantInteraction"; }

    // These functions are defined in the base class
    std::vector<std::pair<PID, PID>> InitialStates() const override { return m_states; }
    InteractionResults CrossSection(const Particle &, const Particle &) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           Random &) const override;

    // These functions are for testing only
    void AddInteraction(const std::pair<PID, PID> &state, const InteractionResults &results) {
        if(m_interactions.find(state) != m_interactions.end()) {
            auto msg =
                fmt::format("Initial state: [{}, {}] already exists", state.first, state.second);
            throw std::runtime_error("Initial state already exists");
        }
        m_interactions[state] = results;
        m_states.push_back(state);
    }

  private:
    ThreeVector MakeMomentum(bool, double, const std::vector<double> &) const;

    // Variables
    std::map<std::pair<PID, PID>, InteractionResults> m_interactions;
    std::vector<std::pair<PID, PID>> m_states;
};

} // namespace achilles

namespace YAML {
template <> struct convert<achilles::InteractionResult> {
    static bool decode(const Node &node, achilles::InteractionResult &result) {
        result.particles = node["Outgoing"].as<std::vector<achilles::PID>>();
        result.cross_section = node["CrossSection"].as<double>();
        return true;
    }
};

template <> struct convert<std::vector<achilles::InteractionResult>> {
    static bool decode(const Node &node, std::vector<achilles::InteractionResult> &results) {
        for(auto result : node) results.push_back(result.as<achilles::InteractionResult>());
        return true;
    }
};

} // namespace YAML
