#pragma once

#include <array>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Achilles/Interpolation.hh"
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
class PID;
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
    virtual std::vector<std::pair<int, int>> InitialStates() const = 0;

    /// Function to determine the cross-section between two particles
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
                                                   const std::vector<double> &rands) const = 0;

    virtual std::string Name() const = 0;

  protected:
    double CrossSectionLab(bool, const double &) const noexcept;
};

/// The InteractionFactory creates a method of generating an interaction model for the
/// intranuclear cascade at runtime by a string in the run card. This allows the user to
/// test how different interaction models may effect results without having to recompile any
/// code.
class InteractionFactory {
  public:
    /// Typedef for easier coding
    using TCreateMethod = std::unique_ptr<Interaction> (*)(const YAML::Node &);

    InteractionFactory() = delete;

    /// Register a new Interactions subclass to the factory
    ///@param name: The name of the subclass
    ///@param funcCreate: The Create function of the subclass
    ///@return bool: True if the class was registered to the factory, False if the
    ///              registration failed
    static bool Register(const std::string &, TCreateMethod);

    /// Create an instance of the desired subclass based on the input string
    ///@param name: Name of the subclass object to be created
    ///@return std::shared_ptr<Interactions>: A pointer to the interaction subclass object
    ///     (**NOTE**: Returns nullptr if name is not registered)
    static std::unique_ptr<Interaction> Create(const YAML::Node &);

    /// Produce a list of all the registered interactions
    static void Display();

  private:
    static auto &methods() {
        static std::unordered_map<std::string, TCreateMethod> map;
        return map;
    }
};

#define REGISTER_INTERACTION(interaction) \
    bool interaction::registered =        \
        InteractionFactory::Register(interaction::GetName(), interaction::Create)

class NasaInteraction : public Interaction {
  public:
    NasaInteraction() = default;
    NasaInteraction(const YAML::Node &){};

    static std::unique_ptr<Interaction> Create(const YAML::Node &data) {
        return std::make_unique<NasaInteraction>(data);
    }

    static std::string GetName() { return "NasaInteractions"; }
    std::string Name() const override { return NasaInteraction::GetName(); }

    // These functions are defined in the base class
    static bool IsRegistered() noexcept { return registered; }
    std::vector<std::pair<int, int>> InitialStates() const override {
        return {{2112, 2112}, {2112, 2212}, {2212, 2212}};
    }
    InteractionResults CrossSection(const Particle &, const Particle &) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           const std::vector<double> &rands) const override;

  private:
    ThreeVector MakeMomentum(bool, double, const std::vector<double> &) const;
    static bool registered;
};

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class GeantInteraction : public Interaction {
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
    static std::unique_ptr<Interaction> Create(const YAML::Node &data) {
        return std::make_unique<GeantInteraction>(data);
    }

    /// Default Destructor
    ~GeantInteraction() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    static std::string GetName() { return "GeantInteractions"; }
    std::string Name() const override { return GeantInteraction::GetName(); }

    // These functions are defined in the base class
    static bool IsRegistered() noexcept { return registered; }
    std::vector<std::pair<int, int>> InitialStates() const override {
        return {{2112, 2112}, {2112, 2212}, {2212, 2212}};
    }
    InteractionResults CrossSection(const Particle &, const Particle &) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           const std::vector<double> &rands) const override;

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
    static bool registered;
};

/*
/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class GeantInteractionsDt : public Interactions {
    public:
        ///@name Constructors and Destructors
        ///@{

        /// Initialize GeantInteractionsDt class. This loads data from an input file
        ///@param filename: The location of the Geant4 hdf5 data file
        GeantInteractionsDt(const YAML::Node&);
        GeantInteractionsDt(const GeantInteractionsDt&) = default;
        GeantInteractionsDt(GeantInteractionsDt&&) = default;
        GeantInteractionsDt& operator=(const GeantInteractionsDt&) = default;
        GeantInteractionsDt& operator=(GeantInteractionsDt&&) = default;

        /// Generate a GeantInteractionsDt object. This is used in the InteractionFactory.
        ///@param data: The location of the data file to load containing the Geant4 cross-sections
        static std::unique_ptr<Interactions> Create(const YAML::Node& data) {
            return std::make_unique<GeantInteractionsDt>(data);
        }

        /// Default Destructor
        ~GeantInteractionsDt() override = default;
        ///@}

        /// Returns the name of the class, used in the InteractionFactory
        ///@return std::string: The name of the class
        static std::string GetName() { return "GeantInteractionsDt"; }
        std::string Name() const override { return GeantInteractionsDt::GetName(); }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2>&) const override { return {}; }
        MomentumPair FinalizeMomentum(const Particle&,
                                      const Particle&,
                                      std::shared_ptr<Potential>) const override;
    private:
        // Functions
        double CrossSectionAngle(bool, const double&, const double&) const { return 0; }
        void LoadData(bool, const std::string&);
        std::vector<double> ReadBlock(std::ifstream &file, size_t nlines=1) const;

        // Variables
        Interp1D m_crossSectionPP, m_crossSectionNP;
        Interp2D m_pdfPP, m_pdfNP;
        static bool registered;
};
*/

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class ConstantInteraction : public Interaction {
  public:
    ///@name Constructors and Destructors
    ///@{

    /// Initialize ConstantInteractions class. This loads data from an input file
    ///@param filename: The location of the Geant4 hdf5 data file
    ConstantInteraction(const YAML::Node &xsec) : m_xsec(xsec["CrossSection"].as<double>()) {}
    ConstantInteraction(double xsec) : m_xsec(xsec) {}
    ConstantInteraction(const ConstantInteraction &) = default;
    ConstantInteraction(ConstantInteraction &&) = default;
    ConstantInteraction &operator=(const ConstantInteraction &) = default;
    ConstantInteraction &operator=(ConstantInteraction &&) = default;

    /// Generate a ConstantInteractions object. This is used in the InteractionFactory.
    ///@param data: The location of the data file to load containing the Geant4 cross-sections
    static std::unique_ptr<Interaction> Create(const YAML::Node &data) {
        return std::make_unique<ConstantInteraction>(data);
    }

    /// Default Destructor
    ~ConstantInteraction() override = default;
    ///@}

    /// Returns the name of the class, used in the InteractionFactory
    ///@return std::string: The name of the class
    static std::string GetName() { return "ConstantInteractions"; }
    std::string Name() const override { return ConstantInteraction::GetName(); }

    // These functions are defined in the base class
    static bool IsRegistered() noexcept { return registered; }
    std::vector<std::pair<int, int>> InitialStates() const override {
        return {{2112, 2112}, {2112, 2212}, {2212, 2212}};
    }
    InteractionResults CrossSection(const Particle &, const Particle &) const override;
    std::vector<Particle> GenerateMomentum(const Particle &part1, const Particle &part2,
                                           const std::vector<PID> &out_pids,
                                           const std::vector<double> &rands) const override;

  private:
    ThreeVector MakeMomentum(bool, double, const std::vector<double> &) const;

    // Variables
    double m_xsec;
    static bool registered;
};

} // namespace achilles
