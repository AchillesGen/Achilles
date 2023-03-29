#ifndef INTERACTIONS_HH
#define INTERACTIONS_HH

#include <array>
#include <map>
#include <unordered_map>
#include <memory>
#include <vector>

#include "Achilles/ThreeVector.hh"
#include "Achilles/Interpolation.hh"

// #include "H5Cpp.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wuseless-cast"
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

/// Base class for implementing interaction models. These interaction models focus on the
/// interactions that occur between hadrons during the intranuclear cascade. This base class
/// is **not** to be used to inherit from for modeling the hard interactions between leptons
/// and the nucleus.
class Interactions {
    public:
        /// @name Constructors and Destructors
        ///@{

        /// Base class constructor
        Interactions() = default;
        Interactions(const Interactions&) = default;
        Interactions(Interactions&&) = default;
        Interactions& operator=(const Interactions&) = default;
        Interactions& operator=(Interactions&&) = default;

        /// Base class constructor for classes that need data files
        // Interactions(const std::string& data) {}

        /// Default destructor
        virtual ~Interactions() = default;
        ///@}

        /// Function to determine the cross-section between two particles
        ///@param part1: The first particle involved with the interaction
        ///@param part2: The second particle involved with the interaction
        ///@return double: The cross-section
        virtual double CrossSection(const Particle&, const Particle&) const = 0;

        /// Function to generate momentum for the particles after an interaction
        ///@param samePID: Used to determine if the two particles are the same type
        ///@param p1CM: The momentum of the first particle in the center of mass frame
        ///@param pcm: The center of mass momentum
        ///@param rans: An array containing two random numbers used to generate phase space
        virtual ThreeVector MakeMomentum(bool, const double&,
                                         const std::array<double, 2>&) const = 0;

        using MomentumPair = std::pair<FourVector, FourVector>;
        virtual MomentumPair FinalizeMomentum(const Particle&,
                                              const Particle&,
                                              std::shared_ptr<Potential>) const;
        virtual std::string Name() const = 0;
    protected:
        double CrossSectionLab(bool, const double&) const noexcept;
};


/// The InteractionFactory creates a method of generating an interaction model for the 
/// intranuclear cascade at runtime by a string in the run card. This allows the user to 
/// test how different interaction models may effect results without having to recompile any
/// code.
class InteractionFactory {
    public:
        /// Typedef for easier coding
        using TCreateMethod = std::unique_ptr<Interactions>(*)(const YAML::Node&);

        InteractionFactory() = delete;

        /// Register a new Interactions subclass to the factory
        ///@param name: The name of the subclass
        ///@param funcCreate: The Create function of the subclass
        ///@return bool: True if the class was registered to the factory, False if the
        ///              registration failed
        static bool Register(const std::string&, TCreateMethod);

        /// Create an instance of the desired subclass based on the input string
        ///@param name: Name of the subclass object to be created
        ///@return std::shared_ptr<Interactions>: A pointer to the interaction subclass object
        ///     (**NOTE**: Returns nullptr if name is not registered)
        static std::unique_ptr<Interactions> Create(const YAML::Node&);

        /// Produce a list of all the registered interactions
        static void Display();

    private:
        static auto &methods() {
            static std::unordered_map<std::string, TCreateMethod> map;
            return map;
        }
};

#define REGISTER_INTERACTION(interaction) \
    bool interaction::registered = InteractionFactory::Register(interaction::GetName(), \
                                                                interaction::Create)
class NasaInteractions : public Interactions {
    public:
        NasaInteractions(const YAML::Node&) {};

        static std::unique_ptr<Interactions> Create(const YAML::Node& data) {
            return std::make_unique<NasaInteractions>(data);
        }

        static std::string GetName() { return "NasaInteractions"; }
        std::string Name() const override { return NasaInteractions::GetName(); }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2>&) const override;
    private:
        static bool registered;
};

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class GeantInteractions : public Interactions {
    public:
        ///@name Constructors and Destructors
        ///@{

        /// Initialize GeantInteractions class. This loads data from an input file
        ///@param filename: The location of the Geant4 hdf5 data file
        GeantInteractions(const YAML::Node&);
        GeantInteractions(const GeantInteractions&) = default;
        GeantInteractions(GeantInteractions&&) = default;
        GeantInteractions& operator=(const GeantInteractions&) = default;
        GeantInteractions& operator=(GeantInteractions&&) = default;

        /// Generate a GeantInteractions object. This is used in the InteractionFactory.
        ///@param data: The location of the data file to load containing the Geant4 cross-sections
        static std::unique_ptr<Interactions> Create(const YAML::Node& data) {
            return std::make_unique<GeantInteractions>(data);
        }

        /// Default Destructor
        ~GeantInteractions() override = default;
        ///@}

        /// Returns the name of the class, used in the InteractionFactory
        ///@return std::string: The name of the class
        static std::string GetName() { return "GeantInteractions"; }
        std::string Name() const override { return GeantInteractions::GetName(); }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2>&) const override;
    private:
        // Functions
        double CrossSectionAngle(bool, const double&, const double&) const;
        void LoadData(bool, const HighFive::Group&);

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
class ConstantInteractions : public Interactions {
    public:
        ///@name Constructors and Destructors
        ///@{

        /// Initialize ConstantInteractions class. This loads data from an input file
        ///@param filename: The location of the Geant4 hdf5 data file
        ConstantInteractions(const YAML::Node& xsec) : m_xsec(xsec["CrossSection"].as<double>()) {}
        ConstantInteractions(const ConstantInteractions&) = default;
        ConstantInteractions(ConstantInteractions&&) = default;
        ConstantInteractions& operator=(const ConstantInteractions&) = default;
        ConstantInteractions& operator=(ConstantInteractions&&) = default;

        /// Generate a ConstantInteractions object. This is used in the InteractionFactory.
        ///@param data: The location of the data file to load containing the Geant4 cross-sections
        static std::unique_ptr<Interactions> Create(const YAML::Node& data) {
            return std::make_unique<ConstantInteractions>(data);
        }

        /// Default Destructor
        ~ConstantInteractions() override = default;
        ///@}

        /// Returns the name of the class, used in the InteractionFactory
        ///@return std::string: The name of the class
        static std::string GetName() { return "ConstantInteractions"; }
        std::string Name() const override { return ConstantInteractions::GetName(); }

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

#endif // end of include guard: INTERACTIONS_HH
