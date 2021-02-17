#ifndef INTERACTIONS_HH
#define INTERACTIONS_HH

#include <array>
#include <map>
#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

#include "H5Cpp.h"

#include "nuchic/Interpolation.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/ParticleInfo.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace nuchic {

class Particle;

double CrossSection(bool, const double&);
double CrossSectionAngle(bool, const double&, const double&);
double CrossSectionLab(bool, const double&) noexcept;
ThreeVector MakeMomentumAngular(bool, const double&, const double&, const std::array<double, 2>&);

enum class FSInteractionType {
    NucleonNucleon,
    NucleonPion,
    NucleonDelta,
    DeltaProduction
};

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

        /// Function that returns the interactions calculated in the class
        virtual FSInteractionType InteractionType() const = 0;

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
        static void ListInteractions();

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

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2>&) const override;

        /// Function that returns the interactions calculated in the class
        FSInteractionType InteractionType() const override { return FSInteractionType::NucleonNucleon; }
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

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2>&) const override;

        /// Function that returns the interactions calculated in the class
        FSInteractionType InteractionType() const override { return FSInteractionType::NucleonNucleon; }
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

        /// Function that returns the interactions calculated in the class
        FSInteractionType InteractionType() const override { return FSInteractionType::NucleonNucleon; }

    private:
        // Variables
        double m_xsec;
        static bool registered;
};

/// Class for implementing an interaction model based on the pion-nucleus cross-section calculations. This
/// interaction model contains information about the angular distribution of Npi
/// interactions that occur during the intranuclear cascade.
class PionNucleon : public Interactions {
    public:
        ///@name Constructors and Destructors
        ///@{

        /// Initialize PionNucleon class. This loads data from an input file
        ///@param filename: The location of the Geant4 hdf5 data file
        PionNucleon(const YAML::Node&);
        PionNucleon(const PionNucleon&) = default;
        PionNucleon(PionNucleon&&) = default;
        PionNucleon& operator=(const PionNucleon&) = default;
        PionNucleon& operator=(PionNucleon&&) = default;

        /// Generate a PionNucleon object. This is used in the InteractionFactory.
        ///@param data: The location of the data file to load containing the Geant4 cross-sections
        static std::unique_ptr<Interactions> Create(const YAML::Node& data) {
            return std::make_unique<PionNucleon>(data);
        }

        /// Default Destructor
        ~PionNucleon() override = default;
        ///@}

        /// Returns the name of the class, used in the InteractionFactory
        ///@return std::string: The name of the class
        static std::string GetName() { return "PionNucleon"; }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override { return 0; }
        ThreeVector MakeMomentum(bool, const double&,
                                 const std::array<double, 2>&) const override { return {0, 0, 0}; }

        /// Function that returns the interactions calculated in the class
        FSInteractionType InteractionType() const override { return FSInteractionType::NucleonPion; }

    private:
        // Structure to hold the cross-section data
        struct cross_section {
            std::vector<std::pair<PID, PID>> m_pids;
            std::vector<Interp2D> m_cross_sections;
        };

        // Helper function to load data
        cross_section LoadData(const std::string&, const std::vector<std::pair<PID, PID>>&) const;
        
        // Variables
        double m_xsec;
        static bool registered;
        std::array<std::string, 6> filenames{"pi0-n.out", "pi0-p.out",
                                             "pim-n.out", "pim-p.out",
                                             "pip-n.out", "pip-p.out"}; 

        cross_section pi0_n, pi0_p, pim_n, pim_p, pip_n, pip_p;
};

}

#endif // end of include guard: INTERACTIONS_HH
