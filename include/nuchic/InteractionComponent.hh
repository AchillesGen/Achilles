#ifndef INTERACTION_COMPONENT_HH
#define INTERACTION_COMPONENT_HH

#include <array>
#include <map>
#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

#include "H5Cpp.h"

#include "nuchic/InteractionComponentEnums.hh"
#include "nuchic/Interpolation.hh"
#include "nuchic/Random.hh"
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

/// Base class for implementing interaction models. These interaction models focus on the
/// interactions that occur between hadrons during the intranuclear cascade. This base class
/// is **not** to be used to inherit from for modeling the hard interactions between leptons
/// and the nucleus.
class InteractionComponent {
    public:
        /// @name Constructors and Destructors
        ///@{

        /// Base class constructor
        InteractionComponent() = default;
        InteractionComponent(const InteractionComponent&) = default;
        InteractionComponent(InteractionComponent&&) = default;
        InteractionComponent& operator=(const InteractionComponent&) = default;
        InteractionComponent& operator=(InteractionComponent&&) = default;

        /// Default destructor
        virtual ~InteractionComponent() = default;
        ///@}

        /// Function to determine the total cross-section between two particles
        ///@param part1: The first particle involved with the interaction
        ///@param part2: The second particle involved with the interaction
        ///@return double: The cross-section
        virtual double CrossSection(const Particle&, const Particle&) const = 0;

        /// Function to calculate all possible cross-sections between two particles
        /// as a vector for each unique final state.
        ///@param part1: The first particle involved in the interaction in CM frame
        ///@param part2: The second particle involved in the interaction in CM frame
        ///@return std::vector<double>: Vector of all allowed cross-sections
        virtual std::vector<double> CrossSections(const Particle&, const Particle&) const = 0;

        /// Function to determine the final state of the calculation given the cross-sections
        /// and a random number.
        ///@param part1: The first particle involved in the interaction in CM frame
        ///@param part2: The second particle involved in the interaction in CM frame
        ///@return std::vector<Particle>: Vector of generated final state particles
        virtual std::vector<Particle> GenerateFinalState(const Particle&, const Particle&) const = 0;

        /// Function that returns the interactions calculated in the class
        virtual InteractionComponentType InteractionType() const = 0;

    protected:
        double CrossSectionLab(bool, const double&) const noexcept;
};


/// The InteractionComponentFactory creates a method of generating an interaction component for the 
/// intranuclear cascade at runtime by a string in the run card. This allows the user to 
/// test how different interaction models may effect results without having to recompile any
/// code.
class InteractionComponentFactory {
    public:
        /// Typedef for easier coding
        using TCreateMethod = std::unique_ptr<InteractionComponent>(*)(const YAML::Node&);

        InteractionComponentFactory() = delete;

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
        static std::unique_ptr<InteractionComponent> Create(const YAML::Node&);

        /// Produce a list of all the registered interactions
        static void ListInteractions();

    private:
        static auto &methods() {
            static std::unordered_map<std::string, TCreateMethod> map;
            return map;
        }
};

#define REGISTER_INTERACTION_COMPONENT(interaction) \
    bool interaction::registered = InteractionComponentFactory::Register(interaction::GetName(), \
                                                                         interaction::Create)

class NasaInteractions : public InteractionComponent {
    public:
        NasaInteractions(const YAML::Node&) {};

        static std::unique_ptr<InteractionComponent> Create(const YAML::Node& data) {
            return std::make_unique<NasaInteractions>(data);
        }

        static std::string GetName() { return "NasaInteractions"; }

        // These functions are defined in the base class
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle&, const Particle&) const override;
        std::vector<double> CrossSections(const Particle &part1,
                                          const Particle &part2) const override {
            return { CrossSection(part1, part2) };
        }
        std::vector<Particle> GenerateFinalState(const Particle&, const Particle&) const override;

        /// Function that returns the interactions calculated in the class
        InteractionComponentType InteractionType() const override { return InteractionComponentType::NucleonNucleon; }
    private:
        static bool registered;
};

/// Class for implementing an interaction model based on the Geant4 cross-section data. This
/// interaction model contains information about the angular distribution of pp, pn, and nn
/// interactions that occur during the intranuclear cascade.
class GeantInteractions : public InteractionComponent {
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
        static std::unique_ptr<InteractionComponent> Create(const YAML::Node& data) {
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
        std::vector<double> CrossSections(const Particle &part1,
                                          const Particle &part2) const override {
            return { CrossSection(part1, part2) };
        }
        std::vector<Particle> GenerateFinalState(const Particle&, const Particle&) const override;

        /// Function that returns the interactions calculated in the class
        InteractionComponentType InteractionType() const override { return InteractionComponentType::NucleonNucleon; }
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
class ConstantInteractions : public InteractionComponent {
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
        static std::unique_ptr<InteractionComponent> Create(const YAML::Node& data) {
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
        std::vector<double> CrossSections(const Particle&, const Particle&) const override {
            return { m_xsec };
        }
        std::vector<Particle> GenerateFinalState(const Particle&, const Particle&) const override;

        /// Function that returns the interactions calculated in the class
        InteractionComponentType InteractionType() const override { return InteractionComponentType::NucleonNucleon; }

    private:
        // Variables
        double m_xsec;
        static bool registered;
};

/// Class for implementing an interaction model based on the pion-nucleus cross-section calculations. This
/// interaction model contains information about the angular distribution of Npi
/// interactions that occur during the intranuclear cascade.
class PionNucleon : public InteractionComponent {
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
        static std::unique_ptr<InteractionComponent> Create(const YAML::Node& data) {
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
        double CrossSection(const Particle&, const Particle&) const override;
        std::vector<double> CrossSections(const Particle&, const Particle&) const override;
        std::vector<Particle> GenerateFinalState(const Particle&, const Particle&) const override;

        /// Function that returns the interactions calculated in the class
        InteractionComponentType InteractionType() const override { return InteractionComponentType::NucleonPion; }

    private:
        // Structure to hold the cross-section data
        struct cross_section {
            std::vector<std::pair<PID, PID>> m_pids;
            std::vector<Interp1D> m_cross_sections;
            std::vector<Interp2D> m_theta_dist;
            std::vector<double> m_energies, m_angles;
            std::vector<double> m_sigma1, m_sigma2;
            std::vector<double> m_theta, m_cdf;

            void CalcCDF(const std::vector<double>&, const std::vector<double>&);
        };

        // Helper functions
        cross_section LoadData(const std::string&, const std::vector<std::pair<PID, PID>>&) const;
        size_t SelectFinalState(const std::pair<PID, PID>&, const double&, const double&) const;
        double GenerateAngle(const std::pair<PID, PID>&, const size_t&, const double&, const double&) const;
        
        // Variables
        static bool registered;
        static constexpr size_t nfiles = 6;
        std::array<std::string, nfiles> filenames{"/pi0-n.out", "/pi0-p.out",
                                                  "/pim-n.out", "/pim-p.out",
                                                  "/pip-n.out", "/pip-p.out"}; 

        std::map<std::pair<PID, PID>, cross_section> m_cross_sections;
};

}

#endif // end of include guard: INTERACTIONS_HH
