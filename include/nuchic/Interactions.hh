#ifndef INTERACTIONS_HH
#define INTERACTIONS_HH

#include <array>
#include <map>
#include <memory>
#include <vector>

#include "H5Cpp.h"

#include "nuchic/Interpolation.hh"
#include "nuchic/ThreeVector.hh"

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
        using TCreateMethod = std::unique_ptr<Interactions>(*)(const std::string&);

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
        static std::shared_ptr<Interactions> Create(const std::string&, const std::string&);

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

#endif // end of include guard: INTERACTIONS_HH
