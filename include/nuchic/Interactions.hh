#ifndef INTERACTIONS_HH
#define INTERACTIONS_HH

#include <array>
#include <unordered_map>
#include <memory>
#include <vector>

#include "nuchic/ThreeVector.hh"

namespace nuchic {

class Particle;

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
        virtual ThreeVector MakeMomentum(bool, const double&, const double&,
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

}

#endif // end of include guard: INTERACTIONS_HH
