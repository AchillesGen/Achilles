#ifndef CASCADE_HH
#define CASCADE_HH

#include <array>
#include <memory>
#include <vector>

#include "Achilles/SymplecticIntegrator.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Random.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Interactions.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

class Nucleus;
class Particle;
class Event;

using Particles = std::vector<Particle>;
using InteractionDistances = std::vector<std::pair<std::size_t, double>>;

/// The Cascade class performs a cascade of the nucleons inside the nucleus. The nucleons that
/// are struck in the hard interaction propagate through the nuclear medium. To determine if an
/// interaction occurs, we calculate the interaction cross-section of Np and Nn, where N is the
/// propagating nucleon.
class Cascade {
    static constexpr int cMaxSteps = 100000;
    public:
        // Probability Enums
        enum ProbabilityType {
            Gaussian,
            Pion,
            Cylinder
        };

        // In-Medium Effects Enums
        enum InMedium {
            None,
            NonRelativistic,
            Relativistic
        };

        /// @name Constructor and Destructor
        ///@{

        /// Create the Cascade object
        ///@param interactions: The interaction model for pp, pn, and np interactions
        ///@param prob: The interaction probability function to be used
        ///@param dist: The maximum distance step to take when propagating
        ///TODO: Should the ProbabilityType be part of the interaction class or the cascade class?
        Cascade() = default;
        Cascade(std::unique_ptr<Interactions>,  const ProbabilityType&,
                const InMedium&, bool potential_prob=false, const double& dist=0.03);
        Cascade(Cascade&&) = default;
        Cascade& operator=(Cascade&&) = default;

        /// Default destructor
        ~Cascade() = default;
        ///@}
        
        /// @name Getters
        ///@{
        
        /// Get the name of the interaction model used
        ///@return std::string: Name of the interaction model
        std::string InteractionModel() const { return m_interactions -> Name(); } 

        /// Get the probability model used
        ///@return std::string: Name of the probability model
        std::string ProbabilityModel() const { return m_probability_name; }

        /// Get InMedium setting used
        ///@return std::string: InMedium setting
        std::string InMediumSetting() const {
            std::string name;
            switch(m_medium) {
                case None:
                    name = "None";
                    break;
                case NonRelativistic:
                    name = "NonRelativistic";
                    break;
                case Relativistic:
                    name = "Relativistic";
                    break;
            }
            return name;
        }

        /// Get potential prop option
        ///@return bool: PotentialProp option
        bool UsePotentialProp() const { return m_potential_prop; }

        /// Get step size
        ///@return double: default step size
        double StepSize() const { return distance; }

        /// @name Functions
        ///@{

        /// Give a random particle in the nucleus a kick. The kick is defined by a input
        /// four vector. To determine whether a proton or a neutron is kicked, the cross-section
        /// for protons and neutrons needs to be supplied.
        ///@param particles: The list of particles inside the nucleus
        ///@param energyTransfer: The energy transfered to the nucleon during the kick
        ///@param sigma: An array representing the cross-section for different kicked nucleons
        void Kick(std::shared_ptr<Nucleus>, const FourVector&, const std::array<double, 2>&);

        /// Reset the cascade internal variables for the next cascade
        void Reset();

        /// Helper function to make a specific nucleon as the kicked nucleon
        ///@param idx: The index of the particle that has been kicked
        void SetKicked(const std::size_t& idx) { kickedIdxs.push_back(idx); }

        /// Simulate the cascade until all particles either escape, are recaptured, or are in
        /// the background.
        ///@param nucleus: The nucleus to evolve
        ///@param maxSteps: The maximum steps to take in the cascade
        void Evolve(std::shared_ptr<Nucleus>, const std::size_t& maxSteps = cMaxSteps);

        /// Simulate the cascade on an event until all particles either escape,
        /// are recaptured, or are in the background.
        ///@param event: The event to run the cascade evolution on
        ///@param maxSteps: The maximum steps to take in the cascade
        void Evolve(achilles::Event*, const std::size_t& maxSteps = cMaxSteps);

        /// Simulate evolution of a kicked particle until it interacts for the
        /// first time with another particle, accumulating the total distance
        /// traveled by the kicked particle before it interacts.
        ///@param nucleus: The nucleus to evolve according to the mean free path calculation
        ///@param maxSteps: The maximum steps to take in the particle evolution
        void MeanFreePath(std::shared_ptr<Nucleus>, const std::size_t& maxSteps = cMaxSteps);

        /// Simulate the cascade until all particles either escape, are recaptured, or are in
        /// the background. This is done according to the NuWro algorithm.
        ///@param nucleus: The nucleus to evolve according to the NuWro method of cascade
        ///@param maxSteps: The maximum steps to take in the particle evolution
        void NuWro(std::shared_ptr<Nucleus>, const std::size_t& maxSteps = cMaxSteps);

        /// Simulate evolution of a kicked particle until it interacts for the
        /// first time with another particle, accumulating the total distance
        /// traveled by the kicked particle before it interacts.
        ///@param nucleus: The nucleus to evolve according to the mean free path calculation
        ///@param maxSteps: The maximum steps to take in the particle evolution
        void MeanFreePath_NuWro(std::shared_ptr<Nucleus>, const std::size_t& maxSteps = cMaxSteps);
        ///@}
    private:
        // Functions
        std::size_t GetInter(Particles&, const Particle&, double& stepDistance);
        void AdaptiveStep(const Particles&, const double&) noexcept;
        bool BetweenPlanes(const ThreeVector&, const ThreeVector&, const ThreeVector&, double) const noexcept;
        const ThreeVector Project(const ThreeVector&, const ThreeVector&, const ThreeVector&) const noexcept;
        const InteractionDistances AllowedInteractions(Particles&, const std::size_t&) const noexcept;
        double GetXSec(const Particle&, const Particle&) const;
        std::size_t Interacted(const Particles&, const Particle&,
                const InteractionDistances&) noexcept;
        void Escaped(Particles&);
        bool FinalizeMomentum(Particle&, Particle&) noexcept;
        bool PauliBlocking(const Particle&) const noexcept;
        void AddIntegrator(size_t, const Particle&);
        void Propagate(size_t, Particle*, double);
        void UpdateIntegrator(size_t, Particle*);

        // Variables
        std::vector<std::size_t> kickedIdxs;
        double distance{}, timeStep{};
        std::unique_ptr<Interactions> m_interactions;
        std::function<double(double, double)> probability;
        std::shared_ptr<Nucleus> localNucleus;
        InMedium m_medium;
        bool m_potential_prop;
        std::map<size_t, SymplecticIntegrator> integrators;
        std::string m_probability_name;
};

}

namespace YAML {
template<>
struct convert<achilles::Cascade> {
    static bool decode(const Node &node, achilles::Cascade &cascade) {
        auto interaction = achilles::InteractionFactory::Create(node["Interaction"]);
        auto probType = node["Probability"].as<achilles::Cascade::ProbabilityType>();
        auto mediumType = node["InMedium"].as<achilles::Cascade::InMedium>();
        auto potentialProp = node["PotentialProp"].as<bool>();
        auto distance = node["Step"].as<double>();
        cascade = achilles::Cascade(std::move(interaction), probType, mediumType, potentialProp, distance);
        return true;
    }
};

template<>
struct convert<achilles::Cascade::ProbabilityType> {
    static bool decode(const Node &node, achilles::Cascade::ProbabilityType &type) {
        if(node.as<std::string>() == "Gaussian")
            type = achilles::Cascade::ProbabilityType::Gaussian;
        else if(node.as<std::string>() == "Pion")
            type = achilles::Cascade::ProbabilityType::Pion;
        else if(node.as<std::string>() == "Cylinder")
            type = achilles::Cascade::ProbabilityType::Cylinder;
        else
            return false;
        return true;
    }
};

template<>
struct convert<achilles::Cascade::InMedium> {
    static bool decode(const Node &node, achilles::Cascade::InMedium &type) {
        if(node.as<std::string>() == "None")
            type = achilles::Cascade::InMedium::None;
        else if(node.as<std::string>() == "NonRelativistic")
            type = achilles::Cascade::InMedium::NonRelativistic;
        else if(node.as<std::string>() == "Relativistic")
            type = achilles::Cascade::InMedium::Relativistic;
        else
            return false;
        return true;
    }
};
}

#endif // end of include guard: CASCADE_HH
