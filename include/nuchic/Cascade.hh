#ifndef CASCADE_HH
#define CASCADE_HH

#include <array>
#include <memory>
#include <vector>

#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Random.hh"
#include "nuchic/Interpolation.hh"


namespace nuchic {

class Nucleus;
class Particle;
class Interactions;

using Particles = std::vector<Particle>;
using InteractionDistances = std::vector<std::pair<std::size_t, double>>;

/// The Cascade class performs a cascade of the nucleons inside the nucleus. The nucleons that
/// are struck in the hard interaction propagate through the nuclear medium. To determine if an
/// interaction occurs, we calculate the interaction cross-section of Np and Nn, where N is the
/// propagating nucleon.
class Cascade {
    public:
        // Probability Enums
        enum ProbabilityType {
            Gaussian,
            Pion
        };

        /// @name Constructor and Destructor
        ///@{

        /// Create the Cascade object
        ///@param interactions: The interaction model for pp, pn, and np interactions
        ///@param prob: The interaction probability function to be used
        ///@param dist: The maximum distance step to take when propagating
        ///TODO: Should the ProbabilityType be part of the interaction class or the cascade class?
        Cascade(const std::shared_ptr<Interactions>,  const ProbabilityType&, const double&);
        Cascade(const Cascade&) = default;
        Cascade(Cascade&&) = default;
        Cascade& operator=(const Cascade&) = default;
        Cascade& operator=(Cascade&&) = default;

        /// Default destructor
        ~Cascade() = default;
        ///@}

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
        void Evolve(std::shared_ptr<Nucleus>, const std::size_t& maxSteps);

        /// Simulate evolution of a kicked particle until it interacts for the 
        /// first time with another particle, accumulating the total distance
        /// traveled by the kicked particle before it interacts.
        ///@param nucleus: The nucleus to evolve according to the mean free path calculation
        ///@param maxSteps: The maximum steps to take in the particle evolution
        void MeanFreePath(std::shared_ptr<Nucleus>, const std::size_t& maxSteps);

        /// Simulate the cascade until all particles either escape, are recaptured, or are in 
        /// the background. This is done according to the NuWro algorithm.
        ///@param nucleus: The nucleus to evolve according to the NuWro method of cascade
        ///@param maxSteps: The maximum steps to take in the particle evolution
	void NuWro(std::shared_ptr<Nucleus>, const std::size_t& maxSteps);

	/// Simulate evolution of a kicked particle until it interacts for the 
        /// first time with another particle, accumulating the total distance
        /// traveled by the kicked particle before it interacts.
        ///@param nucleus: The nucleus to evolve according to the mean free path calculation
        ///@param maxSteps: The maximum steps to take in the particle evolution
        void MeanFreePath_NuWro(std::shared_ptr<Nucleus>, const std::size_t& maxSteps);
        ///@}
    private:
        // Functions
        std::size_t GetInter(Particles&, const Particle&, double& stepDistance); 
        void AdaptiveStep(const Particles&, const double&) noexcept;
        bool BetweenPlanes(const ThreeVector&, const ThreeVector&, const ThreeVector&) const noexcept;
        const ThreeVector Project(const ThreeVector&, const ThreeVector&, const ThreeVector&) const noexcept;
        const InteractionDistances AllowedInteractions(Particles&, const std::size_t&) const noexcept;
        double GetXSec(const Particle&, const Particle&) const;
        std::size_t Interacted(const Particles&, const Particle&,
                const InteractionDistances&) noexcept;
        void Escaped(Particles&);
        bool FinalizeMomentum(Particle&, Particle&) noexcept;
        bool PauliBlocking(const Particle&) const noexcept;

        // Variables
        std::vector<std::size_t> kickedIdxs;
        double distance, timeStep{};
        std::shared_ptr<Interactions> m_interactions;
        std::function<double(double, double)> probability;
        randutils::mt19937_rng rng;
        std::shared_ptr<Nucleus> localNucleus;
};

}

#endif // end of include guard: CASCADE_HH
