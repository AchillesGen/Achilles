#ifndef CASCADE_HH
#define CASCADE_HH

#include <array>
#include <memory>
#include <vector>

#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Random.hh"

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
        /// @name Constructor and Destructor
        ///@{

        /// Create the Cascade object
        ///@param interactions: The interaction model for pp, pn, and np interactions
        ///@param dist: The maximum distance step to take when propagating
        Cascade(const std::shared_ptr<Interactions> interactions, const double& dist = 0.05)
            : distance(dist), m_interactions(interactions) {}

        /// Default destructor
        ~Cascade() {}
        ///@}

        /// @name Functions
        ///@{

        /// Give a random particle in the nucleus a kick. The kick is defined by a input
        /// four vector. To determine whether a proton or a neutron is kicked, the cross-section
        /// for protons and neutrons needs to be supplied.
        ///@param particles: The list of particles inside the nucleus
        ///@param energyTransfer: The energy transfered to the nucleon during the kick
        ///@param sigma: An array representing the cross-section for different kicked nucleons
        Particles Kick(const Particles&, const FourVector&, const std::array<double, 2>&);

        /// Reset the cascade internal variables for the next cascade
        void Reset();

        /// Simulate the cascade until all particles either escape, are recaptured, or are in 
        /// the background.
        ///@param particles: The list of particles in the nucleus
        ///@param kf: The Fermi Momentum to use for Pauli Blocking
        ///@param radius2: The squared radius denoting the edge of the nucleus
        ///@param maxSteps: The maximum steps to take in the cascade
        Particles operator()(const Particles&, const double&, const double&,
                const std::size_t& maxSteps=1000000);

        /// Helper function to make a specific nucleon as the kicked nucleon
        ///@param idx: The index of the particle that has been kicked
        void SetKicked(const std::size_t& idx) {kickedIdxs.push_back(idx);}

        Particles Evolve(const Nucleus& nuc, const std::size_t& maxSteps);
        ///@}

    private:
        // Functions
        void AdaptiveStep(const Particles&, const double&) noexcept;
        bool BetweenPlanes(const ThreeVector&, const ThreeVector&, const ThreeVector&) const noexcept;
        const ThreeVector Project(const ThreeVector&, const ThreeVector&, const ThreeVector&) const noexcept;
        const InteractionDistances AllowedInteractions(Particles&, const std::size_t&) const noexcept;
        const double GetXSec(const Particle&, const Particle&) const;
        int Interacted(const Particles&, const Particle&,
                const InteractionDistances&) noexcept;
        bool FinalizeMomentum(Particle&, Particle&) noexcept;
        bool PauliBlocking(const FourVector&) const noexcept;

        // Variables
        std::vector<std::size_t> kickedIdxs;
        double distance, timeStep, fermiMomentum;
        std::shared_ptr<Interactions> m_interactions;
        randutils::mt19937_rng rng;
};

}

#endif // end of include guard: CASCADE_HH
