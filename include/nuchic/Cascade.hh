#ifndef CASCADE_HH
#define CASCADE_HH

#include <array>
#include <memory>
#include <vector>

#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/Random.hh"

class Particle;
class Interactions;

using Particles = std::vector<Particle>;
using InteractionDistances = std::vector<std::pair<std::size_t, double>>;

class Cascade {
    public:
        // Constructor and Destructor
        Cascade(const std::shared_ptr<Interactions> interactions, const double& dist = 0.05) 
            : distance(dist), m_interactions(interactions) {}
        ~Cascade() {}

        // Functions
        Particles Kick(const Particles&, const FourVector&, const std::array<double, 2>&);
        void Reset();
        Particles operator()(const Particles&, const double&, const double&,
                const std::size_t& maxSteps=10000); 

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

#endif // end of include guard: CASCADE_HH
