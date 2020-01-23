#ifndef NUCLEUS_HH
#define NUCLEUS_HH

#include <functional>
#include <iosfwd>
#include <map>
#include <memory>
#include <vector>

#include "nuchic/Random.hh"

class Particle;
class ThreeVector;

using Particles = std::vector<Particle>;

class Nucleus {
    public:
        Nucleus(const int&, const int&, const double&, const double&,
                const std::function<Particles()>&);
        ~Nucleus() {};

        // Setters
        void SetNucleons(Particles& _nucleons) noexcept;
        void SetBindingEnergy(const double& energy) noexcept {binding = energy;}
        void SetFermiMomentum(const double& mom) noexcept {fermiMomentum = mom;}
        void SetPotential(const double& energy) noexcept {potential = energy;}
        void SetDensity(const std::function<Particles()>& _density) noexcept {density = _density;}

        // Getters
        Particles Nucleons() const noexcept {return nucleons;}
        Particles Protons() const noexcept {return protons;}
        Particles Neutrons() const noexcept {return neutrons;}
        const int NNucleons() const noexcept {return nucleons.size();}
        const int NProtons() const noexcept {return protons.size();}
        const int NNeutrons() const noexcept {return neutrons.size();}
        const double& BindingEnergy() const noexcept {return binding;}
        const double& FermiMomentum() const noexcept {return fermiMomentum;}
        const double& PotentialEnergy() const noexcept {return potential;}

        // Functions
        bool Escape(Particle&) noexcept;
        Particles GenerateConfig();
        const std::array<double, 3> GenerateMomentum() noexcept;
        const std::string ToString() const noexcept;

        // Nucleus maker
        static Nucleus MakeNucleus(const std::string&, const double&,
                                   const double&, const std::function<Particles()>&);

        // Stream Operators
        friend std::ostream& operator<<(std::ostream&, const Nucleus&);
        friend std::istream& operator>>(std::istream&, Nucleus&);

    private:
        Particles nucleons, protons, neutrons;
        double binding, fermiMomentum, radius, potential;
        std::function<Particles()> density;

        static const std::map<int, std::string> ZToName;
        static const std::map<std::string, int> NameToZ;

        randutils::mt19937_rng rng;
};



#endif // end of include guard: NUCLEUS_HH
