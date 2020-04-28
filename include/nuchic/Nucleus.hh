#ifndef NUCLEUS_HH
#define NUCLEUS_HH

#include <functional>
#include <iosfwd>
#include <map>
#include <memory>
#include <vector>

#include "nuchic/Constants.hh"
#include "nuchic/Interpolation.hh"
#include "nuchic/Random.hh"

namespace nuchic {

class Particle;
class ThreeVector;

using Particles = std::vector<Particle>;

/// The Nucleus class implements the physics needed to describe an arbitrary nucleus. It provides
/// the ability to generate configurations of nucleons for the cascade, as well as perform checks
/// on if nucleons are captured in the potential or escape.
class Nucleus {
    public:
        // Fermigas Model
        enum FermiGasType {
            Local,
            Global
        };

        /// @name Constructors and Destructors
        /// @{

        /// Create a nucleus
        ///@param Z: The number of protons
        ///@param A: The number of nucleons
        ///@param binding: The binding energy of the nucleus
        ///@param densityFile: The file containing density information for Pauli Blocking
        ///TODO: This should be added to the Nucleus class when we refactor to have the Nucleus
        ///      passed in as an object
        ///@param density: A function that generates nucleon configurations according 
        ///                to the density profile
        Nucleus(const std::size_t&, const std::size_t&, const double&, const double&,
                const std::string&, const FermiGasType&, std::function<Particles()>);
        Nucleus(const Nucleus&) = default;
        Nucleus(Nucleus&&) = default;
        Nucleus& operator=(const Nucleus&) = default;
        Nucleus& operator=(Nucleus&&) = default;

        /// Default destructor
        ~Nucleus() = default;
        ///@}

        /// @name Setters
        /// @{
        /// These functions provide access to setting the parameters of the Nucleus object

        /// Set the nucleons, protons, and neutrons of the nucleus
        ///@param nucleons: The nucleons to be used for the nucleus
        void SetNucleons(Particles& _nucleons) noexcept;

        /// Set the binding energy of the nucleus in MeV
        ///@param energy: The binding energy to be set in MeV
        void SetBindingEnergy(const double& energy) noexcept { binding = energy; }

        /// Set the Fermi Momentum of the nucleus in MeV
        ///@param mom: The Fermi Momentum to be set in MeV
        void SetFermiMomentum(const double& mom) noexcept { fermiMomentum = mom; }

        /// Set the potential energy of the nucleus. This value defaults to:
        /// @f[
        ///     \sqrt{m_{N}^2+k_{f}^2} - m_{N} + 8 \text{MeV}.
        /// @f]
        ///@param energy: The potential energy to be set in MeV
        void SetPotential(const double& energy) noexcept { potential = energy; }

        /// Set the density function to use for configuration generation
        ///@param density: The function to be use for generating nucleons
        void SetDensity(const std::function<Particles()>& _density) noexcept { density = _density; }
        ///@}

        /// Set the radius of the nucleus in fm
        /// @param radius: The radius to be set in fm
        void SetRadius(const double& _radius) noexcept { radius = _radius; }

        /// @name Getters
        /// @{
        /// These functions provide get specific features from the Nucleus object

        /// Return a vector of the current nucleons
        ///@return Particles: The current nucleons generated for the nucleus
        Particles& Nucleons() noexcept { return nucleons; }

        /// Return a vector of the current protons
        ///@return Particles: The current protons generated for the nucleus
        Particles& Protons() noexcept { return protons; }

        /// Return a vector of the current neutrons
        ///@return Particles: The current neutrons generated for the nucleus
        Particles& Neutrons() noexcept { return neutrons; }

        /// Return the number of nucleons in the nucleus
        ///@return int: The number of nucleons in the nucleus
        std::size_t NNucleons() const noexcept { return nucleons.size(); }

        /// Return the number of protons in the nucleus
        ///@return int: The number of protons in the nucleus
        std::size_t NProtons() const noexcept { return protons.size(); }

        /// Return the number of neutrons in the nucleus
        ///@return int: The number of neutrons in the nucleus
        std::size_t NNeutrons() const noexcept { return neutrons.size(); }

        /// Return the current binding energy of the nucleus
        ///@return double: The binding energy in MeV
        const double& BindingEnergy() const noexcept { return binding; }

	/// Return the current Fermi Momentum of the nucleus
        ///@return double: The Fermi Momentum in MeV
        //const double& FermiMomentum() const noexcept {return fermiMomentum;}

        /// Return the phenomenological potential
	///@return double: The potential in MeV	
        double Potential(const double&) const noexcept;

        /// Return the current potential energy of the nucleus
        ///@return double: The potential energy in MeV
        const double& PotentialEnergy() const noexcept { return potential; }

        /// Return the radius cutoff of the nucleus used for the cascade
        ///@return double: The radius in femtometers
        const double& Radius() const noexcept { return radius; }

        /// Return the density of the nucleus at a given location
        ///@param position: The radius to calculate the density at
        ///@return double: The density at the input radius
        double Rho(const double &position) const noexcept { return rhoInterp(position); }
        ///@}
	
        /// Return the Fermi momentum according to a given FG model
	///@param position: The radius to calculate the density
        double FermiMomentum(const double&) const noexcept;	//
	///@}

        /// @name Functions
        /// @{

        /// Determine if a particle escapes from the nucleus or is recaptured
        ///@param particle: The particle that is attempting to escape
        ///@return bool: True if the particle escapes, False if recaptured
        bool Escape(Particle&) noexcept;

        /// Generate a configuration of the nucleus based on the density function
        void GenerateConfig();

        /// Generate a random momentum for a nucleon in the nucleus
        ///@return std::array<double, 3>: Random momentum generated using the Fermi momentum
        const std::array<double, 3> GenerateMomentum(const double&) noexcept;

        /// Return a string representation of the nucleus
        ///@return std::string: a string representation of the nucleus
        const std::string ToString() const noexcept;
        /// @}

        // Nucleus maker
        /// Generate a nucleus based on an input string, a binding energy, a Fermi Momentum,
        /// and a density function. The format of the string has to be of the form, AX, where
        /// A represents the number of nucleons, and X represents the chemical symbol for the 
        /// element. Currently, the following elements are available in the code:
        ///
        /// @rst
        /// 
        /// +-----+-----+
        /// |  Z  |  X  |
        /// +=====+=====+
        /// |  1  |  H  |
        /// +-----+-----+
        /// |  2  | He  |
        /// +-----+-----+
        /// |  3  | Li  |
        /// +-----+-----+
        /// |  6  |  C  |
        /// +-----+-----+
        /// |  8  |  O  |
        /// +-----+-----+
        /// | 13  | Al  |
        /// +-----+-----+
        /// | 18  | Ar  |
        /// +-----+-----+
        /// | 20  | Ca  |
        /// +-----+-----+
        /// | 26  | Fe  |
        /// +-----+-----+
        ///
        /// @endrst
        ///@param name: Nucleus in format discussed above
        ///@param binding: The binding energy of the nucleus
        ///@param densityFile: The file containing density information for Pauli Blocking
        ///TODO: This should be added to the Nucleus class when we refactor to have the Nucleus
        ///      passed in as an object
        ///@param density: The density function to use to generate configurations with
        static Nucleus MakeNucleus(const std::string&, const double&, const double&,
                                   const std::string&, const FermiGasType&,
                                   const std::function<Particles()>&);

        /// @name Stream Operators
        /// @{
        /// These functions are currently **not implemented** for this class

        friend std::ostream& operator<<(std::ostream&, const Nucleus&);
        friend std::istream& operator>>(std::istream&, Nucleus&);
        /// @}

    private:
        Particles nucleons, protons, neutrons;
        double binding, fermiMomentum, radius, potential{};
        FermiGasType fermiGas;
        std::function<Particles()> density;
        Interp1D rhoInterp;	

        static const std::map<std::size_t, std::string> ZToName;
        static std::size_t NameToZ(const std::string&);

        randutils::mt19937_rng rng;
};

}

#endif // end of include guard: NUCLEUS_HH
