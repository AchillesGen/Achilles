#ifndef NUCLEUS_HH
#define NUCLEUS_HH

#include <filesystem>
#include <functional>
#include <iosfwd>
#include <map>
#include <memory>
#include <vector>

#include "Achilles/Achilles.hh"
#include "Achilles/Configuration.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Potential.hh"
#include "Achilles/Random.hh"
#include "Achilles/System.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace fs = std::filesystem;

namespace achilles {

class PID;
class Particle;
class ThreeVector;

using Particles = std::vector<Particle>;

/// The Nucleus class implements the physics needed to describe an arbitrary nucleus. It provides
/// the ability to generate configurations of nucleons for the cascade, as well as perform checks
/// on if nucleons are captured in the potential or escape.
class Nucleus {
  public:
    // Fermigas Model
    enum FermiGasType { Local, Global };
    struct FermiGas {
        FermiGasType type{FermiGasType::Local};
        bool correlated = false;
        double SRCfraction{0.2}, lambdaSRC{2.75};
    };

    /// @name Constructors and Destructors
    /// @{

    /// Create a nucleus
    ///@param Z: The number of protons
    ///@param A: The number of nucleons
    ///@param binding: The binding energy of the nucleus
    ///@param densityFile: The file containing density information for Pauli Blocking
    /// TODO: This should be added to the Nucleus class when we refactor to have the Nucleus
    ///      passed in as an object
    ///@param density: A function that generates nucleon configurations according
    ///                to the density profile
    Nucleus() = default;
    Nucleus(const std::size_t &, const std::size_t &, const double &, const double &,
            const std::string &, const std::string &, const FermiGas &, std::unique_ptr<Density>);
    Nucleus(const Nucleus &) = delete;
    Nucleus(Nucleus &&) = default;
    Nucleus &operator=(const Nucleus &) = delete;
    Nucleus &operator=(Nucleus &&) = default;

    /// Default destructor
    MOCK ~Nucleus() = default;
    ///@}

    /// @name Setters
    /// @{
    /// These functions provide access to setting the parameters of the Nucleus object

    /// Set the binding energy of the nucleus in MeV
    ///@param energy: The binding energy to be set in MeV
    void SetBindingEnergy(const double &energy) noexcept { binding = energy; }

    /// Set the Fermi Momentum of the nucleus in MeV
    ///@param mom: The Fermi Momentum to be set in MeV
    void SetFermiMomentum(const double &mom) noexcept { fermiMomentum = mom; }

    /// Set the density function to use for configuration generation
    ///@param density: The function to be use for generating nucleons
    void SetDensity(std::unique_ptr<Density> _density) noexcept { density = std::move(_density); }

    /// Set the potential function used for propagation in the nucleus
    ///@param potential: The potential to use
    void SetPotential(std::unique_ptr<Potential> _potential) noexcept {
        potential = std::move(_potential);
    }

    ///@}

    /// Set the radius of the nucleus in fm
    /// @param radius: The radius to be set in fm
    void SetRadius(const double &_radius) noexcept { radius = _radius; }

    /// Setup basic properties for nuclei, needed to setup hydrogen correctly
    ///@param protons: Number of protons
    ///@param nucleons: Number of nucleons
    void Initialize(size_t protons, size_t nucleons);

    /// @name Getters
    /// @{
    /// These functions provide get specific features from the Nucleus object

    /// Return the nuclear PID
    ///@return PID: The pid for the given nucleus
    MOCK PID ID() const { return m_pid; }

    /// Return a particle with initial momentum for at rest nucleus and correct PID
    ///@return Particle: Particle representing the initial nucleus
    Particle InitParticle() const {
        double mass = ParticleInfo(ID()).Mass();
        return Particle(ID(), {mass, 0, 0, 0}, {}, ParticleStatus::target);
    }

    /// Return the number of nucleons in the nucleus
    ///@return int: The number of nucleons in the nucleus
    MOCK std::size_t NNucleons() const noexcept { return nnucleons; }

    /// Return the number of protons in the nucleus
    ///@return int: The number of protons in the nucleus
    MOCK std::size_t NProtons() const noexcept { return nprotons; }

    /// Return the number of neutrons in the nucleus
    ///@return int: The number of neutrons in the nucleus
    MOCK std::size_t NNeutrons() const noexcept { return nneutrons; }

    /// Return the current binding energy of the nucleus
    ///@return double: The binding energy in MeV
    const double &BindingEnergy() const noexcept { return binding; }

    /// Return the current Fermi Momentum of the nucleus
    ///@return double: The Fermi Momentum in MeV
    const double &FermiMomentum() const noexcept { return fermiMomentum; }

    /// Return the phenomenological potential
    ///@return std::shared_ptr<Potential>: The potential of the nucleus
    MOCK std::shared_ptr<Potential> GetPotential() const noexcept { return potential; }

    /// Return the radius cutoff of the nucleus used for the cascade
    ///@return double: The radius in femtometers
    MOCK const double &Radius() const noexcept { return radius; }

    /// Return the density of the nucleus at a given location
    ///@param position: The radius to calculate the density at
    ///@return double: The density at the input radius
    MOCK double Rho(const double &position) const noexcept {
        // return position > rhoInterp.max() ? 0 : rhoInterp(position);
        return position > radius ? 0 : protonrhoInterp(position);
    }
    ///@}

    /// Return the density of the nucleus at a given location
    ///@param position: The radius to calculate the density at
    ///@return double: The density at the input radius
    MOCK double ProtonRho(const double &position) const noexcept {
        // return position > rhoInterp.max() ? 0 : rhoInterp(position);
        return position > radius ? 0 : protonrhoInterp(position);
    }
    ///@}

    /// Return the density of the nucleus at a given location
    ///@param position: The radius to calculate the density at
    ///@return double: The density at the input radius
    MOCK double NeutronRho(const double &position) const noexcept {
        // return position > rhoInterp.max() ? 0 : rhoInterp(position);
        return position > radius ? 0 : neutronrhoInterp(position);
    }
    ///@}

    /// Return the Fermi momentum according to a given FG model
    ///@param position: The radius to calculate the density
    double FermiMomentum(const double &) const noexcept;

    void SetRecoil(const FourVector recoil) { m_recoil = recoil; }
    ///@}

    /// @name Functions
    /// @{

    /// Generate a configuration of the nucleus based on the density function
    MOCK Particles GenerateConfig();

    /// Generate a random momentum for a nucleon in the nucleus
    ///@return std::array<double, 3>: Random momentum generated using the Fermi momentum
    const std::array<double, 3> GenerateMomentum(const double &) noexcept;

    double SampleMagnitudeMomentum(const double &position) noexcept;

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
    /// TODO: This should be added to the Nucleus class when we refactor to have the Nucleus
    ///      passed in as an object
    ///@param density: The density function to use to generate configurations with
    static Nucleus MakeNucleus(const std::string &, const double &, const double &,
                               const std::string &, const std::string &, const FermiGas &,
                               std::unique_ptr<Density>);

    /// @name Stream Operators
    /// @{
    /// These functions are currently **not implemented** for this class

    friend std::ostream &operator<<(std::ostream &, const Nucleus &);
    friend std::istream &operator>>(std::istream &, Nucleus &);
    /// @}

  private:
    size_t nnucleons, nprotons, nneutrons;
    double binding{}, fermiMomentum{}, radius{};
    FermiGas fermi_gas{};
    std::unique_ptr<Density> density;
    Interp1D protonrhoInterp;
    Interp1D neutronrhoInterp;

    static const std::map<std::size_t, std::string> ZToName;
    static std::size_t NameToZ(const std::string &);
    PID m_pid;

    FourVector m_recoil{};
    std::shared_ptr<Potential> potential;
    bool is_hydrogen{false};
};

} // namespace achilles

namespace YAML {
template <> struct convert<achilles::Nucleus::FermiGas> {
    static bool decode(const Node &node, achilles::Nucleus::FermiGas &fg) {
        if(node["Type"].as<std::string>() == "Local")
            fg.type = achilles::Nucleus::FermiGasType::Local;
        else if(node["Type"].as<std::string>() == "Global")
            fg.type = achilles::Nucleus::FermiGasType::Global;
        else
            return false;

        if(node["Correlated"]) fg.correlated = node["Correlated"].as<bool>();
        if(fg.correlated) {
            if(node["SRCfraction"]) fg.SRCfraction = node["SRCfraction"].as<double>();
            if(node["LambdaSRC"]) fg.lambdaSRC = node["LambdaSRC"].as<double>();
        }

        return true;
    }
};

template <> struct convert<achilles::Nucleus> {
    static bool decode(const Node &node, achilles::Nucleus &nuc) {
        std::string name = node["Name"].as<std::string>();

        // Special case for hydrogen
        if(name == "1H") {
            spdlog::info("Nucleus: creating hydrogen nucleus.");
            nuc.Initialize(1, 1);
            return true;
        }

        auto binding = node["Binding"].as<double>();
        auto kf = node["Fermi Momentum"].as<double>();

        auto fermi_gas = node["FermiGas"].as<achilles::Nucleus::FermiGas>();
        auto protondensityFile = node["Density"]["ProtonFile"].as<std::string>();
        auto neutrondensityFile = node["Density"]["NeutronFile"].as<std::string>();
        auto configs = std::make_unique<achilles::DensityConfiguration>(
            node["Density"]["Configs"].as<std::string>());
        nuc = achilles::Nucleus::MakeNucleus(name, binding, kf, protondensityFile,
                                             neutrondensityFile, fermi_gas, std::move(configs));

        return true;
    }
};

template <> struct convert<std::vector<std::shared_ptr<achilles::Nucleus>>> {
    static bool decode(const Node &node, std::vector<std::shared_ptr<achilles::Nucleus>> &nuclei) {
        for(YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
            YAML::Node nucleus_node = (*it)["Nucleus"];
            auto nucleus =
                std::make_shared<achilles::Nucleus>(nucleus_node.as<achilles::Nucleus>());

            if(nucleus->ID() != achilles::PID::hydrogen()) {
                // Set potential for the nucleus
                auto potential_name = nucleus_node["Potential"]["Name"].as<std::string>();
                auto potential = achilles::PotentialFactory::Initialize(potential_name, nucleus,
                                                                        nucleus_node["Potential"]);
                nucleus->SetPotential(std::move(potential));
            }

            nuclei.push_back(nucleus);
        }
        return true;
    }
};
} // namespace YAML

template <> struct std::hash<achilles::Nucleus> {
    std::size_t operator()(const achilles::Nucleus &nuc) const {
        return std::hash<int>{}(nuc.ID());
    }
};

#endif // end of include guard: NUCLEUS_HH
