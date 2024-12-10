#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <iosfwd>
#include <utility>
#include <vector>

#include "fmt/core.h"
#include "fmt/format.h"
#include "spdlog/fmt/ostr.h"

#include "Achilles/FourVector.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/ThreeVector.hh"

namespace achilles {

// The codes given here are to match to the NuHepMC standard
enum class ParticleStatus : int {
    any = 0,
    final_state = 1,
    decayed = 2,
    initial_state = 3,
    beam = 4,
    target = 11,
    internal_test = 22,
    external_test = 23,
    propagating = 24,
    background = 25,
    captured = 26,
    residue = 27,
    spectator = 28,
    interacted = 29,
    absorption_partner = 99,
};
inline auto format_as(achilles::ParticleStatus s) {
    return fmt::underlying(s);
}

/// The Particle class provides a container to handle information about the particle.
/// The information includes the particle identification (PID), the momentum of the particle,
/// the position of the particle, the status code associated with the particle, the particle's
/// mother particles, the particle's daughter particles, and any other needed information.
/// Currently, the status codes are as follows:
///
/// @rst
///
/// +-------------+-----------------------+
/// | Status Code | Meaning               |
/// +=============+=======================+
/// |     -1      |  Propagating Particle |
/// +-------------+-----------------------+
/// |      0      |  Background Particle  |
/// +-------------+-----------------------+
/// |      1      |  Escaped Particle     |
/// +-------------+-----------------------+
///
/// @endrst
class Particle {
  public:
    /// @name Constructors and Destructors
    ///@{

    /// Create a particle
    ///@param pid: The PID of the particle (default = 0)
    ///@param momentum: The momentum of the particle (default = FourVector())
    ///@param position: The position of the particle (default = ThreeVector())
    ///@param status: The status code of the particle (default = 0)
    ///@param mothers: The mother particles of the particle (default = Empty)
    ///@param daughters: The daughter particles of the particle (default = Empty)
    explicit Particle(const PID &pid = PID{0}, FourVector mom = FourVector(),
                      ThreeVector pos = ThreeVector(),
                      const ParticleStatus &_status = ParticleStatus::background,
                      std::vector<int> _mothers = std::vector<int>(),
                      std::vector<int> _daughters = std::vector<int>()) noexcept
        : info(pid), momentum(std::move(mom)), position(std::move(pos)), status(_status),
          mothers(std::move(_mothers)), daughters(std::move(_daughters)) {
        formationZone = 0;
    }

    explicit Particle(const long int &pid, const FourVector &mom = FourVector(),
                      ThreeVector pos = ThreeVector(), const int &_status = 0,
                      std::vector<int> _mothers = std::vector<int>(),
                      std::vector<int> _daughters = std::vector<int>()) noexcept
        : info(pid), momentum(mom), position(std::move(pos)),
          status(static_cast<ParticleStatus>(_status)), mothers(std::move(_mothers)),
          daughters(std::move(_daughters)) {
        formationZone = 0;
    }

    explicit Particle(ParticleInfo _info, const FourVector &mom = FourVector(),
                      ThreeVector pos = ThreeVector(),
                      ParticleStatus _status = ParticleStatus::background,
                      std::vector<int> _mothers = std::vector<int>(),
                      std::vector<int> _daughters = std::vector<int>()) noexcept
        : info(std::move(_info)), momentum(std::move(mom)), position(std::move(pos)),
          status(std::move(_status)), mothers(std::move(_mothers)),
          daughters(std::move(_daughters)) {
        formationZone = 0;
    }

    Particle(const Particle &other)
        : info{other.info}, momentum{other.momentum}, position{other.position},
          status{other.status}, mothers{other.mothers}, formationZone{other.formationZone} {}

    Particle(Particle &&) = default;
    Particle &operator=(const Particle &) = default;
    Particle &operator=(Particle &&) = default;

    /// Default destructor
    ~Particle() = default;
    ///@}

    /// @name Setters
    /// @{
    /// These functions provide access to setting the parameters of the Particle object

    /// Set the position of the particle
    ///@param ThreeVector: The position to be set
    void SetPosition(const ThreeVector &pos) noexcept { position = pos; }

    /// Set the momentum of the particle
    ///@param FourVector: The momentum to be set
    void SetMomentum(const FourVector &mom) noexcept { momentum = mom; }

    /// Set the mother particles of the given particle
    ///@param std::vector<int>: A vector containing information about the mother particles
    void SetMothers(const std::vector<int> &_mothers) noexcept { mothers = _mothers; }

    /// Set the daughter particles of the given particle
    ///@param std::vector<int>: A vector containing information about the daughter particles
    void SetDaughters(const std::vector<int> &_daughters) noexcept { daughters = _daughters; }

    /// Add a new mother particle to an existing particle
    ///@param idx: The index of the mother particle to be set
    void AddMother(const int &idx) noexcept { mothers.push_back(idx); }

    /// Add a new daughter particle to an existing particle
    ///@param idx: The index of the daughter particle to be set
    void AddDaughter(const int &idx) noexcept { daughters.push_back(idx); }

    /// Set the formation zone of the particle. The formation zone is a time in which
    /// the particle is not allowed to interact. The formation zone is discussed in detail in:
    ///
    /// * L. Stodolsky, Formation Zone Description in Multiproduction, 1975
    /// * Phys. Rev. C. 86.015505.
    ///
    /// The equation is given by:
    /// @f[
    ///    t_{f} = \frac{E}{|p\cdot q|},
    /// @f]
    ///  where E is the energy of the incoming hadron, p and q are the momentum
    ///  of the incoming (outgoing) hadron.
    ///@param p: The momentum of the hadron before the interaction
    ///@param q: The momentum of the hadron after the interaction
    void SetFormationZone(const FourVector &, const FourVector &) noexcept;

    /// Update the formation zone time by a given time step
    ///@param timeStep: The time step to update the formation zone by
    void UpdateFormationZone(const double &timeStep) noexcept { formationZone -= timeStep; }
    ///@}

    /// @name Getters
    /// @{
    /// These functions provide get specific features from the Particle object

    /// Return the pid of the particle
    ///@return int: PID of the particle
    PID ID() const noexcept { return info.IntID(); }

    ParticleInfo Info() const noexcept { return info; }

    /// Returns the position of the particle
    ///@return ThreeVector: The position of the particle
    ThreeVector &Position() noexcept { return position; }
    const ThreeVector &Position() const noexcept { return position; }

    /// Returns the momentum of the particle
    ///@return FourVector: The momentum of the particle
    const FourVector &Momentum() const noexcept { return momentum; }
    FourVector &Momentum() noexcept { return momentum; }

    /// Gets the velocity / boost vector of a given particle
    ///@return ThreeVector: Velocity of the particle in units of c
    const ThreeVector Beta() const noexcept { return momentum.BoostVector(); }

    /// Gets the current particle status
    ///@return int: The status of the particle
    const ParticleStatus &Status() const noexcept { return status; }

    /// Set the status of the particle (See Particle description for details)
    ///@param int: The status to be set
    ParticleStatus &Status() noexcept { return status; }

    /// Return a vector of the mother particle indices
    ///@return std::vector<int>: A vector of indices referring to the mother particle
    const std::vector<int> &Mothers() const noexcept { return mothers; }

    /// Return a vector of the daughter particle indices
    ///@return std::vector<int>: A vector of indices referring to the daughter particle
    const std::vector<int> &Daughters() const noexcept { return daughters; }

    /// Return the current time remaining in the formation zone
    ///@return double: Time left in formation zone
    const double &FormationZone() const noexcept { return formationZone; }

    /// Return the mass of the given particle
    ///@return double: The mass of the particle
    double Mass() const noexcept { return info.Mass(); }

    /// Return the momentum in the x-direction
    ///@return double: Value of momentum in x-direction
    double Px() const noexcept { return momentum.Px(); }

    /// Return the momentum in the y-direction
    ///@return double: Value of momentum in y-direction
    double Py() const noexcept { return momentum.Py(); }

    /// Return the momentum in the z-direction
    ///@return double: Value of momentum in z-direction
    double Pz() const noexcept { return momentum.Pz(); }

    /// Return the energy
    ///@return double: Value of energy
    double E() const noexcept { return momentum.E(); }
    double Radius() const noexcept { return position.Magnitude(); }
    ///@}

    /// @name Functions
    /// @{

    /// Check to see if the particle is in the formation zone
    ///@return bool: True if in formation zone, False otherwise
    bool InFormationZone() const noexcept { return formationZone > 0; }

    /// Check to see if the particle is a background particle
    ///@return bool: True if a background particle, False otherwise
    bool IsBackground() const noexcept { return status == ParticleStatus::background; }

    /// Check to see if the particle is a propagating particle in the nucleus
    ///@return bool: True if a propagating particle, False otherwise
    bool IsPropagating() const noexcept {
        return status == ParticleStatus::propagating || status == ParticleStatus::external_test ||
               status == ParticleStatus::internal_test;
    }

    /// Check to see if the particle is a final state particle
    ///@return bool: True if a final state particle, False otherwise
    bool IsFinal() const noexcept { return status == ParticleStatus::final_state; }

    /// Check to see if the particle is a final state particle
    ///@return bool: True if a final state particle, False otherwise
    bool IsInitial() const noexcept { return status == ParticleStatus::initial_state; }

    /// Check to see if the particle should be written out
    ///@return bool: True if the particle should be written out, False otherwise
    bool IsExternal() const noexcept { return IsInitial() || IsFinal(); }

    /// Propagate the particle according to its momentum by a given time step
    ///@param timeStep: The amount of time to propagate the particle for
    void Propagate(const double &) noexcept;

    void SpacePropagate(const double &) noexcept;

    double &DistanceTraveled() { return distanceTraveled; }

    /// Returns the distance travelled by the particle
    /// @return double: the distance travelled by the particle
    double GetDistanceTraveled() const { return distanceTraveled; }

    /// Propagate a particle back in time. Useful for testing purposes
    ///@param timeStep: The amount of time to propagate a particle back in time for
    void BackPropagate(const double &) noexcept;

    /// Rotate the particle's spatial momentum using a 3x3 rotation matrix
    ///@param rot_mat: the 3x3 rotation matrix
    void Rotate(const std::array<double, 9> &) noexcept;

    /// Return a string representation of the particle
    ///@return std::string: a string representation of the particle
    std::string ToString() const noexcept;

    /// Determine if two particles are the same particle
    ///@param other: Particle to compare against
    ///@return bool: True if the particles are the same, otherwise False
    bool operator==(const Particle &) const noexcept;

    /// Determine if two particles are not the same particle
    ///@param other: Particle to compare against
    ///@return bool: False if the particles are the same, otherwise True
    bool operator!=(const Particle &other) const noexcept { return !(*this == other); }
    ///@}

    /// @name Stream Operators
    /// @{
    /// Stream operators for writing to and reading from streams

    /// Write out a particle to an output stream
    ///@param ostr: Output stream to write to
    ///@param part: The particle to be written out
    template <typename OStream> friend OStream &operator<<(OStream &os, const Particle &particle) {
        os << "Particle(" << particle.info.IntID() << ", " << particle.momentum << ", "
           << particle.position << ", " << static_cast<int>(particle.status) << ")";
        return os;
    }

    /// Write in a particle to an input stream
    ///@param istr: Input stream to read from
    ///@param part: The particle to be read into
    friend std::istream &operator>>(std::istream &, Particle &);
    /// @}

  private:
    ParticleInfo info;
    FourVector momentum;
    ThreeVector position;
    ParticleStatus status;
    std::vector<int> mothers, daughters;
    double formationZone;
    double distanceTraveled = 0.0;
};

/// Find the time to closest approach between two particles
/// @param p1: The first particle
/// @param p2: The second particle
/// @return double: The time to closest approach between the two particles
double ClosestApproach(const Particle &, const Particle &);

// Comparisons for std::reference_wrapper<Particle> with Particle
bool operator==(const std::reference_wrapper<Particle> &, const Particle &);
bool operator==(const Particle &, const std::reference_wrapper<Particle> &);

// Filtering functions
using vParticles = std::vector<Particle>;
using refParticles = std::vector<std::reference_wrapper<Particle>>;
using crefParticles = std::vector<std::reference_wrapper<const Particle>>;

template <class UnaryPred>
crefParticles FilterParticles(const vParticles &particles, UnaryPred pred) {
    crefParticles result;
    std::copy_if(particles.begin(), particles.end(), std::back_inserter(result), pred);
    return result;
}

template <class UnaryPred> refParticles FilterParticles(vParticles &particles, UnaryPred pred) {
    refParticles result;
    std::copy_if(particles.begin(), particles.end(), std::back_inserter(result), pred);
    return result;
}

} // namespace achilles

namespace fmt {

template <> struct formatter<achilles::Particle> {
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        return ctx.begin();
    }

    auto format(const achilles::Particle &particle, format_context &ctx) const
        -> format_context::iterator {
        return format_to(ctx.out(), "Particle[{}, {}, {}, {}]", particle.ID(), particle.Status(),
                         particle.Momentum(), particle.Position());
    }
};

template <> struct formatter<std::reference_wrapper<achilles::Particle>> {
    template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const std::reference_wrapper<achilles::Particle> &particle,
                FormatContext &ctx) const {
        return format_to(ctx.out(), "{}", particle.get());
    }
};

template <> struct formatter<achilles::ParticleStatus> {
    template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const achilles::ParticleStatus &status, FormatContext &ctx) const {
        switch(status) {
        case achilles::ParticleStatus::internal_test:
            return format_to(ctx.out(), "internal_test({})", static_cast<int>(status));
        case achilles::ParticleStatus::external_test:
            return format_to(ctx.out(), "external_test({})", static_cast<int>(status));
        case achilles::ParticleStatus::propagating:
            return format_to(ctx.out(), "propagating({})", static_cast<int>(status));
        case achilles::ParticleStatus::background:
            return format_to(ctx.out(), "background({})", static_cast<int>(status));
        case achilles::ParticleStatus::initial_state:
            return format_to(ctx.out(), "initial_state({})", static_cast<int>(status));
        case achilles::ParticleStatus::final_state:
            return format_to(ctx.out(), "final_state({})", static_cast<int>(status));
        case achilles::ParticleStatus::captured:
            return format_to(ctx.out(), "captured({})", static_cast<int>(status));
        case achilles::ParticleStatus::decayed:
            return format_to(ctx.out(), "decayed({})", static_cast<int>(status));
        case achilles::ParticleStatus::beam:
            return format_to(ctx.out(), "beam({})", static_cast<int>(status));
        case achilles::ParticleStatus::target:
            return format_to(ctx.out(), "target({})", static_cast<int>(status));
        case achilles::ParticleStatus::spectator:
            return format_to(ctx.out(), "spectator({})", static_cast<int>(status));
        case achilles::ParticleStatus::interacted:
            return format_to(ctx.out(), "cascade({})", static_cast<int>(status));
        default:
            return format_to(ctx.out(), "Unknown achilles::ParticleStatus({}) ",
                             static_cast<int>(status));
        }
    }
};

} // namespace fmt

#endif // end of include guard: PARTICLE_HH
