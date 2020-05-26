#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <iosfwd>
#include <utility>
#include <vector>

#include "spdlog/fmt/ostr.h"

#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ParticleInfo.hh"

namespace nuchic {

enum class ParticleStatus : int {
    internal_test = -3,
    external_test = -2,
    propagating = -1,
    background = 0,
    escaped = 1
};

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
        Particle(const PID& pid = PID{0}, const FourVector& mom = FourVector(),
                 ThreeVector  pos = ThreeVector(), const ParticleStatus& _status = static_cast<ParticleStatus>(0),
                 std::vector<int>  _mothers = std::vector<int>(),
                 std::vector<int>  _daughters = std::vector<int>()) noexcept :
            info(pid), momentum(mom), position(std::move(pos)), status(_status),
            mothers(std::move(_mothers)), daughters(std::move(_daughters)) { formationZone = 0;}
        Particle(const long int& pid, const FourVector& mom = FourVector(),
                 ThreeVector  pos = ThreeVector(), const int& _status = 0,
                 std::vector<int>  _mothers = std::vector<int>(),
                 std::vector<int>  _daughters = std::vector<int>()) noexcept :
            info(pid), momentum(mom), position(std::move(pos)), status(static_cast<ParticleStatus>(_status)),
            mothers(std::move(_mothers)), daughters(std::move(_daughters)) {formationZone = 0;}
        Particle(const Particle&) = default;
        Particle(Particle&&) = default;
        Particle& operator=(const Particle&) = default;
        Particle& operator=(Particle&&) = default;

        /// Default destructor
        ~Particle() = default;
        ///@}

        /// @name Setters
        /// @{
        /// These functions provide access to setting the parameters of the Particle object

        /// Set the position of the particle
        ///@param ThreeVector: The position to be set
        void SetPosition(const ThreeVector& pos) noexcept {position = pos;}

        /// Set the momentum of the particle
        ///@param FourVector: The momentum to be set
        void SetMomentum(const FourVector& mom) noexcept {momentum = mom;}

        /// Set the status of the particle (See Particle description for details)
        ///@param int: The status to be set
        void SetStatus(const ParticleStatus& _status) noexcept {status = _status;}

        /// Set the mother particles of the given particle
        ///@param std::vector<int>: A vector containing information about the mother particles
        void SetMothers(const std::vector<int>& _mothers) noexcept {mothers = _mothers;}

        /// Set the daughter particles of the given particle
        ///@param std::vector<int>: A vector containing information about the daughter particles
        void SetDaughters(const std::vector<int>& _daughters) noexcept {daughters = _daughters;}

        /// Add a new mother particle to an existing particle
        ///@param idx: The index of the mother particle to be set
        void AddMother(const int& idx) noexcept {mothers.push_back(idx);}

        /// Add a new daughter particle to an existing particle
        ///@param idx: The index of the daughter particle to be set
        void AddDaughter(const int& idx) noexcept {daughters.push_back(idx);}

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
        void SetFormationZone(const FourVector&, const FourVector&) noexcept;

        /// Update the formation zone time by a given time step
        ///@param timeStep: The time step to update the formation zone by
        void UpdateFormationZone(const double& timeStep) noexcept {formationZone -= timeStep;}
        ///@}

        /// @name Getters
        /// @{
        /// These functions provide get specific features from the Particle object

        /// Return the pid of the particle
        ///@return int: PID of the particle
        PID ID() const noexcept { return info.ID(); }

        /// Returns the position of the particle
        ///@return ThreeVector: The position of the particle
        const ThreeVector& Position() const noexcept {return position;}

        /// Returns the momentum of the particle
        ///@return FourVector: The momentum of the particle
        const FourVector& Momentum() const noexcept {return momentum;}

        /// Gets the velocity / boost vector of a given particle
        ///@return ThreeVector: Velocity of the particle in units of c
        const ThreeVector Beta() const noexcept {return momentum.BoostVector();}

        /// Gets the current particle status
        ///@return int: The status of the particle
        const ParticleStatus& Status() const noexcept {return status;}

        /// Return a vector of the mother particle indices
        ///@return std::vector<int>: A vector of indices referring to the mother particle
        const std::vector<int>& Mothers() const noexcept {return mothers;}

        /// Return a vector of the daughter particle indices
        ///@return std::vector<int>: A vector of indices referring to the daughter particle
        const std::vector<int>& Daughters() const noexcept {return daughters;}

        /// Return the current time remaining in the formation zone
        ///@return double: Time left in formation zone
        const double& FormationZone() const noexcept {return formationZone;}

        /// Return the mass of the given particle
        ///@return double: The mass of the particle
        double Mass() const noexcept {return momentum.M();}

        /// Return the momentum in the x-direction
        ///@return double: Value of momentum in x-direction
        double Px() const noexcept {return momentum.Px();}

        /// Return the momentum in the y-direction
        ///@return double: Value of momentum in y-direction
        double Py() const noexcept {return momentum.Py();}

        /// Return the momentum in the z-direction
        ///@return double: Value of momentum in z-direction
        double Pz() const noexcept {return momentum.Pz();}

        /// Return the energy
        ///@return double: Value of energy
        double E() const noexcept {return momentum.E();}
        double Radius() const noexcept {return position.Magnitude();}
        ///@}

        /// @name Functions
        /// @{

        /// Check to see if the particle is in the formation zone
        ///@return bool: True if in formation zone, False otherwise
        bool InFormationZone() const noexcept {return formationZone > 0;}

        /// Check to see if the particle is a background particle
        ///@return bool: True if a background particle, False otherwise
        bool IsBackground() const noexcept {return status == ParticleStatus::background;}

        /// Check to see if the particle is a propagating particle in the nucleus
        ///@return bool: True if a propagating particle, False otherwise
        bool IsPropagating() const noexcept {return status == ParticleStatus::propagating;}

        /// Check to see if the particle is a final state particle
        ///@return bool: True if a final state particle, False otherwise
        bool IsFinal() const noexcept {return status == ParticleStatus::escaped;}

        /// Propagate the particle according to its momentum by a given time step
        ///@param timeStep: The amount of time to propagate the particle for
        void Propagate(const double&) noexcept;

	void SpacePropagate(const double&) noexcept; 

        /// Returns the distance travelled by the particle
        /// @return double: the distance travelled by the particle
        double GetDistanceTraveled() {return distanceTraveled;}

        /// Propagate a particle back in time. Useful for testing purposes
        ///@param timeStep: The amount of time to propagate a particle back in time for
        void BackPropagate(const double&) noexcept;

        /// Return a string representation of the particle
        ///@return std::string: a string representation of the particle
        std::string ToString() const noexcept;

        /// Determine if two particles are the same particle
        ///@param other: Particle to compare against
        ///@return bool: True if the particles are the same, otherwise False
        bool operator==(const Particle&) const noexcept;

        /// Determine if two particles are not the same particle
        ///@param other: Particle to compare against
        ///@return bool: False if the particles are the same, otherwise True
        bool operator!=(const Particle& other) const noexcept {return !(*this == other);}
        ///@}

        /// @name Stream Operators
        /// @{
        /// Stream operators for writing to and reading from streams

        /// Write out a particle to an output stream
        ///@param ostr: Output stream to write to
        ///@param part: The particle to be written out
        template<typename OStream>
        friend OStream& operator<<(OStream &os, const Particle &particle) {
            os << "Particle(" << particle.info.IntID() << ", " << particle.momentum << ", "
               << particle.position << ", " << static_cast<int>(particle.status) << ")";
            return os;
        }

        /// Write in a particle to an input stream
        ///@param istr: Input stream to read from
        ///@param part: The particle to be read into
        friend std::istream& operator>>(std::istream&, Particle&);
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

}

#endif // end of include guard: PARTICLE_HH
