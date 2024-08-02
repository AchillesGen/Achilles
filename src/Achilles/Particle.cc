#include <cmath>
#include <iostream>

#include "Achilles/Constants.hh"
#include "Achilles/Particle.hh"

using namespace achilles;

void Particle::SetFormationZone(const FourVector &p1, const FourVector &p2) noexcept {
    formationZone = p1.E() / std::abs(Constant::mN * Constant::mN - p1 * p2);
}

void Particle::Propagate(const double &time) noexcept {
    const double dist = momentum.P() / momentum.E() * time * Constant::HBARC;

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist * std::sin(theta) * std::cos(phi);
    const double propDistY = dist * std::sin(theta) * std::sin(phi);
    const double propDistZ = dist * std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    if(status == ParticleStatus::internal_test) distanceTraveled += propDist.Magnitude();
    if(formationZone > 0) formationZone -= time;
    position += propDist;
}

void Particle::SpacePropagate(const double &dist) noexcept {
    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist * std::sin(theta) * std::cos(phi);
    const double propDistY = dist * std::sin(theta) * std::sin(phi);
    const double propDistZ = dist * std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    if(status == ParticleStatus::internal_test) distanceTraveled += propDist.Magnitude();
    position += propDist;
}

void Particle::BackPropagate(const double &time) noexcept {
    const double dist = momentum.P() / momentum.E() * time * Constant::HBARC;

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist * std::sin(theta) * std::cos(phi);
    const double propDistY = dist * std::sin(theta) * std::sin(phi);
    const double propDistZ = dist * std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    position -= propDist;
}

void Particle::Rotate(const std::array<double, 9> &rot_mat) noexcept {
    momentum = momentum.Rotate(rot_mat);
}

bool Particle::operator==(const Particle &other) const noexcept {
    if(info != other.info) return false;
    if(status != other.status) return false;
    if(formationZone != other.formationZone) return false;
    if(momentum != other.momentum) return false;
    if(position != other.position) return false;
    if(mothers != other.mothers) return false;
    if(daughters != other.daughters) return false;
    return true;
}

std::string Particle::ToString() const noexcept {
    return "Particle(" + std::to_string(info.IntID()) + ", " + momentum.ToString() + ", " +
           position.ToString() + ", " + std::to_string(static_cast<int>(status)) + ")";
}

namespace achilles {

std::istream &operator>>(std::istream &is, Particle &particle) {
    std::string head(9, ' '), sep1(2, ' '), sep2(2, ' '), sep3(1, ' '), tail(1, ' ');
    int pid, status;
    FourVector momentum;
    ThreeVector position;
    is.read(&head[0], 9);
    is >> pid;
    is.read(&sep1[0], 2);
    is >> momentum;
    is.read(&sep2[0], 2);
    is >> position;
    is.read(&sep3[0], 1);
    is >> status;
    is.read(&tail[0], 1);
    if(head == "Particle(" && sep1 == ", " && sep2 == ", " && sep3 == "," && tail == ")")
        particle = Particle(pid, momentum, position, status);
    return is;
}

double ClosestApproach(const Particle &particle1, const Particle &particle2) {
    auto position = particle2.Position() - particle1.Position();
    auto velocity = particle1.Beta() * Constant::HBARC;
    return position.Dot(velocity) / velocity.Magnitude2();
}

bool operator==(const std::reference_wrapper<Particle> &lhs, const Particle &rhs) {
    return lhs.get() == rhs;
}

bool operator==(const Particle &lhs, const std::reference_wrapper<Particle> &rhs) {
    return lhs == rhs.get();
}

} // namespace achilles
