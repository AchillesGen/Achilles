// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>
#include <iostream>

#include "Achilles/Particle.hh"
#include "Achilles/PhysicalUnits.hh"

using namespace achilles;

namespace {
// Displacement of length `dist` along the particle's direction of travel.
achilles::ThreePosition PropagationStep(const achilles::FourVector &momentum,
                                        achilles::units::Length dist) noexcept {
    const double theta = momentum.Theta();
    const double phi = momentum.Phi();
    return {dist * std::sin(theta) * std::cos(phi), dist * std::sin(theta) * std::sin(phi),
            dist * std::cos(theta)};
}
} // namespace

void Particle::SetFormationZone(const FourVector &p1, const FourVector &p2) noexcept {
    // E / |m_N^2 - p1.p2| is already an inverse energy, i.e. a time: in natural
    // units the hbar*c the old expression carried is 1.
    formationZone = p1.E() / units::abs(Constant::mN2() - p1 * p2);
}

void Particle::Propagate(units::Time time) noexcept {
    const units::Length dist = (momentum.P() / momentum.E()).native() * time;

    const ThreePosition propDist = PropagationStep(momentum, dist);
    if(status == ParticleStatus::internal_test) distanceTraveled += propDist.Magnitude();
    if(formationZone > units::Time{}) formationZone -= time;
    position += propDist;
}

void Particle::SpacePropagate(units::Length dist) noexcept {
    const ThreePosition propDist = PropagationStep(momentum, dist);
    if(status == ParticleStatus::internal_test) distanceTraveled += propDist.Magnitude();
    position += propDist;
}

void Particle::BackPropagate(units::Time time) noexcept {
    const units::Length dist = (momentum.P() / momentum.E()).native() * time;

    const ThreePosition propDist = PropagationStep(momentum, dist);
    if(status == ParticleStatus::internal_test) distanceTraveled -= propDist.Magnitude();
    if(formationZone != units::Time{}) formationZone += time;
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
    ThreePosition position;
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

units::Time ClosestApproach(const Particle &particle1, const Particle &particle2) {
    auto position = particle2.Position() - particle1.Position();
    auto velocity = particle1.Beta();
    return position.Dot(velocity) / velocity.Magnitude2();
}

bool operator==(const std::reference_wrapper<Particle> &lhs, const Particle &rhs) {
    return lhs.get() == rhs;
}

bool operator==(const Particle &lhs, const std::reference_wrapper<Particle> &rhs) {
    return lhs == rhs.get();
}

} // namespace achilles
