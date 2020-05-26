#include <iostream>
#include <cmath>

#include "nuchic/Constants.hh"
#include "nuchic/Particle.hh"

using namespace nuchic;

void Particle::SetFormationZone(const FourVector& p1, const FourVector& p2) noexcept {
    formationZone = p1.E()/std::abs(Constant::mN*Constant::mN-p1*p2);
}

void Particle::Propagate(const double& time) noexcept {
    const double dist = momentum.P()/momentum.E()*time*Constant::HBARC;

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist*std::sin(theta)*std::cos(phi);
    const double propDistY = dist*std::sin(theta)*std::sin(phi);
    const double propDistZ = dist*std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    if (status == ParticleStatus::internal_test) distanceTraveled += propDist.Magnitude();
    position += propDist;
}


void Particle::SpacePropagate(const double& dist) noexcept {

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist*std::sin(theta)*std::cos(phi);
    const double propDistY = dist*std::sin(theta)*std::sin(phi);
    const double propDistZ = dist*std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    if (status == ParticleStatus::internal_test) distanceTraveled += propDist.Magnitude();
    position += propDist;
}


void Particle::BackPropagate(const double& time) noexcept {
    const double dist = momentum.P()/momentum.E()*time*Constant::HBARC;

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist*std::sin(theta)*std::cos(phi);
    const double propDistY = dist*std::sin(theta)*std::sin(phi);
    const double propDistZ = dist*std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    position -= propDist;
}

bool Particle::operator==(const Particle& other) const noexcept {
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
    return "Particle(" + std::to_string(info.IntID()) + ", " + momentum.ToString()
        + ", " + position.ToString() + ", " + std::to_string(static_cast<int>(status)) + ")";
}

namespace nuchic {

std::istream& operator>>(std::istream& is, Particle& particle) {
    std::string head(9, ' '), sep1(1, ' '), sep2(1, ' '), sep3(1, ' '), tail(1, ' ');
    int pid, status;
    FourVector momentum;
    ThreeVector position;
    is.read(&head[0], 9);
    is >> pid;
    is.read(&sep1[0], 1);
    is >> momentum;
    is.read(&sep2[0], 1);
    is >> position;
    is.read(&sep3[0], 1);
    is >> status;
    is.read(&tail[0], 1);
    if(head == "Particle(" &&
       sep1 == "," && sep2 == "," &&
       sep3 == "," && tail == ")")
        particle = Particle(pid, momentum, position, status);
    return is;
}

}
