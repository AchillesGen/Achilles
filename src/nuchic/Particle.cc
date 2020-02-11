#include <iostream>
#include <cmath>

#include "nuchic/Particle.hh"

#define mN 938
#define HBARC 200

void nuchic::Particle::SetFormationZone(const nuchic::FourVector& p1, const nuchic::FourVector& p2) noexcept {
    formationZone = p1.E()/std::abs(mN*mN-p1*p2);
}

void nuchic::Particle::Propagate(const double& time) noexcept {
    const double dist = momentum.P()/momentum.E()*time*HBARC;

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist*std::sin(theta)*std::cos(phi);
    const double propDistY = dist*std::sin(theta)*std::sin(phi);
    const double propDistZ = dist*std::cos(theta);

    const nuchic::ThreeVector propDist(propDistX, propDistY, propDistZ);

    position += propDist;
}

void nuchic::Particle::BackPropagate(const double& time) noexcept {
    const double dist = momentum.P()/momentum.E()*time*HBARC;

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist*std::sin(theta)*std::cos(phi);
    const double propDistY = dist*std::sin(theta)*std::sin(phi);
    const double propDistZ = dist*std::cos(theta);

    const nuchic::ThreeVector propDist(propDistX, propDistY, propDistZ);
    position -= propDist;
}

bool nuchic::Particle::operator==(const nuchic::Particle& other) const noexcept {
    if(pid != other.pid) return false;
    if(status != other.status) return false;
    if(formationZone != other.formationZone) return false;
    if(momentum != other.momentum) return false;
    if(position != other.position) return false;
    if(mothers != other.mothers) return false;
    if(daughters != other.daughters) return false;
    return true;
}

const std::string nuchic::Particle::ToString() const noexcept {
    return "Particle(" + std::to_string(pid) + ", " + momentum.ToString()
        + ", " + position.ToString() + ", " + std::to_string(status) + ")";
}

namespace nuchic {

std::ostream& operator<<(std::ostream& os, const nuchic::Particle& particle) {
    os << "Particle(" << particle.pid << ", " << particle.momentum << ", " 
       << particle.position << ", " << particle.status << ")";
    return os;
}

std::istream& operator>>(std::istream& is, nuchic::Particle& particle) {
    std::string head(9, ' '), sep1(1, ' '), sep2(1, ' '), sep3(1, ' '), tail(1, ' ');
    int pid, status;
    nuchic::FourVector momentum;
    nuchic::ThreeVector position;
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
        particle = nuchic::Particle(pid, momentum, position, status);
    return is;
}

}
