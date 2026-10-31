// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>
#include <iostream>

#include "Achilles/Constants.hh"
#include "Achilles/Particle.hh"

using namespace achilles;

void Particle::SetFormationZone(const FourVector &p1, const FourVector &p2) noexcept {
    formationZone = p1.E() / std::abs(Constant::mN * Constant::mN - p1 * p2) * Constant::HBARC;
}

void Particle::Propagate(const double &time) noexcept {
    const double dist = momentum.P() / momentum.E() * time;

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist * std::sin(theta) * std::cos(phi);
    const double propDistY = dist * std::sin(theta) * std::sin(phi);
    const double propDistZ = dist * std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    if(status == ParticleStatus::internal_test) distanceTraveled += propDist.Magnitude();
    if(formationZone > 0) formationZone -= time;
    Move(propDist);
}

void Particle::SpacePropagate(const double &dist) noexcept {
    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist * std::sin(theta) * std::cos(phi);
    const double propDistY = dist * std::sin(theta) * std::sin(phi);
    const double propDistZ = dist * std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    if(status == ParticleStatus::internal_test) distanceTraveled += propDist.Magnitude();
    Move(propDist);
}

void Particle::BackPropagate(const double &time) noexcept {
    const double dist = momentum.P() / momentum.E() * time;

    const double theta = momentum.Theta();
    const double phi = momentum.Phi();

    const double propDistX = dist * std::sin(theta) * std::cos(phi);
    const double propDistY = dist * std::sin(theta) * std::sin(phi);
    const double propDistZ = dist * std::cos(theta);

    const ThreeVector propDist(propDistX, propDistY, propDistZ);
    if(status == ParticleStatus::internal_test) distanceTraveled -= propDist.Magnitude();
    if(formationZone != 0) formationZone += time;
    Move(-propDist);
}

void Particle::Rotate(const std::array<double, 9> &rot_mat) noexcept {
    momentum = momentum.Rotate(rot_mat);
	for(Particle& part:sub_parts)
		part.Rotate(rot_mat);
}
void Particle::RotateBack(const std::array<double, 9> &rot_mat) noexcept {
    momentum = momentum.RotateBack(rot_mat);
	for(Particle& part:sub_parts)
		part.RotateBack(rot_mat);
}

bool Particle::equalsIgnoreVertices(const Particle &other) const noexcept {
    if(info != other.info) return false;
    if(status != other.status) return false;
    if(formationZone != other.formationZone) return false;
    if(momentum != other.momentum) return false;
    if(position != other.position) return false;
    if(mothers != other.mothers) return false;
    return true;
}

bool Particle::operator==(const Particle& other) const noexcept {
	return equalsIgnoreVertices(other)
		&&start_vertex==other.start_vertex
		&&end_vertex==other.end_vertex;
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
    auto velocity = particle1.Beta();
    return position.Dot(velocity) / velocity.Magnitude2();
}

bool operator==(const std::reference_wrapper<Particle> &lhs, const Particle &rhs) {
    return lhs.get() == rhs;
}

bool operator==(const Particle &lhs, const std::reference_wrapper<Particle> &rhs) {
    return lhs == rhs.get();
}

void Particle::Move(const ThreeVector& dx) {
	position+=dx;
	for(Particle& part:sub_parts)
		part.Move(dx);
}


std::pair<Particle,Particle> Particle::splitSingle(size_t partIndex) {
	switch(sub_parts.size()) {
		case 0:
			throw std::runtime_error("Particle::splitSingle(): Cannot split a noncomposite particle");
		case 1:
			spdlog::debug("Composite Particle represents only One Piece, returning that piece");
			return std::pair<Particle,Particle>(sub_parts[partIndex],Particle());
		case 2:
			return std::pair<Particle,Particle>(sub_parts[partIndex],sub_parts[1-partIndex]);
	}
	vParticles remaining;
	remaining.reserve(sub_parts.size()-1);
	for(size_t i=0;i<sub_parts.size();i++)
		if(i!=partIndex)
			remaining.push_back(sub_parts[i]);
	return std::pair<Particle,Particle>(sub_parts[partIndex],
		NuclearRemnant(remaining,status,momentum-sub_parts[partIndex].Momentum(),position));
}
std::pair<Particle,Particle> Particle::splitSingleRandom(PID type) {
	std::vector<size_t> candidates;
	for(size_t i=0;i<sub_parts.size();i++)
		if(sub_parts[i].ID()==type)
			candidates.push_back(i);
	if(candidates.empty())
		throw std::runtime_error(fmt::format("NuclearRemnant::splitSingleRandom(): No particle with PID {} was found",type));
	return splitSingle(Random::Instance().Pick(candidates));
}

} // namespace achilles