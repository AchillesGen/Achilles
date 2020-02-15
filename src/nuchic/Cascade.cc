#include <random>
#include <iostream>

#include "spdlog/spdlog.h"

#include "nuchic/Cascade.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Utilities.hh"
#include "nuchic/Interactions.hh"

#define HBARC 200
#define mN 938

nuchic::Particles nuchic::Cascade::Kick(const nuchic::Particles& particles, const nuchic::FourVector& energyTransfer,
        const std::array<double, 2>& sigma) {
    std::vector<std::size_t> indices;

    auto ddSigma = {sigma[0], sigma[1]};
    std::size_t index = rng.variate<int, std::discrete_distribution>(ddSigma);

    int interactPID = index == 0 ? 2212 : 2112;

    for(std::size_t i = 0; i < particles.size(); ++i) {
        if(particles[i].PID() == interactPID) indices.push_back(i);
    }

    kickedIdxs.push_back(rng.pick(indices));
    Particles result = particles;
    result[kickedIdxs.back()].SetStatus(-1);
    result[kickedIdxs.back()].SetMomentum(result[kickedIdxs.back()].Momentum() + energyTransfer);

    return result;
}

void nuchic::Cascade::Reset() {
    kickedIdxs.resize(0);
}

nuchic::Particles nuchic::Cascade::operator()(const nuchic::Particles& _particles, const double& kf, const double& radius2,
        const std::size_t& maxSteps) {

    nuchic::Particles particles = _particles;
    fermiMomentum = kf;
    for(std::size_t step = 0; step < maxSteps; ++step) {
        SPDLOG_DEBUG("Step number: %d" step);
        // Stop loop if no particles are propagating
        if(kickedIdxs.size() == 0) break;
        
        // Adapt time step
        AdaptiveStep(particles, distance);

        // Make local copy of kickedIdxs
        auto newKicked = kickedIdxs;
        for(auto idx : kickedIdxs) {
            Particle* kickNuc = &particles[idx];
            // Update formation zones
            if(kickNuc -> InFormationZone()) {
                kickNuc -> UpdateFormationZone(timeStep);
                continue;
            }

            // Get allowed interactions
            auto dist2 = AllowedInteractions(particles, idx);
            if(dist2.size() == 0) continue;

            // Get interaction
            auto hitIdx = Interacted(particles, *kickNuc, dist2);
            if(hitIdx == -1) continue;
            Particle* hitNuc = &particles[hitIdx];

            // Finalize Momentum
            bool hit = FinalizeMomentum(*kickNuc, *hitNuc);

            if(hit) newKicked.push_back(hitIdx);
        }

        // Replace kicked indices with new list
        kickedIdxs = newKicked;

        // After step checks
        //std::cout << "CURRENT STEP: " << step << std::endl;
        for(auto it = kickedIdxs.begin() ; it != kickedIdxs.end(); ) {
            // Nucleon outside nucleus (will check if captured or escaped after cascade)
            auto particle = &particles[*it];
            //if(particle -> Status() == -2) {
            //    std::cout << particle -> Position().Pz() << " " << sqrt(radius2) << std::endl;
            //    std::cout << *particle << std::endl;
            //}
            if(particle -> Position().Magnitude2() > radius2 && particle -> Status() != -2) {
                particle -> SetStatus(1);
                it = kickedIdxs.erase(it);
            } else if(particle -> Status() == -2 && particle -> Position().Pz() > sqrt(radius2)) {
                it = kickedIdxs.erase(it);
            } else {
                ++it;
            }
        }

        for(auto particle : particles) {
            if(particle.Position()[0] != particle.Position()[0])
                throw;
        }
    }

    for(auto particle : particles) {
        if(particle.Status() == -1) {
            for(auto p : particles) std::cout << p << std::endl;
            throw std::runtime_error("Cascade has failed. Insufficient max steps.");
        }
    }

    return particles;
}

void nuchic::Cascade::AdaptiveStep(const Particles& particles, const double& distance) noexcept {
    double beta = 0;
    for(auto idx : kickedIdxs) {
        if(particles[idx].Beta().Magnitude() > beta)
            beta = particles[idx].Beta().Magnitude();
    }

    timeStep = distance/(beta*HBARC);
}

bool nuchic::Cascade::BetweenPlanes(const nuchic::ThreeVector& position,
                                    const nuchic::ThreeVector& point1,
                                    const nuchic::ThreeVector& point2) const noexcept {
    // Get distance between planes
    const nuchic::ThreeVector dist = point2 - point1;

    // Determine if point is between planes
    return ((position - point1).Dot(dist) >= 0 && (position - point2).Dot(dist) <=0);
}

const nuchic::ThreeVector nuchic::Cascade::Project(const nuchic::ThreeVector& position,
                                                   const nuchic::ThreeVector& planePt,
                                                   const nuchic::ThreeVector& planeVec) const noexcept {
    // Project point onto a plane
    const nuchic::ThreeVector projection = (position - planePt).Dot(planeVec)*planeVec; 
    return position - projection;
}

const nuchic::InteractionDistances nuchic::Cascade::AllowedInteractions(nuchic::Particles& particles,
                                                                        const std::size_t& idx) const noexcept {
    nuchic::InteractionDistances results;

    // Build planes
    const nuchic::ThreeVector point1 = particles[idx].Position();
    particles[idx].Propagate(timeStep);
    const nuchic::ThreeVector point2 = particles[idx].Position();
    auto normedMomentum = particles[idx].Momentum().Vec3().Unit();

    // Build results vector
    for(std::size_t i = 0; i < particles.size(); ++i) {
        // TODO: Should particles propagating be able to interact with 
        //       other propagating particles?
        if(particles[i].Status() == -1 || particles[i].Status() == -2) continue;
        // if(i == idx) continue;
        // if(particles[i].InFormationZone()) continue;
        if(!BetweenPlanes(particles[i].Position(), point1, point2)) continue;
        auto projectedPosition = Project(particles[i].Position(), point1, normedMomentum);
        double dist2 = (projectedPosition - point1).Magnitude2();

        results.push_back(std::make_pair(i, dist2));
    }

    // Sort array by distances
    std::sort(results.begin(), results.end(), sortPairSecond);

    return results;
}

const double nuchic::Cascade::GetXSec(const nuchic::Particle& particle1,
                                      const nuchic::Particle& particle2) const {
    return m_interactions -> CrossSection(particle1, particle2);
}

int nuchic::Cascade::Interacted(const Particles& particles, const nuchic::Particle& kickedParticle,
                                const InteractionDistances& dists) noexcept {
    for(auto dist : dists) {
        const double xsec = GetXSec(kickedParticle, particles[dist.first]);
        const double prob = 1.0/(2.0*M_PI)*exp(-dist.second/(2*xsec/10.));
        if(rng.uniform(0.0, 1.0) < prob) return dist.first;
    }

    return -1;
}

bool nuchic::Cascade::FinalizeMomentum(nuchic::Particle& particle1,
                                       nuchic::Particle& particle2) noexcept {
    // Boost to center of mass frame
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM), p2CM = p2Lab.Boost(-boostCM);

    // Generate outgoing momentum
    bool samePID = particle1.PID() == particle2.PID(); 
    double ecm = (p1CM + p2CM).E();
    const double pcm = particle1.Momentum().Vec3().Magnitude() * mN / ecm;
    std::array<double, 2> rans;
    rng.generate(rans, 0.0, 1.0);
    ThreeVector momentum = m_interactions -> MakeMomentum(samePID, p1CM.Vec3().Magnitude(), pcm, rans);

    FourVector p1Out = FourVector(momentum[0], momentum[1], momentum[2], p1CM.E()); 
    FourVector p2Out = FourVector(-momentum[0], -momentum[1], -momentum[2], p1CM.E()); 
   
    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    // Check for Pauli Blocking
    bool hit = !(PauliBlocking(p1Out) || PauliBlocking(p2Out));

    if(hit) {
        // Assign formation zone
        particle1.SetFormationZone(particle1.Momentum(), p1Out);
        particle2.SetFormationZone(particle2.Momentum(), p2Out);

        // Assign momenta to particles
        particle1.SetMomentum(p1Out);
        particle2.SetMomentum(p2Out);

        // Hit nucleon is now propagating
        particle2.SetStatus(-1);
    }

    return hit;
}

bool nuchic::Cascade::PauliBlocking(const nuchic::FourVector& momentum) const noexcept {
    return momentum.Vec3().Magnitude() < fermiMomentum;
}
