#include <random>
#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "spdlog/spdlog.h"

#include "nuchic/Constants.hh"
#include "nuchic/Cascade.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Utilities.hh"
#include "nuchic/Interactions.hh"

using namespace nuchic;

Cascade::Cascade(const std::shared_ptr<Interactions> interactions,
                 const ProbabilityType& prob,
                 const double& dist = 0.05) 
        : distance(dist), m_interactions(interactions) {

    switch(prob) {
        case ProbabilityType::Gaussian:
            probability = [](const double &b2, const double &sigma) -> double {
                return exp(-M_PI*b2/sigma);
            };
            break;
        case ProbabilityType::Pion:
            probability = [](const double &b2, const double &sigma) -> double {
                double b = sqrt(b2);
                //return (135_MeV*sigma)/Constant::HBARC/(2*M_PI*b)*exp(-135_MeV*b/Constant::HBARC); 
		return exp(-sqrt(2*M_PI/sigma)*b);
            };
            break;
    }
}

void Cascade::Kick(std::shared_ptr<Nucleus> nucleus, const FourVector& energyTransfer,
                   const std::array<double, 2>& sigma) {
    std::vector<std::size_t> indices;

    auto ddSigma = {sigma[0], sigma[1]};
    auto index = rng.variate<std::size_t, std::discrete_distribution>(ddSigma);

    auto interactPID = index == 0 ? PID::proton() : PID::neutron();

    for(std::size_t i = 0; i < nucleus -> Nucleons().size(); ++i) {
        if(nucleus -> Nucleons()[i].ID() == interactPID) indices.push_back(i);
    }

    kickedIdxs.push_back(rng.pick(indices));
    auto kicked = &nucleus -> Nucleons()[kickedIdxs.back()];
    kicked -> SetStatus(ParticleStatus::propagating);
    kicked -> SetMomentum(kicked -> Momentum() + energyTransfer);
}

std::size_t  Cascade::GetInter(Particles& particles, const Particle& particle1, double& stepDistance) {
    std::vector<std::size_t> index_same, index_diff;

    for(std::size_t i = 0; i < particles.size(); ++i) {
        if (particles[i].Status() != ParticleStatus::background) continue;
	if (particles[i].ID() == particle1.ID()) index_same.push_back(i);
        else index_diff.push_back(i);
    }

    if(index_diff.size()==0 && index_same.size()==0) return SIZE_MAX; 
    double position = particle1.Position().Magnitude();	
    auto mom = localNucleus -> GenerateMomentum(position);
    double energy = Constant::mN*Constant::mN;
    for(auto p : mom) energy += p*p;

    std::size_t idx1 = SIZE_MAX;
    double xsec1 = 0;
    if(index_same.size() != 0) {
        idx1 = rng.pick(index_same);
        particles[idx1].SetMomentum(FourVector(mom[0], mom[1], mom[2], sqrt(energy)));
        xsec1 = GetXSec(particle1, particles[idx1]);
    }

    std::size_t idx2 = SIZE_MAX;
    double xsec2 = 0;
    if(index_diff.size() != 0) {
        idx2 = rng.pick(index_diff);
        particles[idx2].SetMomentum(FourVector(mom[0], mom[1], mom[2], sqrt(energy)));
        xsec2 = GetXSec(particle1, particles[idx2]);
    }
 
    double rho = localNucleus -> Rho(position);
    if(rho <= 0.0) return SIZE_MAX;
    double lambda_tilde = 1.0/(xsec1/10*rho+xsec2/10*rho);
    double lambda = -log(rng.uniform(0.0, 1.0))*lambda_tilde;
    
    if(lambda > stepDistance) return SIZE_MAX;
    
    stepDistance = lambda;
    double ichoice = rng.uniform(0.0, 1.0);
    if(ichoice < xsec1/(xsec1+xsec2)) {
       particles[idx1].SetPosition(particle1.Position());	    
       return idx1;
    } 

    particles[idx2].SetPosition(particle1.Position());
    return idx2;
}    
	
void Cascade::Reset() {
    kickedIdxs.resize(0);
}

void Cascade::Evolve(std::shared_ptr<Nucleus> nucleus, const std::size_t& maxSteps = 100000) {
    localNucleus = nucleus;
    Particles particles = nucleus -> Nucleons();

    for(std::size_t step = 0; step < maxSteps; ++step) {
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
                kickNuc -> Propagate(timeStep);
                continue;
            }

            // Get allowed interactions
            auto dist2 = AllowedInteractions(particles, idx);
            if(dist2.size() == 0) continue;

            // Get interaction
            auto hitIdx = Interacted(particles, *kickNuc, dist2);
            if(hitIdx == SIZE_MAX) continue;
            Particle* hitNuc = &particles[hitIdx];

            // Finalize Momentum
            bool hit = FinalizeMomentum(*kickNuc, *hitNuc);
            if(hit) {
                newKicked.push_back(hitIdx);
                hitNuc->SetStatus(ParticleStatus::propagating);
            }
        }

        // Replace kicked indices with new list
        kickedIdxs = newKicked;

        // After step checks
        Escaped(particles);
    }

    for(auto particle : particles) {
        if(particle.Status() == ParticleStatus::propagating) {
            std::cout << "\n";
            for(auto p : particles) spdlog::error("{}", p);
            throw std::runtime_error("Cascade has failed. Insufficient max steps.");
        }
    }

    nucleus -> Nucleons() = particles;
}

void Cascade::NuWro(std::shared_ptr<Nucleus> nucleus, const std::size_t& maxSteps = 5000000) { 
    localNucleus = nucleus;
    Particles particles = nucleus -> Nucleons();

    for(std::size_t step = 0; step < maxSteps; ++step) {
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
	//	timeStep= distance/(kickNuc -> Beta().Magnitude()*Constant::HBARC);			   
                kickNuc -> UpdateFormationZone(timeStep);
                kickNuc -> SpacePropagate(distance);		
                continue;
            }

	    double step_prop = distance;
            auto hitIdx = GetInter(particles,*kickNuc,step_prop); 
            kickNuc -> SpacePropagate(step_prop);
            if (hitIdx == SIZE_MAX) continue;
	    Particle* hitNuc = &particles[hitIdx];
            bool hit = FinalizeMomentum(*kickNuc, *hitNuc);
            if(hit) {
                newKicked.push_back(hitIdx);
                hitNuc->SetStatus(ParticleStatus::propagating);
            }
        }

        // Replace kicked indices with new list
        kickedIdxs = newKicked;

        Escaped(particles);
    }

    for(auto particle : particles) {
        if(particle.Status() == ParticleStatus::propagating) {
            for(auto p : particles) std::cout << p << std::endl;
            throw std::runtime_error("Cascade has failed. Insufficient max steps.");
        }
    }

    nucleus -> Nucleons() = particles;
}

void Cascade::MeanFreePath(std::shared_ptr<Nucleus> nucleus, const std::size_t& maxSteps = 10000) {
    localNucleus = nucleus;
    Particles particles = nucleus -> Nucleons();

    auto idx = kickedIdxs[0];
    Particle* kickNuc = &particles[idx];
    
    if (kickedIdxs.size() != 1) {
        std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }
    if (kickNuc -> Status() != ParticleStatus::internal_test) {
        std::runtime_error(
            "MeanFreePath: kickNuc must have status -3 "
            "in order to accumulate DistanceTraveled."
            );
    }
    bool hit = false;
    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Are we already outside nucleus?
        if (kickNuc -> Position().Magnitude() >= nucleus -> Radius()) {
            kickNuc -> SetStatus(ParticleStatus::escaped);
            break;
        }
        AdaptiveStep(particles, distance);
        // Identify nearby particles which might interact
        auto nearby_particles = AllowedInteractions(particles, idx);
        if (nearby_particles.size() == 0) continue;
        // Did we hit?
        auto hitIdx = Interacted(particles, *kickNuc, nearby_particles);
        if (hitIdx == SIZE_MAX) continue;
        Particle* hitNuc = &particles[hitIdx];
        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        hit = FinalizeMomentum(*kickNuc, *hitNuc);
        // Stop as soon as we hit anything
        if (hit) break;
    }

    nucleus -> Nucleons() = particles;
    Reset();
}


void Cascade::MeanFreePath_NuWro(std::shared_ptr<Nucleus> nucleus, const std::size_t& maxSteps = 10000) {
    localNucleus = nucleus;
    Particles particles = nucleus -> Nucleons();

    auto idx = kickedIdxs[0];
    Particle* kickNuc = &particles[idx];
    
    if (kickedIdxs.size() != 1) {
        std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }
    if (kickNuc -> Status() != ParticleStatus::internal_test) {
        std::runtime_error(
            "MeanFreePath: kickNuc must have status -3 "
            "in order to accumulate DistanceTraveled."
            );
    }
    bool hit = false;
    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Are we already outside nucleus?
        if (kickNuc -> Position().Magnitude() >= nucleus -> Radius()) {
            kickNuc -> SetStatus(ParticleStatus::escaped);
            break;
        }
        //AdaptiveStep(particles, distance);
        double step_prop = distance;
        // Did we hit?
        auto hitIdx = GetInter(particles,*kickNuc,step_prop); 
        kickNuc -> SpacePropagate(step_prop);	
       	if (hitIdx == SIZE_MAX) continue;
        Particle* hitNuc = &particles[hitIdx];
        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        hit = FinalizeMomentum(*kickNuc, *hitNuc);
        // Stop as soon as we hit anything
        if (hit) break;
    }

    nucleus -> Nucleons() = particles;
    Reset();
}

// TODO: Rewrite to have most of the logic built into the Nucleus class?
void Cascade::Escaped(Particles &particles) {
    //std::cout << "CURRENT STEP: " << step << std::endl;
    for(auto it = kickedIdxs.begin() ; it != kickedIdxs.end(); ) {
        // Nucleon outside nucleus (will check if captured or escaped after cascade)
        auto particle = &particles[*it];
        if(particle -> Status() == ParticleStatus::background)
            throw std::domain_error("Invalid Particle in kicked list");
        //if(particle -> Status() == -2) {
        //    std::cout << particle -> Position().Pz() << " " << sqrt(radius2) << std::endl;
        //    std::cout << *particle << std::endl;
        //}
        auto radius = localNucleus -> Radius();
        if(particle -> Position().Magnitude2() > pow(radius, 2)
           && particle -> Status() != ParticleStatus::external_test) {
            particle -> SetStatus(ParticleStatus::escaped);
            it = kickedIdxs.erase(it);
        } else if(particle -> Status() == ParticleStatus::external_test && particle -> Position().Pz() > sqrt(radius)) {
            it = kickedIdxs.erase(it);
        } else {
            ++it;
        }
    }
}

void Cascade::AdaptiveStep(const Particles& particles, const double& stepDistance) noexcept {
    double beta = 0;
    for(auto idx : kickedIdxs) {
        if(particles[idx].Beta().Magnitude() > beta)
            beta = particles[idx].Beta().Magnitude();
    }

    timeStep = stepDistance/(beta*Constant::HBARC);
}

bool Cascade::BetweenPlanes(const ThreeVector& position,
                            const ThreeVector& point1,
                            const ThreeVector& point2) const noexcept {
    // Get distance between planes
    const ThreeVector dist = point2 - point1;

    // Determine if point is between planes
    return ((position - point1).Dot(dist) >= 0 && (position - point2).Dot(dist) <=0);
}

const ThreeVector Cascade::Project(const ThreeVector& position,
                                   const ThreeVector& planePt,
                                   const ThreeVector& planeVec) const noexcept {
    // Project point onto a plane
    const ThreeVector projection = (position - planePt).Dot(planeVec)*planeVec; 
    return position - projection;
}

const InteractionDistances Cascade::AllowedInteractions(Particles& particles,
                                                        const std::size_t& idx) const noexcept {
    InteractionDistances results;

    // Build planes
    const ThreeVector point1 = particles[idx].Position();
    particles[idx].Propagate(timeStep);
    const ThreeVector point2 = particles[idx].Position();
    auto normedMomentum = particles[idx].Momentum().Vec3().Unit();

    // Build results vector
    for(std::size_t i = 0; i < particles.size(); ++i) {
        // TODO: Should particles propagating be able to interact with 
        //       other propagating particles?
        if (particles[i].Status() < ParticleStatus::background) continue;
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

double Cascade::GetXSec(const Particle& particle1, const Particle& particle2) const {
    return m_interactions -> CrossSection(particle1, particle2);
}

std::size_t Cascade::Interacted(const Particles& particles, const Particle& kickedParticle,
                                const InteractionDistances& dists) noexcept {
    for(auto dist : dists) {
        const double xsec = GetXSec(kickedParticle, particles[dist.first]);
        const double prob = probability(dist.second, xsec/10);
        if(rng.uniform(0.0, 1.0) < prob) return dist.first;
    }

    return SIZE_MAX;
}

bool Cascade::FinalizeMomentum(Particle& particle1, Particle& particle2) noexcept {
    // Boost to center of mass frame
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM), p2CM = p2Lab.Boost(-boostCM);

    // Generate outgoing momentum
    bool samePID = particle1.ID() == particle2.ID(); 
    double ecm = (p1CM + p2CM).M();
    const double pcm = particle1.Momentum().Vec3().Magnitude() * Constant::mN / ecm;
    std::array<double, 2> rans{};
    rng.generate(rans, 0.0, 1.0);
    ThreeVector momentum = m_interactions -> MakeMomentum(samePID,
                                                          p1CM.Vec3().Magnitude(),
                                                          pcm,
                                                          rans);

    FourVector p1Out = FourVector(momentum[0], momentum[1], momentum[2], p1CM.E()); 
    FourVector p2Out = FourVector(-momentum[0], -momentum[1], -momentum[2], p1CM.E()); 
   
    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    // Assign momenta to particles
    particle1.SetMomentum(p1Out);
    particle2.SetMomentum(p2Out);

    // Check for Pauli Blocking
    bool hit = !(PauliBlocking(particle1) || PauliBlocking(particle2));

    if(hit) {
        // Assign formation zone
        particle1.SetFormationZone(p1Lab, p1Out);
        particle2.SetFormationZone(p2Lab, p2Out);

        // Hit nucleon is now propagating
        // Users are responsibile for updating the status externally as desired
    } else {
        // Assign momenta to particles
        particle1.SetMomentum(p1Lab);
        particle2.SetMomentum(p2Lab);
    }

    return hit;
}

// TODO: Rewrite to have most of the logic built into the Nucleus class?
bool Cascade::PauliBlocking(const Particle& particle) const noexcept { 
    double position = particle.Position().Magnitude();	
    return particle.Momentum().Vec3().Magnitude() < localNucleus -> FermiMomentum(position);
}
