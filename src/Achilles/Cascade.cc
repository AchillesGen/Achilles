// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include <map>
#include <string>
#include <vector>

#include "Achilles/System.hh"
#include "spdlog/spdlog.h"

#include "Achilles/Cascade.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Exception.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/NuclearMass.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/ResonanceHelper.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

Cascade::Cascade(InteractionHandler interactions, const ProbabilityType &prob, Algorithm alg,
                 const InMedium &medium, std::string_view decay_file, bool potential_prop,
                 double dist)
    : distance(dist), m_interactions(std::move(interactions)),
      m_decays(Filesystem::FindFile(std::string(decay_file), "Cascade")), m_medium(medium),
      m_potential_prop(potential_prop) {
    switch(alg) {
    case Algorithm::Base:
        algorithm = [&](Cascade *cascade,Particle& prop,vParticles& bg) -> size_t {
            return cascade->BaseAlgorithm(prop,bg);
        };
        break;
    case Algorithm::MFP:
        algorithm = [&](Cascade *cascade,Particle& prop,vParticles& bg) -> size_t {
            return cascade->MFPAlgorithm(prop,bg);
        };
        break;
    }

    switch(prob) {
    case ProbabilityType::Gaussian:
        m_probability_name = "Gaussian";
        probability = [](const double &b2, const double &sigma) -> double {
            return exp(-M_PI * b2 / sigma);
        };
        break;
    case ProbabilityType::Pion:
        m_probability_name = "Pion";
        probability = [](const double &b2, const double &sigma) -> double {
            double b = sqrt(b2);
            // TODO: This does not work, we need to rethink this
            // return
            // (135_MeV*sigma)/Constant::HBARC/(2*M_PI*b)*exp(-135_MeV*b/Constant::HBARC);
            return exp(-sqrt(2 * M_PI / sigma) * b);
        };
        break;
    case ProbabilityType::Cylinder:
        m_probability_name = "Cylinder";
        probability = [](const double &b2, const double &sigma) -> double {
            return b2 < sigma / M_PI ? 1 : 0;
        };
    }

    spdlog::debug("Cascade initialized with distance: {}, medium: {}, potential_prop: {}, "
                  "probability: {}, algorithm: {}",
                  distance, ToString(m_medium), m_potential_prop, m_probability_name,
                  ToString(alg));
}

/*void Cascade::Kick(Event &event, const FourVector &energyTransfer,
                   const std::array<double, 2> &sigma) {
    std::vector<std::size_t> indices;

    // Interact with protons or neutrons according to their total cross section
    auto ddSigma = {sigma[0], sigma[1]};
    size_t index = Random::Instance().SelectIndex(ddSigma);

    PID interactPID = index == 0 ? PID::proton() : PID::neutron();

    // Restrict to particles of chosen species
    for(std::size_t i = 0; i < event.Hadrons().size(); ++i) {
        if(event.Hadrons()[i].ID() == interactPID) indices.push_back(i);
    }

    // Kick a single particle from the list
    kickedIdxs.insert(Random::Instance().Pick(indices));
    Particle* kicked = &event.Hadrons()[*kickedIdxs.begin()];
    kicked->Status() = ParticleStatus::propagating;
    kicked->SetMomentum(kicked->Momentum() + energyTransfer);
}*/

std::size_t Cascade::GetInter(vParticles &, const Particle &, double &) {
    throw std::logic_error("Cascade: MFPAlgorithm is currently not working with pions");
    /*
        std::vector<std::size_t> index_same, index_diff;

        for(std::size_t i = 0; i < particles.size(); ++i) {
            if(particles[i].Status() != ParticleStatus::background) continue;
            if(particles[i].ID() == kickedPart.ID())
                index_same.push_back(i);
            else
                index_diff.push_back(i);
        }

        if(index_diff.size() == 0 && index_same.size() == 0) return SIZE_MAX;

        double position = kickedPart.Position().Magnitude();
        auto p1 = kickedPart.Momentum();
        double mass = kickedPart.Info().Mass();

        auto mom = m_nucleus->GenerateMomentum(position);
        double energy = pow(kickedPart.Info().Mass(), 2);
        for(auto p : mom) energy += p * p;
        std::size_t idxSame = SIZE_MAX;
        double xsecSame = 0;
        if(index_same.size() != 0) {
            idxSame = Random::Instance().Pick(index_same);
            particles[idxSame].SetMomentum(FourVector(mom[0], mom[1], mom[2], sqrt(energy)));

            auto p2 = particles[idxSame].Momentum();
            double fact = 1.0;
            if(m_medium == InMedium::NonRelativistic)
                fact = m_nucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position);

            xsecSame = GetXSec(kickedPart, particles[idxSame]) * fact;
        }

        auto otherMass = kickedPart.ID() == PID::proton() ? ParticleInfo(PID::neutron()).Mass()
                                                          : ParticleInfo(PID::proton()).Mass();
        energy = otherMass * otherMass;
        for(auto p : mom) energy += p * p;
        std::size_t idxDiff = SIZE_MAX;
        double xsecDiff = 0;
        if(index_diff.size() != 0) {
            idxDiff = Random::Instance().Pick(index_diff);
            particles[idxDiff].SetMomentum(FourVector(mom[0], mom[1], mom[2], sqrt(energy)));

            auto p2 = particles[idxSame].Momentum();
            double fact = 1.0;
            if(m_medium == InMedium::NonRelativistic)
                fact = m_nucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position);

            xsecDiff = GetXSec(kickedPart, particles[idxDiff]) * fact;
        }

        double rhoSame = 0.0;
        double rhoDiff = 0.0;
        if(position < m_nucleus->Radius()) {
            // TODO: Adjust below to handle non-isosymmetric nuclei
            rhoSame = m_nucleus->Rho(position) * 2 * static_cast<double>(index_same.size()) /
                      static_cast<double>(particles.size());
            rhoDiff = m_nucleus->Rho(position) * 2 * static_cast<double>(index_diff.size()) /
                      static_cast<double>(particles.size());
        }
        if(rhoSame <= 0.0 && rhoDiff <= 0.0) return SIZE_MAX;
        double lambda_tilde = 1.0 / (xsecSame / 10 * rhoSame + xsecDiff / 10 * rhoDiff);
        double lambda = -log(Random::Instance().Uniform(0.0, 1.0)) * lambda_tilde;

        if(lambda > stepDistance) return SIZE_MAX;

        stepDistance = lambda;
        double ichoice = Random::Instance().Uniform(0.0, 1.0);
        if(ichoice < xsecSame / (xsecSame + xsecDiff)) {
            particles[idxSame].SetPosition(kickedPart.Position());
            return idxSame;
        }

        particles[idxDiff].SetPosition(kickedPart.Position());
        return idxDiff;
    */
}

void Cascade::Reset() {
    integrators.clear();
    m_removal_debt.clear();
    m_track_recoil = false;
}

void Cascade::TransferDebt(int from, int to) {
    auto it = m_removal_debt.find(from);
    if(it == m_removal_debt.end()) return;
    m_removal_debt[to] = it->second;
    m_removal_debt.erase(it);
}

bool Cascade::PayRemovalEnergy(Particle& particle) {
	// No remnant = no remaining potential = no need to pay?
	if(!m_remnant) {
		spdlog::debug("Cascade::PayRemovalEnergy(): Cannot pay, no remnant present");
		return true;
	}
    auto it = m_removal_debt.find(particle.cascadeID());
    if(it == m_removal_debt.end()) return true;

    const bool is_proton = it->second.is_proton;

    double separation{};
    if(m_track_recoil) {
        // Charge against the composition settled so far, not the current bound pool, so
        // the charges telescope however the knockouts and escapes interleave.
        separation = SeparationEnergy(m_paid_Z, m_paid_A, is_proton);
    } else {
        int nZ = m_remnant->ID().NuclearZ(), nA = m_remnant->ID().NuclearA();
        separation = SeparationEnergy(nZ + (is_proton ? 1 : 0), nA + 1, is_proton);
    }
    const double base = separation + it->second.kinetic;

    const double mass = particle.Info().Mass();
    const ThreeVector p_old = particle.Momentum().Vec3();
    const ThreeVector dir = p_old.Unit();

    // Recoil: every knockout since the last settlement shoved the remnant by -p_i, and
    // this escape shoves it again along the nucleon's direction. Energy telescoping to
    // the ground-state mass leaves no room for that kinetic energy, so without this term
    // it comes out of the remnant's invariant mass and E_x drifts down ~0.5 MeV per
    // knockout. The escape recoil depends on the debit, so iterate; it converges at once.
    double debit = base;
    if(m_track_recoil) {
        const double rem_mass = m_remnant->Momentum().M();
        const ThreeVector rem_p = m_remnant->Momentum().Vec3();
        const double paid = RecoilKinetic(rem_mass, m_recoil_paid);
        for(int iter = 0; iter < 4; ++iter) {
            const double energy = particle.Momentum().E() - debit;
            if(energy <= mass) return false;
            const double p_new = std::sqrt(energy * energy - mass * mass);
            const ThreeVector rem_after = rem_p + dir * (p_old.Magnitude() - p_new);
            debit = base + RecoilKinetic(rem_mass, rem_after) - paid;
        }
    }

    const double energy = particle.Momentum().E() - debit;
    if(energy <= mass) return false;

    const double momentum = std::sqrt(energy * energy - mass * mass);
    particle.SetMomentum({dir * momentum, energy});
    if(m_track_recoil) {
        m_recoil_paid = m_remnant->Momentum().Vec3() + (p_old - particle.Momentum().Vec3());
        m_paid_Z -= is_proton ? 1 : 0;
        m_paid_A -= 1;
    }
    m_removal_debt.erase(it);
    return true;
}

void Cascade::InitializeIntegrator() {
    for(Particle& part:m_propagating) {
        if(m_potential_prop && m_nucleus->GetPotential()->Hamiltonian(
                                   part.Momentum().P(),
                                   part.Position().P()) < Constant::mN) {
            part.Status() = ParticleStatus::captured;
            m_removal_debt.erase(part.cascadeID());
        } else {
            AddIntegrator(part);
        }
    }
}

void Cascade::UpdateKicked() {
    for(Particle& part:m_propagating) {
        if(m_potential_prop &&
           m_nucleus->GetPotential()->Hamiltonian(
               part.Momentum().P(), part.Position().P()) < Constant::mN) {
            part.Status() = ParticleStatus::captured;
            m_removal_debt.erase(part.cascadeID());
        } else {
            spdlog::trace("Adding particle: {}", part);
            AddIntegrator(part);
        }
    }
}

void Cascade::UpdateTracked(Event& event) {
	m_propagating=event.PropagatingHadrons();
	m_remnant=event.Remnant();
}

void Cascade::Validate(Event &event) {
    spdlog::trace("End of Event Status:");
	if(!m_propagating.empty()) {
		spdlog::error("Propagating Hadrons Remaining");
        for(Particle& p:m_propagating) spdlog::error("{}", p);
        achilles::PrintVisitor visitor;
        event.WalkHistory(visitor);
        spdlog::error("Event history: {}", visitor.data);
        throw AchillesCascadeError("Cascade has failed. Insufficient max steps.");
	}
    for(size_t idx = 0; idx < event.Particles().size(); ++idx) {
        Particle& particle = event.Particles()[idx];
        if(particle.isInFinalState() && particle.Info().IsResonance()
				&& particle.Status() == ParticleStatus::final_state) {
            static constexpr size_t nwarns = 10;
            static size_t iwarns = 0;
            if(iwarns < nwarns) {
                spdlog::warn("Cascade: Resonance did not decay before escaping, decaying now");
            } else if(iwarns == nwarns) {
                spdlog::warn(
                    "Cascade: Reached maximum escaped resonance warnings, suppressing the rest");
            }
            iwarns++;

            // Force decay
            vParticles particles_out = m_decays.Decay(particle);
            particle.Status() = ParticleStatus::decayed;

			size_t decay_vertex=event.addVertex(EventVertex::StatusCode::decay,particle.Position());
			particle.endVertex()=decay_vertex;

            // Ensure outgoing particles are propagating and add to list of particles in event
            // and assign formation zone
            Particle remnant_out;
			if(m_remnant) {
				// If there's no tracked remnant, remnant_out will be ignored later
				remnant_out=*m_remnant;
				m_remnant->SubParts().clear();
			}
            bool repaid = false;
            for(Particle& out : particles_out) {
                out.Status() = ParticleStatus::final_state;
                out.Position() = particle.Position();
                // This path bypasses Escaped, so settle here or the debt is never paid.
                // The particle is already outside the nucleus, so it cannot be recaptured.
                if(out.Info().IsBaryon()) {
					// replaces TransferDebt()
                    out.cascadeID()=particle.cascadeID();
                    const bool owes = m_removal_debt.count(out.cascadeID()) > 0;
                    const FourVector before = out.Momentum();
                    if(!PayRemovalEnergy(out)) {
                        spdlog::debug("Cascade: escaped resonance product cannot repay its "
                                      "removal energy, releasing it unpaid");
                        m_removal_debt.erase(out.cascadeID());
                    } else if(owes && m_remnant != nullptr) {
                        // What the product paid goes to the remnant, as in Escaped.
                        remnant_out.Momentum()+=before-out.Momentum();
                        repaid = true;
                    }
                }
				event.addParticle(out,decay_vertex,Particle::OPEN_END);
            }
			// Event was edited, must account for vector reallocation
			if(m_remnant)
				m_remnant=event.Remnant();

            // When a product repaid removal energy the remnant took it up,
			// so it enters and leaves this vertex to keep it balanced.
            if(repaid&&m_remnant!=nullptr) {
                m_remnant->endVertex()=decay_vertex;
				event.addParticle(remnant_out,decay_vertex,Particle::OPEN_END);
				m_remnant=&event.Particles().back();
            }
        }
    }
}

/*void Cascade::PropagateAll(vParticles &particles, double step) const {
    for(Particle& particle : particles) {
        if(particle.IsPropagating()) { particle.Propagate(step); }
    }
}*/

void Cascade::Evolve(achilles::Event &event, Nucleus *nucleus,
                     [[maybe_unused]] const std::size_t &maxSteps) {
	UpdateTracked(event);
	if(!m_remnant) {
		// If there's no remnant, that *might* be okay because the
		// cascade just doesn't need to run in the first place
		spdlog::debug("Cascade::Evolve(): No remnant particle found in event");
		return;
	}

    // Run the cascade
    currentTime = 0;
    m_nucleus = nucleus;
    // The recoil the remnant already has from the primary vertex is paid for by the
    // spectral function; only recoil the cascade adds from here on is charged.
    m_track_recoil = m_remnant!=nullptr;
    if(m_track_recoil) {
        m_recoil_paid = m_remnant->Momentum().Vec3();
        m_paid_Z = m_remnant->ID().NuclearZ();
        m_paid_A = m_remnant->ID().NuclearA();
    }
    InitializeIntegrator();

	// Commented this block bc m_time_steps is never used -Hayden
	/*vParticles& background=m_remnant->SubParts();
    for(size_t i_prop=0;i_prop<m_propagating.size();i_prop++) {
		Particle& prop=m_propagating[i_prop].get();
        for(size_t i_bg=0;i_bg<background.size();i_bg++) {
			Particle& bg=background[i_bg];
            double closest = ClosestApproach(prop,bg);
            spdlog::debug("Closest approach time({}, {}) = {}", prop, bg, closest);
            if(closest > 0) { m_time_steps.push({closest, {i_prop, i_bg}}); }
        }
    }*/

    // while(!m_time_steps.empty()) {
    //     timeStep = m_time_steps.top().time - currentTime;
    //     auto [kicked, hit] = m_time_steps.top().idxs;
    //     m_time_steps.pop();
    //     PropagateAll(particles, timeStep);
    //     currentTime += timeStep;
    //     if(HasInteraction(particles, kicked, hit)) {
    //         spdlog::trace("Kicked = {}", kicked);
    //         spdlog::trace("Hit = {}", hit);
    //         FinalizeMomentum(event, particles, kicked, hit);
    //         // tmp = m_time_steps;
    //         // pop_println("After Interaction", tmp);
    //         // if(m_time_steps.empty()) throw;
    //     }
    //     Escaped(particles);
    // }

    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Stop loop if no particles are propagating
		UpdateTracked(event);
        if(m_propagating.size() == 0) break;

        // Adapt time step
        AdaptiveStep(m_propagating, distance);

		// Iterate with indices over the whole particles list;
		// avoids concurrent vector modification issues
        for(size_t i=0;i<event.Particles().size();i++) {
			Particle& part=event.Particles()[i];
			// Iterate only over propagating particles
			if(!(part.isInFinalState()&&part.IsPropagating())) continue;

            spdlog::trace("Kicked Particle: {}", part);

            // Attempt decay of unstable particles
            if(!part.Info().IsStable() && Decay(event, part))
				continue;

            // Update formation zones
            if(part.InFormationZone() && !part.Info().IsPion()) {
                Propagate(part);
                continue;
            }

            size_t hitIdx = algorithm(this,part,m_remnant->SubParts());
            if(hitIdx == SIZE_MAX) { continue; }
            spdlog::trace("Kicked = {}", part);
            spdlog::trace("Hit = {}", hitIdx);
            FinalizeMomentum(event,part,hitIdx);
        }
        // Updated the kicked list
        // NOTE: Needs to be here in case there are two interactions in the same time step
        UpdateKicked();

        // After step checks
        Escaped(event);
    }

    Validate(event);
    Reset();
}

size_t Cascade::BaseAlgorithm(Particle& propagating,vParticles& background) {
    // Get allowed interactions
    InteractionDistances dist2 = AllowedInteractions(propagating,background);
    if(dist2.size() == 0) { return SIZE_MAX; }

    // Get interaction
    return Interacted(propagating,dist2);
}

size_t Cascade::MFPAlgorithm(Particle&,vParticles&) {
    throw std::logic_error("Cascade: MFPAlgorithm is currently not working with pions");
    /*double step_prop = distance;
    Particle kickNuc = hadrons[idx];
    size_t hitIdx = GetInter(hadrons, kickNuc, step_prop);
    PropagateSpace(kickNuc, step_prop);
    return hitIdx;*/
}

void Cascade::AddIntegrator(Particle& part) {
    // Only add if not already there
    if(part.cascadeID()!=-1)
		return;

    // Otherwise setup integrator
    static constexpr double omega = 20;
    auto dHamiltonian_dr = [&](const ThreeVector &q, const ThreeVector &p,
                               std::shared_ptr<Potential> potential) {
        auto vals = potential->operator()(p.P(), q.P());
        auto dpot_dp = potential->derivative_p(p.P(), q.P());

        auto mass_eff =
            achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
        double numerator = (vals.rscalar + achilles::Constant::mN) * dpot_dp.rscalar + p.P();
        double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
        return numerator / denominator * p / p.P() + dpot_dp.rvector * p / p.P();
    };
    auto dHamiltonian_dp = [&](const ThreeVector &q, const ThreeVector &p,
                               std::shared_ptr<Potential> potential) {
        auto vals = potential->operator()(p.P(), q.P());
        auto dpot_dp = potential->derivative_p(p.P(), q.P());

        auto mass_eff =
            achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
        double numerator = (vals.rscalar + achilles::Constant::mN) * dpot_dp.rscalar + p.P();
        double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
        return numerator / denominator * p / p.P() + dpot_dp.rvector * p / p.P();
    };

	/// Using this instead of integrators.size() in case integrators can be deleted,
	/// which would cause ID conflicts
	static int last_id=0;
	part.cascadeID()=last_id++;

    integrators[part.cascadeID()] =
        SymplecticIntegrator(part.Position(), part.Momentum().Vec3(), m_nucleus->GetPotential(),
                             dHamiltonian_dr, dHamiltonian_dp, omega);
}

void Cascade::Propagate(Particle& kickNuc) {
    if(m_potential_prop) {
		int idx=kickNuc.cascadeID();
        integrators[idx].Step<2>(timeStep);
        double energy = sqrt(pow(kickNuc.Info().Mass(), 2) + integrators[idx].P().P2());
        FourVector mom{integrators[idx].P(), energy};
        kickNuc.SetMomentum(mom);
        ThreeVector pos_old = kickNuc.Position();
        kickNuc.SetPosition(integrators[idx].Q());
        ThreeVector& pos_new = kickNuc.Position();
        kickNuc.DistanceTraveled() += (pos_new - pos_old).Magnitude();
    } else {
        kickNuc.Propagate(timeStep);
    }
}

void Cascade::PropagateSpace(Particle& kickNuc, double step) {
    double beta = kickNuc.Beta().Magnitude();
    timeStep = step / beta;
    Propagate(kickNuc);
}

// TODO: Switch to using reference_wrapper to reduce copying
void Cascade::MeanFreePath(Event &event, Nucleus *nucleus, const std::size_t &maxSteps) {
    m_nucleus = nucleus;

	m_propagating=event.PropagatingHadrons();
    if(m_propagating.size() != 1) {
        throw std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }
    Particle& kickNuc_original = m_propagating[0].get();
	Particle kickNuc=kickNuc_original;

    // Initialize symplectic integrator
    AddIntegrator(kickNuc_original);
    if(m_potential_prop && m_nucleus->GetPotential()->Hamiltonian(
                               kickNuc.Momentum().P(), kickNuc.Position().P()) < Constant::mN) {
        kickNuc.Status() = ParticleStatus::captured;
        kickNuc_original=kickNuc;
        Reset();
        return;
    }

    if((kickNuc.Status() != ParticleStatus::internal_test) &&
       (kickNuc.Status() != ParticleStatus::external_test)) {
        throw std::runtime_error("MeanFreePath: kickNuc must have status -3 "
                                 "in order to accumulate DistanceTraveled.");
    }
    for(std::size_t step = 0; step < maxSteps; ++step) {
        AdaptiveStep(m_propagating, distance);

        // if(kickNuc->IsUnstable()) {
        // Check if decay
        // }

        if(kickNuc.InFormationZone()) {
            Propagate(kickNuc);
            kickNuc.UpdateFormationZone(timeStep);
            continue;
        }

        // Are we already outside nucleus?
        if(kickNuc.Position().Magnitude() >= nucleus->Radius() &&
           kickNuc.Status() == ParticleStatus::external_test) {
            if(kickNuc.Position().Pz() > nucleus->Radius()) break;

        } else if(kickNuc.Position().Magnitude() >= nucleus->Radius()) {
            kickNuc.Status() = ParticleStatus::final_state;
            break;
        }
        AdaptiveStep(m_propagating, distance);
        // Identify nearby particles which might interact
        InteractionDistances nearby_particles = AllowedInteractions(kickNuc,m_remnant->SubParts());
        if(nearby_particles.size() == 0) continue;

        // Did we hit?
        size_t hitIdx = Interacted(kickNuc, nearby_particles);
        if(hitIdx == SIZE_MAX) continue;

        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        // Stop as soon as we hit anything
        if(FinalizeMomentum(event, kickNuc_original, hitIdx)) break;
    }
    kickNuc_original=kickNuc;
    Reset();
}

// TODO: Switch to using reference_wrapper to reduce copying
void Cascade::MeanFreePath_NuWro(Event &, Nucleus *, const std::size_t &) {
    throw std::logic_error("Cascade: MFPAlgorithm is currently not working with pions");
    /*m_nucleus = nucleus;

	m_propagating=event.PropagatingHadrons();
    if(m_propagating.size() != 1) {
        std::runtime_error("MeanFreePath: only one particle should be kicked.");
    }

    Particle& kickNuc_original = m_propagating[0].get();
	Particle kickNuc=kickNuc_original;

    // Initialize symplectic integrator
    AddIntegrator(kickNuc_original);
    if(m_potential_prop && m_nucleus->GetPotential()->Hamiltonian(
                               kickNuc.Momentum().P(), kickNuc.Position().P()) < Constant::mN) {
        kickNuc.Status() = ParticleStatus::captured;
        kickNuc_original=kickNuc;
        Reset();
        return;
    }

    if(kickNuc.Status() != ParticleStatus::internal_test) {
        std::runtime_error("MeanFreePath: kickNuc must have status -3 "
                           "in order to accumulate DistanceTraveled.");
    }

    for(std::size_t step = 0; step < maxSteps; ++step) {
        // Are we already outside nucleus?
        if(kickNuc.Position().Magnitude() >= nucleus->Radius()) {
            kickNuc.Status() = ParticleStatus::final_state;
            break;
        }
        // AdaptiveStep(particles, distance);
        double step_prop = distance;
        // Did we hit?
        size_t hitIdx = GetInter(particles, kickNuc, step_prop);
        PropagateSpace(kickNuc, step_prop);
        if(hitIdx == SIZE_MAX) continue;
        // Did we *really* hit? Finalize momentum, check for Pauli blocking.
        FinalizeMomentum(event, particles, idx, hitIdx);
        // Stop as soon as we hit anything
        if(particles[idx].Status() == ParticleStatus::interacted) break;
    }

    hadrons = particles;
    Reset();*/
}

// TODO: Rewrite to have the logic built into the Nucleus class
void Cascade::Escaped(Event &event) {
	// The remnant controls whether any energy needs to be paid (and how much).
	// If there's no remnant remaining, all escaping particles are free.
	if(!m_remnant) {
		for(Particle& part:m_propagating)
			part.Status()=ParticleStatus::final_state;
		return;
	}
	// Otherwise, proceed with debt calculation; m_remnant is guaranteed not null
    const double radius2 = pow(m_nucleus->Radius(),2);

	std::vector<size_t> captured_idxs;
	// Cannot iterate over m_propagating, or even Enhanced-For over Particles()
	// because this loop may concurrently modify Particles()
    for(size_t i=0;i<event.Particles().size();i++) {
		Particle& part=event.Particles()[i];
		// Restrict iteraction to *effectively* be over m_propagating
		if(!(part.isInFinalState()&&part.IsPropagating())) continue;

        if(part.Status() == ParticleStatus::external_test) {
            if(pow(part.Position().Z(),2) < radius2) continue;
        }
        if(part.Position().Magnitude2() <= radius2) continue;

        spdlog::debug("Particle Escaping: {}",part);
        // Resonances leave carrying their debt, and are never captured for inability to
        // pay: the remnant tracks nucleons, so absorbing a resonance would corrupt its
        // composition. They settle through their decay products, in Decay or Validate.
        if(!part.Info().IsNucleon()) {
            part.Status() = ParticleStatus::final_state;
            continue;
        }

        // A nucleon must repay the nucleus what it took to get free. If it cannot afford
        // the separation energy plus the Fermi kinetic energy it was given, it stays bound.
        const bool owes = m_removal_debt.count(part.cascadeID()) > 0;
        Particle after=part;
        if(PayRemovalEnergy(after)) {
            // Losing the removal energy at the surface is a physical step, so record it as
            // its own vertex. The remnant takes up the energy and momentum the nucleon gave
            // back, which is what makes the vertex balance.
            if(owes) {
				after.Status() = ParticleStatus::final_state;
				Particle remNew=*m_remnant;
				remNew.Momentum()+=part.Momentum()-after.Momentum();

				size_t escape_vertex=event.addVertex(EventVertex::StatusCode::escape,part.Position());
				part.endVertex()=escape_vertex;
				m_remnant->endVertex()=escape_vertex;
				m_remnant->SubParts().clear();

				event.addParticle(after,escape_vertex,Particle::OPEN_END);
				event.addParticle(remNew,escape_vertex,Particle::OPEN_END);

				// Event has been modified
				UpdateTracked(event);
            } else
				part.Status()=ParticleStatus::final_state;
        } else {
            // It cannot pay, so it rejoins the bound pool and owes nothing after all.
            m_removal_debt.erase(part.cascadeID());
            part.Status()=ParticleStatus::captured;
			captured_idxs.push_back(i);
        }
    }

	// Now, any particles labeled "captured" need to be reabsorbed

	if(captured_idxs.empty())
		return;

	// Create interaction and new remnant particle
	size_t vertex=event.addVertex(EventVertex::StatusCode::capture,m_remnant->Position());
	vParticles remNew_sub=m_remnant->SubParts();
	FourVector p_rem=m_remnant->Momentum();

	// Old remnant is now "dead", can free up memory from composite
	m_remnant->SubParts().clear();
	m_remnant->endVertex()=vertex;

	// Capture all particles and add them to the new remnant
	for(const size_t& idx:captured_idxs) {
		Particle& part=event.Particles()[idx];
		part.endVertex()=vertex;
		remNew_sub.push_back(part);
		p_rem+=part.Momentum();
	}

	// Add new remnant to the event
	event.addParticle(Particle::NuclearRemnant(remNew_sub,m_remnant->Status(),p_rem,m_remnant->Position()),
			vertex,Particle::OPEN_END);

	// Event has now been changed
	UpdateTracked(event);
}

/// Convert a time step in [fm] to [1/MeV].
/// timeStep = distance / max("betas of all kicked particles") / hbarc
void Cascade::AdaptiveStep(const refParticles& propagating, const double &stepDistance) noexcept {
    double beta=0.,beta_max=0.;
    for(const Particle& part:propagating) {
		beta=part.Beta().Magnitude();
        if(beta>beta_max)
			beta_max=beta;
    }

    timeStep = stepDistance / beta_max;
}

/// Determine whether "test_point" is between two parallel planes defined by
/// the normal vector "plane_normal". With the point "on_plane" on the closer plane
/// and dist the distance from the start plane to the end plane.
bool Cascade::BetweenPlanes(const ThreeVector &test_point, const ThreeVector &on_plane,
                            const ThreeVector &plane_normal, double dist) const noexcept {
    // Get signed distance to plane
    double signed_dist = plane_normal.Dot(test_point - on_plane);

    // Determine if point is between planes
    return signed_dist > 0 && signed_dist < dist;
}

/// Project the vector "position" onto the plane containing the point "planePt"
/// and orthogonal to the vector "planeVec" by subtracting off the projection
/// of "position - planePt" onto the normal vector "planeVec".
const ThreeVector Cascade::Project(const ThreeVector &position, const ThreeVector &planePt,
                                   const ThreeVector &planeVec) const noexcept {
    // Project point onto a plane
    const ThreeVector projection = (position - planePt).Dot(planeVec) * planeVec;
    return position - projection;
}

/// Get a sorted list of allowed InteractionDistances, i.e., of pairs
/// (index, distance). Allowed interactions are those between the planes
/// orthogonal to the momentum of the specified particle and located at
/// (i) the particle's current position and
/// (ii) the particle's position after propagating for "timeStep".
/// For example, in the diagram below, interaction with A is allowed, while
/// interactions with B and C are not allowed.
///    plane1  plane2
///      |   A   |
///      |       |
///      x---->  |
///      |       |   B
///  C   |       |
const InteractionDistances Cascade::AllowedInteractions(Particle& prop,const vParticles& background) noexcept {
    InteractionDistances results;

    // Build planes
    const ThreeVector point1 = prop.Position();
    Propagate(prop);
    const ThreeVector point2 = prop.Position();
    ThreeVector normedMomentum = prop.Momentum().Vec3().Unit();
    double distance2 = (point2 - point1).Dot(normedMomentum);

    // Build results vector
    for(std::size_t i = 0; i < background.size(); ++i) {
        // TODO: Should particles propagating be able to interact with
        //       other propagating particles? Currently, no.
        // if(i == idx) continue;
        //  if(particles[i].InFormationZone()) continue;
        if(!BetweenPlanes(background[i].Position(), point1, normedMomentum, distance2)) continue;
        ThreeVector projectedPosition = Project(background[i].Position(), point1, normedMomentum);
        // (Squared) distance in the direction orthogonal to the momentum
        double dist2 = (projectedPosition - point1).Magnitude2();

        results.push_back(std::make_pair(i, dist2));
    }

    // Sort array by distances
    std::sort(results.begin(), results.end(), sortPairSecond);

    return results;
}

// Appears to be dead code -Hayden
/*bool Cascade::HasInteraction(Particle& kicked,Particle& hit) const {
    if(!kicked.IsPropagating() ||
       hit.Status() != ParticleStatus::background || kicked.InFormationZone())
        return false;
    // Get distance between points
    const double distance2 = (kicked.Position() - hit.Position()).Magnitude2();
    const double xsec = GetXSec(kicked, hit);
    spdlog::trace("dist2[{}, {}] = {}", kicked, hit, distance2);
    const double prob = probability(distance2, xsec / 10);
    if(Random::Instance().Uniform(0.0, 1.0) < prob) {
        spdlog::debug("Accepted interaction with: {}, {} => distance2: {}, xsec: {}, prob: {}",
                      kicked, hit, distance2, xsec, prob);
        return true;
    }
    return false;
}*/

double Cascade::GetXSec(Particle& kicked,size_t rem_idx) const {
    // TODO: Clean this up so we don't have to calculate the cross section multiple times
    // TODO: Handle the in-medium effects for different processes
    double fact = InMediumCorrection(kicked,m_remnant->SubParts()[rem_idx]);
    return m_interactions.TotalCrossSection(kicked,m_remnant,rem_idx,m_nucleus) * fact;
}

/// Decide whether or not an interaction occured.
/// The total probability is normalized to the cross section "sigma" when
/// integrated over the plane.
std::size_t Cascade::Interacted(Particle& kicked,
                                const InteractionDistances &dists) noexcept {
    for(auto dist : dists) {
        // Cross section in mb
        const double xsec = GetXSec(kicked, dist.first);
        // 1 barn = 100 fm^2, so 1 mb = 0.1 fm^2.
        // Thus: (xsec [mb]) x (0.1 [fm^2]/ 1 [mb]) = 0.1 xsec [fm^2]
        // dist.second is fm^2; factor of 10 converts mb to fm^2
        const double prob = probability(dist.second, xsec / 10);
        spdlog::trace("dist2[{}, {}] = {}", kicked, dist.first, dist.second);
        if(Random::Instance().Uniform(0.0, 1.0) < prob) return dist.first;
    }

    return SIZE_MAX;
}

bool Cascade::FinalizeMomentum(Event& event,Particle& particle1,size_t rem_idx) noexcept {
	// While building the cascade, particle2 can move between the remnant
	// and the actual history; using a pointer allows it to be easily located
	// again after moving it
	Particle* particle2=&m_remnant->SubParts()[rem_idx];

    // NOTE: No need to deal with in-medium effects since they are just an overall scaling
    InteractionResults modes = m_interactions.CrossSection(particle1,m_remnant,rem_idx,m_nucleus);
    InteractionResult mode = m_interactions.SelectChannel(modes, Random::Instance().Uniform(0.0, 1.0));

    spdlog::debug("Selected Mode: {}, {} -> {}", particle1.ID(), particle2->ID(),
                  fmt::join(mode.particles, ","));

    vParticles particles_out =
        m_interactions.GenerateMomentum(particle1, *particle2, mode.particles, Random::Instance());

    spdlog::debug("Outgoing particles:");
    for(const Particle& part : particles_out) spdlog::debug("- {}", part);

    // Check for Pauli Blocking
    bool hit = true;
    // bool pionIS = false;
    // bool baryonFS = true;

    // Did we start with pions
    // if(particle1.Info().IsPion() || particle2.Info().IsPion()) pionIS = true;

    if(m_nucleus->ID() != PID::hydrogen() && m_nucleus->ID() != PID::free_neutron()) {
        for(const Particle& part : particles_out) {
            hit &= !PauliBlocking(part);
            // Are there any mesons in the final state?
            // if(!part.Info().IsBaryon()) baryonFS = false;
        }
    }

    // If there are pions in in the initial state
    // and no mesons in the final state
    /*if(pionIS && baryonFS) {
        hit = true;
        spdlog::debug("Checking pion abs mom");
        spdlog::debug("Particle1 ({}): {}, Particle2 ({}): {}", particle1.ID(),
                      particle1.Momentum().Vec3().Magnitude(), particle2.ID(),
                      particle2.Momentum().Vec3().Magnitude());
        spdlog::debug("Pion abs FS mom: {}, {}", particles_out[0].Momentum().Vec3().Magnitude(),
                      particles_out[1].Momentum().Vec3().Magnitude());
        if(PauliBlocking(particles_out[0]) || PauliBlocking(particles_out[1])) {
            spdlog::info("Pauliblocked in absorption");
            spdlog::info("outgoing part 1 ({}): {}",
    particles_out[0].ID(),particles_out[0].Momentum().Vec3().Magnitude()); spdlog::info("part 1
    Fermi momentum = {}",m_nucleus->FermiMomentum(particles_out[0].Position().Magnitude()));
            spdlog::info("outgoing part 2 ({}): {}",
    particles_out[1].ID(),particles_out[1].Momentum().Vec3().Magnitude()); spdlog::info("part 2
    Fermi momentum = {}",m_nucleus->FermiMomentum(particles_out[1].Position().Magnitude()));
        }
    }*/

    if(!hit) {
		// Reject interaction and return if Pauli blocked
	    for(Particle& part : m_remnant->SubParts()) {
	        if(part.Status() == ParticleStatus::absorption_partner)
	            part.Status() = ParticleStatus::background;
	    }
		return false;
	}

    // Collect the removal debts this vertex passes on, before any status is rewritten below.
    // A participant already out of the bound pool forwards what it owed; one still in the
    // background opens a new debt for the nucleon the nucleus is about to lose.
    std::vector<RemovalDebt> debts;
	// "consumed" represents indices within the remnant's sublist,
	// of particles which need to be freed so they can interact.
    std::set<size_t> consumed;
    auto forward = [&](int id) {
        if(auto it = m_removal_debt.find(id); it != m_removal_debt.end()) {
            debts.push_back(it->second);
            m_removal_debt.erase(it);
        }
    };
    auto open = [&](size_t idx) {
		Particle& part=m_remnant->SubParts()[idx];
        debts.push_back({part.Momentum().E() - part.Info().Mass(), part.ID() == PID::proton()});
        consumed.insert(idx);
    };
    forward(particle1.cascadeID());
    forward(particle2->cascadeID());
	// Updating from previous if statements:
	// particle1 is never background, particle2 is always background
    open(rem_idx);
    for(size_t i=0;i<m_remnant->SubParts().size();i++) {
        if(m_remnant->SubParts()[i].Status() == ParticleStatus::absorption_partner) open(i);
    }

    // Nucleons the cascade pulls out of the bound pool leave the remnant here, so the
    // vertex balances: projectile + remnant_in -> products + remnant_out. Absorption
    // partners are consumed the same way, which is what makes those vertices conserve.
    // `consumed` was filled above, before the statuses below were rewritten.
    //const bool threaded = static_cast<size_t>(m_remnant->ID().NuclearA()) > consumed.size();
    //if(threaded) {

	FourVector p_rem=m_remnant->Momentum();
	vParticles remNew_sub;
	remNew_sub.reserve(m_remnant->SubParts().size()-consumed.size());

	// List of indices within event.Particles().
	// Must use indices instead of refs in case of vector resizing.
	std::vector<size_t> freed;
	freed.reserve(consumed.size());

	// Interaction for separating all consumed nucleons from the remnant
	size_t vertex_split_rem=event.addVertex(EventVertex::StatusCode::split,m_remnant->Position());

	for(size_t i=0;i<m_remnant->SubParts().size();i++) {
		auto it=std::find(consumed.begin(),consumed.end(),i);
		if(it==consumed.end()) {
			// Particle stays in background => move to new remnant
			remNew_sub.push_back(m_remnant->SubParts()[i]);
		} else {
			// Particle is interactable => move to final state
			if(event.addParticle(m_remnant->SubParts()[i],vertex_split_rem,Particle::OPEN_END)) {
				// Reallocation, must find remnant again
				m_remnant=event.Remnant();
			}
			if(i==rem_idx) {
				// We've moved the particle we wanted to track
				particle2=&event.Particles().back();
			}
			p_rem-=m_remnant->SubParts()[i].Momentum();
			// Last particle in event is now freed
			freed.push_back(event.Particles().size()-1);
		}
	}
	m_remnant->endVertex()=vertex_split_rem;
	// Sublist in previous remnant is now "dead"
	m_remnant->SubParts().clear();
	if(!remNew_sub.empty()) {
		// Create new remnant containing remaining background nucleons
		event.addParticle(Particle::NuclearRemnant(remNew_sub,m_remnant->Status(),p_rem,m_remnant->Position()),
			vertex_split_rem,Particle::OPEN_END);
		// Particle that was just added is the new remnant
		m_remnant=&event.Particles().back();
	} else {
		// Remnant has been fully dissolved
		m_remnant=nullptr;
	}
    //}

	// The actual hit is now treated as a separate interaction -Hayden
    // TODO: What do we use for the position? (How about average positions?)
    ThreeVector average_position = (particle1.Position() + particle2->Position()) / 2.0;
	size_t vertex_hit=event.addVertex(EventVertex::StatusCode::cascade,average_position);

    particle1.Status() = ParticleStatus::interacted;
    particle1.endVertex() = vertex_hit;
    particle2->Status() = ParticleStatus::interacted;
    particle2->endVertex() = vertex_hit;
    for(const size_t& i:freed) {
		Particle& part=event.Particles()[i];
        if(part.Status() == ParticleStatus::absorption_partner)
            part.endVertex()=vertex_hit;
    }

    // Ensure outgoing particles are propagating and add to list of particles in event
    // and assign formation zone
    size_t idebt = 0;
    for(Particle& part : particles_out) {
        part.Status() = ParticleStatus::propagating;
        if(part.Info().IsNucleon())
            part.SetFormationZone(particle1.Momentum(), part.Momentum());
        event.addParticle(part,vertex_hit,Particle::OPEN_END);
        // Baryon number is conserved, so there is always a baryon for every debt.
        if(part.Info().IsBaryon() && idebt < debts.size())
            m_removal_debt[part.cascadeID()] = debts[idebt++];

        // TODO: Work out the queue to minimize number of calls to AllowedInteractions
        // Update queue of closest approach times
        // for(size_t i = 0; i < particles.size(); ++i) {
        //     if(particles[i].Status() != ParticleStatus::background) continue;
        //     double closest = ClosestApproach(part, particles[i]);
        //     spdlog::trace("Closest approach time({}) = {}", i, closest+currentTime);
        //     if(closest > 0) { m_time_steps.push({closest+currentTime, {particles.size()-1,
        //     i}}); }
        // }
    }

	// Event has been modified
	UpdateTracked(event);
	return true;
}

// TODO: Rewrite to have most of the logic built into the Nucleus class?
bool Cascade::PauliBlocking(const Particle &particle) const noexcept {
    if(!particle.Info().IsNucleon()) return false;
    double position = particle.Position().Magnitude();
    return particle.Momentum().Vec3().Magnitude() <
           m_nucleus->FermiMomentum(position, particle.ID());
}

double Cascade::InMediumCorrection(const Particle &particle1, const Particle &particle2) const {
    if(m_medium != InMedium::NonRelativistic) return 1;

    if(!particle1.Info().IsNucleon() || !particle2.Info().IsNucleon()) return 1;

    const FourVector& p1 = particle1.Momentum();
    const FourVector& p2 = particle2.Momentum();
    double mass = particle1.Info().Mass();
    const ThreeVector& pos_p1 = particle1.Position();
    const ThreeVector& pos_p2 = particle2.Position();
    double position1 = pos_p1.Magnitude();
    double position2 = pos_p2.Magnitude();
    double position3 = (pos_p1 + pos_p2).Magnitude();
    return m_nucleus->GetPotential()->InMediumCorrectionNonRel(p1, p2, mass, position1, position2,
                                                               position3);
}

bool Cascade::Decay(Event& event,Particle& part) {
    double beta = part.Beta().Magnitude();
    double gamma = 1. / sqrt(1. - beta * beta);
    double lifetime = gamma * Constant::HBARC / part.Info().Width();
    double survival_prob = exp(-timeStep / lifetime);

    if(part.Info().IsDelta()) {
        const double mn =
            (ParticleInfo(PID::proton()).Mass() + ParticleInfo(PID::neutron()).Mass()) / 2 / 1_GeV;
        const double mpi =
            (2 * ParticleInfo(PID::pionp()).Mass() + ParticleInfo(PID::pion0()).Mass()) / 3 / 1_GeV;

        double running_width = resonance::GetEffectiveWidth(part.Info().ID(),
                                                          part.Momentum().M() / 1_GeV, mpi, mn, 1) *
                             1_GeV;
        double running_lifetime = gamma * Constant::HBARC / running_width;

        spdlog::debug("Delta minv = {}, mass = {}, fixed width = {}", part.Momentum().M(),
                      part.Info().Mass(), part.Info().Width());
        spdlog::debug("Running width = {}", running_width);
        spdlog::debug("Running survival prob = {}", exp(-timeStep / running_lifetime));
        spdlog::debug("Fixed survival prob = {}", exp(-timeStep / lifetime));

        survival_prob = exp(-timeStep / running_lifetime);
    }

    // Should we attempt a decay in this time step
    spdlog::debug("survival prob = {}, timestep = {}, lifetime = {}", survival_prob, timeStep,
                  lifetime);
    if(Random::Instance().Uniform(0.0, 1.0) < survival_prob) return false;

    // Look up in decay handler
    vParticles particles_out = m_decays.Decay(part);

    for(const Particle& out : particles_out) {
        if(PauliBlocking(out)) return false;
        if(std::isnan(out.Momentum()[0])) {
            spdlog::error("Nan momentum in decay");
            spdlog::error("Pin = {}, Pout = [{}, {}]", part, particles_out[0], particles_out[1]);
            throw AchillesCascadeError("Nan momentum in decay");
        }
    }

	// From this point on, assume the decay IS happening
	// (so editing the event is now allowed)
	size_t decay_vertex=event.addVertex(EventVertex::StatusCode::decay,part.Position());

    part.Status() = ParticleStatus::decayed;
	part.endVertex()=decay_vertex;

    // Add decay results to the event as new propagating hadrons
    std::vector<Particle> final;
	for(Particle& out : particles_out) {
		out.Status() = ParticleStatus::propagating;
		out.SetFormationZone(out.Momentum(), part.Momentum());
		// The baryon inherits whatever the resonance owed; mesons owe nothing.
		// Edit by Hayden: This simply means that the integrator/debt ID is
		// now carried by the decayed baryon.
		if(out.Info().IsBaryon()) {
			out.cascadeID()=part.cascadeID();
		}
		if(event.addParticle(out,decay_vertex,Particle::OPEN_END)) {
			UpdateTracked(event);
		}
	}

    return true;
}