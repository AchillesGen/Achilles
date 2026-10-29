// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Event.hh"

using achilles::BeamMapper;

void BeamMapper::GeneratePoint(Event& event, const std::vector<double> &rans) {
    // TODO: Should Masses().back() be the mass of the final state hadronic system or the initial?
	const PID& pid=event.processInfo().m_leptonic.first;
    FourVector p_init_lep=m_beam->Flux(*m_beam->BeamIDs().begin(), rans, (Smin() - Masses()[1]) / (2 * sqrt(Masses()[1])));

	Particle init_lep(pid,p_init_lep,{},ParticleStatus::initial_state);
	double Ebeam=m_beam->MaxEnergy();
	Particle beam(pid,{Ebeam,0,0,Ebeam},{},ParticleStatus::beam);

	// Primary must always be added first
	size_t primary=event.addVertex(EventVertex::StatusCode::primary);
	size_t beamVertex=event.addVertex(EventVertex::StatusCode::beam);
	event.addParticle(beam,Particle::OPEN_END,beamVertex);
	event.addParticle(init_lep,beamVertex,primary);

#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
#endif
}

double BeamMapper::GenerateWeight(const Event& event, std::vector<double> &rans) {
	const FourVector& p_init_lep=event.leptonIn().Momentum();

    // TODO: Should Masses().back() be the mass of the final state hadronic system or the initial?
    double wgt = m_beam->GenerateWeight(*m_beam->BeamIDs().begin(), p_init_lep, rans,
                                      (Smin() - Masses()[1]) / (2 * sqrt(Masses()[1])));
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
    spdlog::trace("  Beam weight = {}", wgt);
#endif
    return 1.0 / wgt;
}

size_t BeamMapper::NDims() const {
    return static_cast<size_t>(m_beam->NVariables());
}