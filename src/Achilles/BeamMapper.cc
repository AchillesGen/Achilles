// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/BeamMapper.hh"
#include "Achilles/Beams.hh"
#include "Achilles/FourVector.hh"

using achilles::BeamMapper;

void BeamMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) {
    auto beam_id = *m_beam->BeamIDs().begin();
    // TODO: Should Masses().back() be the mass of the final state hadronic system or the initial?
    // Masses() and Smin() are squared masses in MeV^2, so this is an energy.
    const units::Energy emin{(Smin() - Masses()[1]) / (2 * sqrt(Masses()[1]))};
    point[m_idx] = m_beam->Flux(beam_id, rans, emin);
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
#endif
}

double BeamMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) {
    auto beam_id = *m_beam->BeamIDs().begin();
    // TODO: Should Masses().back() be the mass of the final state hadronic system or the initial?
    const units::Energy emin{(Smin() - Masses()[1]) / (2 * sqrt(Masses()[1]))};
    auto wgt = m_beam->GenerateWeight(beam_id, point[m_idx], rans, emin);
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  Beam weight = {}", wgt);
#endif
    return 1.0 / wgt;
}

size_t BeamMapper::NDims() const {
    return static_cast<size_t>(m_beam->NVariables());
}
