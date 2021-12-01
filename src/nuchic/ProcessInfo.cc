#include "nuchic/ProcessInfo.hh"

size_t nuchic::Process_Info::Multiplicity() const {
    return m_ids.size() + m_states.begin() -> first.size() + m_states.begin() -> second.size();
}

std::vector<double> nuchic::Process_Info::Masses() const {
    std::vector<double> masses;

    // Get Leptonic current masses
    for(size_t i = 1; i < m_ids.size(); ++i) {
        masses.push_back(pow(ParticleInfo(m_ids[i]).Mass(), 2));
    }

    // Get Hadronic current masses
    for(size_t i = 0; i < m_states.begin() -> second.size(); ++i) {
        masses.push_back(pow(ParticleInfo(m_states.begin() -> second[i]).Mass(), 2));
    }

    return masses;
}
