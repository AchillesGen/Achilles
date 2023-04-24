#include "Achilles/ProcessInfo.hh"

size_t achilles::ProcessInfo::Multiplicity() const {
    return 1 + m_leptonic.second.size() + m_hadronic.first.size() + m_hadronic.second.size();
}

std::vector<double> achilles::ProcessInfo::Masses() const {
    std::vector<double> masses;

    // Get Hadronic current masses
    for(const auto &part : m_hadronic.second) {
        masses.push_back(pow(ParticleInfo(part).Mass(), 2));
    }

    // Get Leptonic current masses
    for(const auto &part : m_leptonic.second) {
        masses.push_back(pow(ParticleInfo(part).Mass(), 2));
    }

    return masses;
}

std::vector<long> achilles::ProcessInfo::Ids() const {
    std::vector<long> ids;

    // Get initial hadronic ids
    for(const auto &part : m_hadronic.first) ids.push_back(part.AsInt());

    // Get inital lepton id
    ids.push_back(m_leptonic.first.AsInt());

    // Get remaining hadronic ids
    for(const auto &part : m_hadronic.second) ids.push_back(part.AsInt());

    // Get remaining leptonic ids
    for(const auto &part : m_leptonic.second) ids.push_back(part.AsInt());

    return ids;
}

int achilles::ProcessInfo::LeptonicCharge() const {
    int charge = -ParticleInfo(m_leptonic.first).IntCharge();
    for(const auto &part : m_leptonic.second) { charge += ParticleInfo(part).IntCharge(); }
    return charge / 3;
}
