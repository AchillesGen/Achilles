#include "Achilles/ProcessInfo.hh"

size_t achilles::ProcessInfo::FinalStateMultiplicity() const {
    return m_leptonic.second.size() + m_hadronic.second.size();
}

size_t achilles::ProcessInfo::Multiplicity() const {
    return 1 + m_hadronic.first.size() + FinalStateMultiplicity();
}

std::vector<double> achilles::ProcessInfo::Masses() const {
    std::vector<double> masses;

    // Get Leptonic current masses
    for(const auto &part : m_leptonic.second) {
        masses.push_back(pow(ParticleInfo(part).Mass(), 2));
    }

    // Get Hadronic current masses
    for(const auto &part : m_hadronic.second) {
        masses.push_back(pow(ParticleInfo(part).Mass(), 2));
    }

    return masses;
}

std::vector<long> achilles::ProcessInfo::Ids() const {
    std::vector<long> ids;

    // Get inital lepton id
    ids.push_back(m_leptonic.first.AsInt());

    // Get initial hadronic ids
    for(const auto &part : m_hadronic.first) ids.push_back(part.AsInt());

    // Get remaining leptonic ids
    for(const auto &part : m_leptonic.second) ids.push_back(part.AsInt());

    // Get remaining hadronic ids
    for(const auto &part : m_hadronic.second) ids.push_back(part.AsInt());

    // Get spectator ids
    for(const auto &part : m_spectator) ids.push_back(part.AsInt());

    return ids;
}

int achilles::ProcessInfo::LeptonicCharge() const {
    int charge = -ParticleInfo(m_leptonic.first).IntCharge();
    for(const auto &part : m_leptonic.second) { charge += ParticleInfo(part).IntCharge(); }
    return charge / 3;
}

bool achilles::ProcessInfo::SaveState(std::ostream &os) const {
    os << m_leptonic.first << " " << m_leptonic.second.size() << " ";
    for(const auto &part : m_leptonic.second) os << part << " ";
    os << m_hadronic.first.size() << " ";
    for(const auto &part : m_hadronic.first) os << part << " ";
    os << m_hadronic.second.size() << " ";
    for(const auto &part : m_hadronic.second) os << part << " ";
    os << m_spectator.size() << " ";
    for(const auto &part : m_spectator) os << part << " ";
    return true;
}

bool achilles::ProcessInfo::LoadState(std::istream &is) {
    is >> m_leptonic.first;
    size_t size;
    is >> size;
    m_leptonic.second.resize(size);
    for(auto &part : m_leptonic.second) is >> part;
    is >> size;
    m_hadronic.first.resize(size);
    for(auto &part : m_hadronic.first) is >> part;
    is >> size;
    m_hadronic.second.resize(size);
    for(auto &part : m_hadronic.second) is >> part;
    is >> size;
    m_spectator.resize(size);
    for(auto &part : m_spectator) is >> part;
    return true;
}
