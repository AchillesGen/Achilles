#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Beams.hh"
#include "Achilles/NuclearModel.hh"

using achilles::Event;

Event::Event(std::shared_ptr<Nucleus> nuc,
             std::vector<FourVector> mom, double vwgt)
        : m_nuc{std::move(nuc)}, m_mom{std::move(mom)}, m_vWgt{std::move(vwgt)} {
    m_nuc -> GenerateConfig();
}

// bool Event::ValidateEvent(size_t imatrix) const {
//     spdlog::trace("Number of momentums = {}", m_mom.size());
//     spdlog::trace("Number of states = {}", m_me[imatrix].inital_state.size() + m_me[imatrix].final_state.size());
//     return m_mom.size() == m_me[imatrix].inital_state.size() + m_me[imatrix].final_state.size();
// }

// TODO: Make this cleaner and handle multiple nucleons in initial state
void Event::InitializeLeptons(const Process_Info &process) {
    // Setup leptons
    bool initial_state = true;
    for(const auto &elm : process.m_mom_map) {
        if(ParticleInfo(elm.second).IsLepton() || ParticleInfo(elm.second).IsVector()) {
            m_leptons.emplace_back(ParticleInfo(elm.second), m_mom[elm.first]);
            if(initial_state) {
                m_leptons.back().Status() = ParticleStatus::initial_state;
                initial_state = false;
            } else {
                m_leptons.back().Status() = ParticleStatus::final_state;
            }
        }
    }
}

void Event::InitializeHadrons(const Process_Info &process) {
    // Get number of incoming and ougoing hadrons
    const size_t nhin = process.state.first.size();
    const size_t nhout = process.state.second.size();

    // Check if the process is coherent
    if(ParticleInfo(process.state.first[0]).IsNucleus()) {
        // Remove all nucleons
        m_nuc -> Nucleons().clear(); 

        // Setup initial and final state nucleus
        const auto id = process.state.first[0];
        Particle initial = Particle(id, m_mom[0]);
        initial.Status() = ParticleStatus::initial_state;
        m_nuc -> Nucleons().push_back(initial);
        Particle final_state = Particle(id, m_mom[2]);
        final_state.Status() = ParticleStatus::final_state;
        m_nuc -> Nucleons().push_back(final_state);
        return;
    }

    // Check if the process is MEC or QE-MEC interference
    if(nhin > 1 && nhin == nhout) {
        for(size_t i = 0; i < nhin; ++i) {
            const auto id = process.state.first[i];
            auto part = m_nuc -> SelectNucleon(id);
            part.Status() = ParticleStatus::initial_state;
            part.Momentum() = m_mom[i];
            auto final_state = Particle(process.state.second[i], m_mom[nhin+1+i]);
            final_state.Status() = ParticleStatus::propagating;
            final_state.Position() = part.Position();
            m_nuc -> Nucleons().push_back(final_state);
        }
        return;
    }

    // Handle hadrons
    ThreeVector position;
    for(const auto &elm : process.m_mom_map) {
        if(!ParticleInfo(elm.second).IsHadron()) continue;
        if(elm.first < nhin) {
            auto &part = m_nuc -> SelectNucleon(elm.second);
            part.Status() = ParticleStatus::initial_state;
            part.Momentum() = m_mom[elm.first];
            position = part.Position();
        } else {
            auto part = Particle(elm.second, m_mom[elm.first]); 
            part.Status() = ParticleStatus::propagating;
            part.Position() = position;
            m_nuc -> Nucleons().push_back(part);
        }
    }
}

void Event::Finalize() {
    size_t nA = 0, nZ = 0;
    for(auto it = m_nuc -> Nucleons().begin(); it != m_nuc -> Nucleons().end(); ) {
        if(it -> Status() == ParticleStatus::background) {
            if(it -> ID() == PID::proton()) nZ++;
            nA++;
            it = m_nuc -> Nucleons().erase(it);
        } else {
            ++it;
        }
    }

    m_remnant = NuclearRemnant(nA, nZ);
}

achilles::vParticles Event::Particles() const {
    vParticles result;
    result.insert(result.end(), Hadrons().begin(), Hadrons().end());
    result.insert(result.end(), m_leptons.begin(), m_leptons.end());
    return result;
}

const achilles::vParticles& Event::Hadrons() const {
    return m_nuc -> Nucleons();
}

achilles::vParticles& Event::Hadrons() {
    return m_nuc -> Nucleons();
}

void Event::CalcWeight() {
    // if(!ValidateEvent(0))
    //     throw std::runtime_error("Phase space and Matrix element have different number of particles");
    m_wgt = m_vWgt*m_meWgt;
}

bool Event::CalcTotalCrossSection(const std::vector<double> &xsec) {
    m_meWgt = std::accumulate(xsec.begin(), xsec.end(), 0.0, std::plus<>());
    spdlog::debug("Total xsec = {}", m_meWgt);
    return m_meWgt > 0 ? true : false;
}

// size_t Event::SelectNucleon() const {
//     std::vector<double> probs = EventProbs();
//     double rand = Random::Instance().Uniform(0.0, 1.0);
//     return static_cast<size_t>(std::distance(probs.begin(),
//                                std::lower_bound(probs.begin(), probs.end(), rand)))-1;
// }

// std::vector<double> Event::EventProbs() const {
//     std::vector<double> probs;
//     probs.push_back(0.0);
//     double cumulative = 0.0;
//     for(const auto & m : m_me) {
//         cumulative += m;
//         probs.emplace_back(cumulative / m_meWgt);
//     }
//     return probs;
// }

void Event::Rotate(const std::array<double,9>& rot_mat) {
    for (auto& particle: m_nuc -> Nucleons()){ particle.Rotate(rot_mat); }
    for (auto& particle: m_leptons){ particle.Rotate(rot_mat); }
}
