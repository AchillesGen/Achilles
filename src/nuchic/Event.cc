#include "nuchic/Event.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Beams.hh"

using nuchic::Event;

Event::Event(std::shared_ptr<Nucleus> nuc, std::shared_ptr<Beam> beam,
             const std::vector<double> &rans, double vwgt) 
        : m_nuc{std::move(nuc)}, m_beam{std::move(beam)}, m_vWgt{std::move(vwgt)} {
    m_nuc -> GenerateConfig();
    m_ps.momentum.push_back(m_beam -> Flux(PID::electron(), rans));
}

bool Event::ValidateEvent(size_t imatrix) const {
    return m_ps.momentum.size() == m_me[imatrix].inital_state.size() + m_me[imatrix].final_state.size();
}

void Event::InitializeLeptons(size_t imatrix) {
    if(!ValidateEvent(imatrix))
        throw std::runtime_error("Phase space and Matrix element have different number of particles");

    size_t idx = 0;
    for(auto pid : m_me[imatrix].inital_state) {
        ParticleInfo info(pid);
        if(info.IsLepton()) {
            m_leptons.emplace_back(info, m_ps.momentum[idx]);
            m_leptons.back().SetStatus(ParticleStatus::initial_state);
            ++idx;
        }
    }
    for(auto pid : m_me[imatrix].final_state) {
        ParticleInfo info(pid);
        if(info.IsLepton()) {
            m_leptons.emplace_back(info, m_ps.momentum[idx]);
            m_leptons.back().SetStatus(ParticleStatus::final_state);
            ++idx;
        }
    }
}

void Event::InitializeHadrons(const std::vector<std::array<size_t, 3>> &idxs) {
    for(const auto &idx : idxs) {
        m_nuc -> Nucleons()[idx[0]].SetStatus(ParticleStatus::initial_state);
        m_nuc -> Nucleons()[idx[0]].SetMomentum(m_ps.momentum[idx[1]]);
        Particle outNucleon(m_nuc -> Nucleons()[idx[0]]);
        outNucleon.SetStatus(ParticleStatus::propagating);
        outNucleon.SetMomentum(m_ps.momentum[idx[2]]);
        m_nuc -> Nucleons().push_back(outNucleon);
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

void Event::AddParticle(const Particle &part) {
    if(part.Info().IsHadron()) {
        m_nuc -> Nucleons().push_back(part);
    } else {
        m_leptons.push_back(part);
    }
}

nuchic::vParticles Event::Particles() const {
    vParticles result;
    result.insert(result.end(), Hadrons().begin(), Hadrons().end());
    result.insert(result.end(), m_leptons.begin(), m_leptons.end());
    return result;
}

const nuchic::vParticles& Event::Hadrons() const {
    return m_nuc -> Nucleons();
}

nuchic::vParticles& Event::Hadrons() {
    return m_nuc -> Nucleons();
}

double Event::Weight() const {
    static constexpr double conv = 1e6;

    if(!ValidateEvent(0))
        throw std::runtime_error("Phase space and Matrix element have different number of particles");
    return m_ps.weight*m_vWgt*m_meWgt*conv;
}

bool Event::TotalCrossSection() {
    m_meWgt = std::accumulate(m_me.begin(), m_me.end(), 0.0, AddEvents);
    return m_meWgt > 0 ? true : false;
}

std::vector<double> Event::EventProbs() const {
    std::vector<double> probs;
    probs.push_back(0.0);
    double cumulative = 0.0;
    for(const auto & m : m_me) {
        cumulative += m.weight;
        probs.emplace_back(cumulative / m_meWgt);
    }
    return probs;
}

bool Event::MatrixCompare(const MatrixElementStruct &m, double value) {
    return m.weight < value;
}

double Event::AddEvents(double value, const MatrixElementStruct &m) {
    return value + m.weight;
}
