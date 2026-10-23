#include "Achilles/FoldingCorrection.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/Random.hh"

#include "spdlog/spdlog.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

using achilles::FoldingCorrection;

namespace {
constexpr double pbin_width = 50;  // MeV
constexpr std::size_t npbins = 60; // up to 3 GeV
constexpr double tmax = 30, dt = 0.1;
} // namespace

FoldingCorrection::Mode FoldingCorrection::ParseMode(const std::string &name) {
    if(name == "Full") return Mode::Full;
    if(name == "Limited") return Mode::Limited;
    if(name == "Collided") return Mode::Collided;
    throw std::runtime_error("FoldingCorrection: unknown mode " + name +
                             " (expected Full, Limited or Collided)");
}

FoldingCorrection::FoldingCorrection(Mode mode, double max_delta, std::size_t paths_per_bin,
                                     std::size_t npauli)
    : m_mode{mode}, m_max_delta{max_delta}, m_paths{paths_per_bin}, m_npauli{npauli},
      m_radii{0, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5} {
    const auto nt = static_cast<std::size_t>(tmax / dt) + 1;
    for(std::size_t k = 0; k < nt; ++k) m_times.push_back(static_cast<double>(k) * dt);
    m_offsets.push_back(0);
    for(double offset = 0.2; offset < 4000; offset *= 1.04) m_offsets.push_back(offset);
    m_ensembles.resize(2 * m_radii.size() * npbins);
    m_built.assign(2 * npbins, false);
}

std::size_t FoldingCorrection::Stratum(double radius) const {
    std::size_t k = 0;
    while(k + 1 < m_radii.size() && radius >= m_radii[k + 1]) ++k;
    return k;
}

std::size_t FoldingCorrection::Index(bool proton, std::size_t stratum, std::size_t pbin) const {
    return ((proton ? 1 : 0) * m_radii.size() + stratum) * npbins + pbin;
}

void FoldingCorrection::Build(Cascade &cascade, const std::shared_ptr<Nucleus> &nucleus,
                              bool proton, std::size_t pbin) {
    spdlog::debug("FoldingCorrection: building tables for {} momentum bin {}",
                  proton ? "proton" : "neutron", pbin);
    const std::size_t nt = m_times.size();
    for(std::size_t s = 0; s < m_radii.size(); ++s) {
        auto &ens = m_ensembles[Index(proton, s, pbin)];
        ens.survival.assign(nt, 0);
        ens.cumulative.assign(nt, 0); // accumulates <sqrt S> first
    }
    const PID pid = proton ? PID::proton() : PID::neutron();
    for(std::size_t path = 0; path < m_paths; ++path) {
        // The interactions read in-medium quantities from the event's nucleus
        Event event;
        event.CurrentNucleus() = nucleus;
        event.Hadrons() = nucleus->GenerateConfig();
        std::vector<std::size_t> candidates;
        for(std::size_t i = 0; i < event.Hadrons().size(); ++i)
            if(event.Hadrons()[i].ID() == pid) candidates.push_back(i);
        const std::size_t idx = Random::Instance().Pick(candidates);
        auto &part = event.Hadrons()[idx];
        const double p =
            (static_cast<double>(pbin) + Random::Instance().Uniform(0.0, 1.0)) * pbin_width;
        const double cost = Random::Instance().Uniform(-1.0, 1.0),
                     sint = std::sqrt(1 - cost * cost);
        const double phi = Random::Instance().Uniform(0.0, 2 * M_PI), mass = part.Info().Mass();
        part.SetMomentum({std::sqrt(p * p + mass * mass), p * sint * std::cos(phi),
                          p * sint * std::sin(phi), p * cost});
        part.Status() = ParticleStatus::propagating;

        const auto survival = cascade.SurvivalProfile(event, nucleus.get(), idx, m_times, m_npauli);
        auto &ens = m_ensembles[Index(proton, Stratum(part.Position().Magnitude()), pbin)];
        for(std::size_t k = 0; k < nt; ++k) {
            ens.survival[k] += survival[k];
            ens.cumulative[k] += std::sqrt(survival[k]);
        }
        ens.atom += std::sqrt(survival.back());
        ens.transparency += survival.back();
        ens.n += 1;
    }

    // Cumulative of the folding function of <sqrt S> in |offset|:
    // G(D) = A(inf) + (2/pi) int dt (A(t) - A(inf)) sin(D t)/t
    for(std::size_t s = 0; s < m_radii.size(); ++s) {
        auto &ens = m_ensembles[Index(proton, s, pbin)];
        if(ens.n == 0) continue;
        std::vector<double> amplitude(nt);
        for(std::size_t k = 0; k < nt; ++k) {
            ens.survival[k] /= ens.n;
            amplitude[k] = ens.cumulative[k] / ens.n;
        }
        ens.atom /= ens.n;
        ens.transparency /= ens.n;
        ens.cumulative.assign(m_offsets.size(), 0);
        for(std::size_t k = 0; k < m_offsets.size(); ++k) {
            const double omega = m_offsets[k] / Constant::HBARC;
            double g = 0;
            for(std::size_t i = 1; i < nt; ++i) {
                const double t = (static_cast<double>(i) - 0.5) * dt;
                const double a = 0.5 * (amplitude[i - 1] + amplitude[i]) - ens.atom;
                g += (omega > 0 ? std::sin(omega * t) / t : 0.0) * a * dt;
            }
            ens.cumulative[k] =
                std::max(ens.atom + 2 / M_PI * g, k > 0 ? ens.cumulative[k - 1] : 0.0);
        }
    }
}

double FoldingCorrection::Quantile(const Ensemble &ens, double c) const {
    const auto &G = ens.cumulative;
    if(G.empty() || c <= G[0]) return 0;
    const auto it = std::lower_bound(G.begin(), G.end(), c);
    const auto k = static_cast<std::size_t>(it - G.begin());
    if(k >= G.size()) return m_offsets.back();
    const double f = (c - G[k - 1]) / std::max(1e-15, G[k] - G[k - 1]);
    return m_offsets[k - 1] + f * (m_offsets[k] - m_offsets[k - 1]);
}

double FoldingCorrection::Survival(const Ensemble &ens, double t) const {
    const double x = std::min(t, tmax) / dt;
    const auto i = std::min(static_cast<std::size_t>(x), m_times.size() - 2);
    return ens.survival[i] + (x - static_cast<double>(i)) * (ens.survival[i + 1] - ens.survival[i]);
}

double FoldingCorrection::Cumulative(const Ensemble &ens, double delta) const {
    const auto &G = ens.cumulative;
    if(delta >= m_offsets.back()) return G.back();
    const auto it = std::upper_bound(m_offsets.begin(), m_offsets.end(), delta);
    const auto k = static_cast<std::size_t>(it - m_offsets.begin());
    const double f = (delta - m_offsets[k - 1]) / (m_offsets[k] - m_offsets[k - 1]);
    return G[k - 1] + f * (G[k] - G[k - 1]);
}

const FoldingCorrection::Ensemble *FoldingCorrection::Table(Cascade &cascade,
                                                            const std::shared_ptr<Nucleus> &nucleus,
                                                            const Particle &primary) {
    const bool proton = primary.ID() == PID::proton();
    const auto pbin =
        std::min(static_cast<std::size_t>(primary.Momentum().P() / pbin_width), npbins - 1);
    const std::size_t built = (proton ? 1 : 0) * npbins + pbin;
    if(!m_built[built]) {
        Build(cascade, nucleus, proton, pbin);
        m_built[built] = true;
    }
    const auto &ens = m_ensembles[Index(proton, Stratum(primary.Position().Magnitude()), pbin)];
    return ens.n > 0 ? &ens : nullptr;
}

double FoldingCorrection::Delta(const Ensemble &ens, double t1, double delta_max) const {
    // Comonotone coupling: early collisions get the large offsets, escapes the atom
    const double rank =
        t1 >= 0 ? 1 - Survival(ens, t1)
                : (1 - ens.transparency) + Random::Instance().Uniform(0.0, 1.0) * ens.transparency;
    double c = 1 - rank;
    if(c <= ens.atom || delta_max <= 0) return 0;
    if(delta_max < m_offsets.back() && ens.atom < 1)
        c = ens.atom + (c - ens.atom) * (Cumulative(ens, delta_max) - ens.atom) / (1 - ens.atom);
    const double sign = Random::Instance().Uniform(0.0, 1.0) < 0.5 ? -1 : 1;
    return sign * std::min(Quantile(ens, c), delta_max);
}

double FoldingCorrection::Apply(Cascade &cascade, Event &event, std::size_t idx,
                                const Particle &primary) {
    const auto &nucleus = event.CurrentNucleus();
    if(!nucleus || !primary.Info().IsNucleon()) return 0;
    const double p = primary.Momentum().P();

    // Real part: first moment of the folding function, the potential at the vertex
    double shift = 0;
    if(auto potential = nucleus->GetPotential())
        shift = potential->operator()(p, primary.Position().Magnitude()).rvector;

    const auto &first = cascade.FirstCollisions();
    const auto hit = first.find(idx);
    const double t1 = hit == first.end() ? -1 : hit->second.time;

    double delta_max = m_max_delta;
    if(m_mode != Mode::Full) {
        const auto &leptons = event.Leptons();
        auto with_status = [&](ParticleStatus status) {
            return std::find_if(leptons.begin(), leptons.end(), [&](const Particle &part) {
                return part.Info().IsLepton() && part.Status() == status;
            });
        };
        const auto beam = with_status(ParticleStatus::initial_state);
        const auto out = with_status(ParticleStatus::final_state);
        if(beam == leptons.end() || out == leptons.end()) return 0;
        // Keep 0 < omega + offset < E' - m for both signs of Delta
        const double omega = beam->Momentum().E() - out->Momentum().E();
        const double room = out->Momentum().E() - out->Info().Mass();
        delta_max = std::min(delta_max, std::min(omega + shift, room - shift));
        if(delta_max < 0) return 0;
    }
    if(m_mode == Mode::Collided)
        delta_max =
            hit == first.end() ? 0 : std::min(delta_max, std::max(0.0, hit->second.energy_loss));

    const auto *ens = Table(cascade, nucleus, primary);
    const double offset = shift + (ens ? Delta(*ens, t1, delta_max) : 0);
    return ShiftLepton(event, offset) ? offset : 0;
}

bool FoldingCorrection::ShiftLepton(Event &event, double offset) {
    auto shift = [&](Particle &lepton) {
        const double energy = lepton.Momentum().E() - offset, mass = lepton.Info().Mass();
        if(energy <= mass) return false;
        const auto dir = lepton.Momentum().Vec3() / lepton.Momentum().P();
        const double pmag = std::sqrt(energy * energy - mass * mass);
        lepton.SetMomentum({energy, pmag * dir[0], pmag * dir[1], pmag * dir[2]});
        return true;
    };
    auto is_outgoing_lepton = [](const Particle &part) {
        return part.Info().IsLepton() && part.Status() == ParticleStatus::final_state;
    };

    auto &leptons = event.Leptons();
    auto lepton = std::find_if(leptons.begin(), leptons.end(), is_outgoing_lepton);
    if(lepton == leptons.end()) return false;
    const Particle original = *lepton;
    if(!shift(*lepton)) return false;

    // The writer reads the history: update the same lepton in the primary vertex
    if(auto *primary = event.History().Primary()) {
        for(auto &part : primary->ParticlesOut()) {
            if(part.ID() == original.ID() && part.Momentum().Approx(original.Momentum(), 1e-8)) {
                part.SetMomentum(lepton->Momentum());
                break;
            }
        }
    }
    return true;
}
