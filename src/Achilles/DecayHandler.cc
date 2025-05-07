#include "Achilles/DecayHandler.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Random.hh"

#include "spdlog/spdlog.h"

using achilles::DecayHandler;

DecayHandler::DecayHandler(const std::string &filename, double tolerance) {
    auto node = YAML::LoadFile(filename);

    for(const auto &particle : node) {
        auto pid = particle.first.as<PID>();
        double br_total = 0;
        for(const auto &decay_mode : particle.second) {
            DecayMode decay{decay_mode["BranchingRatio"].as<double>(),
                            decay_mode["AngularMom"].as<size_t>(),
                            decay_mode["OutIDs"].as<std::vector<PID>>()};
            br_total += decay.branching_ratio;
            std::sort(decay.out_ids.begin(), decay.out_ids.end());
            m_decays[pid].push_back(decay);
        }

        if(std::abs(br_total - 1) > tolerance) {
            spdlog::warn(
                "DecayHandler: Branching ratio not within {} of 1 ({}) for PID {}. Rescaling.",
                tolerance, br_total, pid);
            for(auto &mode : m_decays[pid]) { mode.branching_ratio /= br_total; }
        }
    }
}

std::vector<achilles::Particle> DecayHandler::Decay(const Particle &part) const {
    spdlog::debug("Decaying: {}", part);
    auto boost = part.Momentum().BoostVector();
    double m2 = part.Momentum().M2();

    // Select decay channel
    auto ratios = BranchingRatios(part);
    spdlog::trace("Branching Rations: {}", fmt::join(ratios, ","));
    // Remove kinematically disallowed decays
    size_t tmp = 0;
    for(const auto &mode : m_decays.at(part.ID())) {
        if(mode.out_ids.size() == 2) {
            spdlog::trace("Checking decay {} ({}) -> {} ({}) + {} ({})", part.ID(), sqrt(m2),
                          mode.out_ids[0], ParticleInfo(mode.out_ids[0]).Mass(), mode.out_ids[1],
                          ParticleInfo(mode.out_ids[1]).Mass());
            auto mass1 = ParticleInfo(mode.out_ids[0]).Mass();
            auto mass2 = ParticleInfo(mode.out_ids[1]).Mass();
            if(mass1 + mass2 > sqrt(m2)) {
                spdlog::debug("Removing decay {} ({}) -> {} ({}) + {} ({})", part.ID(), sqrt(m2),
                              mode.out_ids[0], mass1, mode.out_ids[1], mass2);
                ratios[tmp] = 0;
            }
        }
        tmp++;
    }
    auto idx = Random::Instance().SelectIndex(ratios);
    auto mode = m_decays.at(part.ID())[idx];
    spdlog::trace("Decaying to mode: {}", mode);

    std::vector<Particle> outgoing;
    if(mode.out_ids.size() == 2)
        outgoing = TwoBodyDecay(m2, mode.out_ids, mode.angular_mom);
    else { throw std::runtime_error("Three body and above decays are not implemented yet"); }

    // Boost back to lab frame
    for(auto &particle : outgoing) { particle.Momentum() = particle.Momentum().Boost(boost); }
    return outgoing;
}

std::vector<double> DecayHandler::BranchingRatios(const Particle &part) const {
    std::vector<double> ratios;
    for(const auto &mode : m_decays.at(part.ID())) { ratios.push_back(mode.branching_ratio); }
    return ratios;
}

std::vector<achilles::DecayMode> DecayHandler::AllowedDecays(PID pid) const {
    return m_decays.at(pid);
}

double DecayHandler::BranchingRatio(PID res, std::vector<PID> out) const {
    spdlog::trace("Getting BranchingRatio for {} -> {}", res, fmt::join(out, ","));
    const auto &decays = m_decays.at(res);
    std::sort(out.begin(), out.end());
    for(const auto &decay : decays) {
        if(decay.out_ids == out) return decay.branching_ratio;
    }
    return 0;
}

std::vector<achilles::Particle>
DecayHandler::TwoBodyDecay(double mass2, const std::vector<PID> &pids, size_t) const {
    // TODO: Add in angular momentum (i.e. learn more about how Sherpa handles this)
    double sqrts = sqrt(mass2);

    auto ma = ParticleInfo(pids[0]).Mass();
    auto mb = ParticleInfo(pids[1]).Mass();

    const double Ea = sqrts / 2. * (1 + ma * ma / mass2 - mb * mb / mass2);
    const double Eb = sqrts / 2. * (1 + mb * mb / mass2 - ma * ma / mass2);
    auto lambda = sqrt(pow(mass2 - ma * ma - mb * mb, 2) - 4 * ma * ma * mb * mb);
    auto pf = lambda / 2 / sqrts;

    std::vector<double> rans(2);
    Random::Instance().Generate(rans);

    double cost = 2 * rans[0] - 1;
    double sint = sqrt(1 - cost * cost);
    double phi = 2 * M_PI * rans[1];

    Particle part1{pids[0], {Ea, pf * sint * cos(phi), pf * sint * sin(phi), pf * cost}};
    Particle part2{pids[1], {Eb, -pf * sint * cos(phi), -pf * sint * sin(phi), -pf * cost}};

    return {part1, part2};
}
