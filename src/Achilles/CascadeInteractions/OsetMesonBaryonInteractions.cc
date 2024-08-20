#include "Achilles/CascadeInteractions/OsetMesonBaryonInteractions.hh"
#include "Achilles/Event.hh"
#include "Achilles/Poincare.hh"
#include "Achilles/Random.hh"

using namespace achilles;

OsetMesonBaryonInteraction::OsetMesonBaryonInteraction(const YAML::Node &) {
    // Maps of pion + nucleon PIDs as key
    // value is vector of pairs of allowed outgoing PIDs
    in_out_states[{PID::pionp(), PID::proton()}] = {{PID::pionp(), PID::proton()}};

    in_out_states[{PID::pionp(), PID::neutron()}] = {{PID::pionp(), PID::neutron()},
                                                     {PID::pion0(), PID::proton()}};

    in_out_states[{PID::pion0(), PID::proton()}] = {{PID::pion0(), PID::proton()},
                                                    {PID::pionp(), PID::neutron()}};

    in_out_states[{PID::pion0(), PID::neutron()}] = {{PID::pion0(), PID::neutron()},
                                                     {-PID::pionp(), PID::proton()}};

    in_out_states[{-PID::pionp(), PID::proton()}] = {{-PID::pionp(), PID::proton()},
                                                     {PID::pion0(), PID::neutron()}};

    in_out_states[{-PID::pionp(), PID::neutron()}] = {{-PID::pionp(), PID::neutron()}};
}

InteractionResults OsetMesonBaryonInteraction::CrossSection(Event &event, size_t part1,
                                                            size_t part2) const {
    const auto &particle1 = event.Hadrons()[part1];
    const auto &particle2 = event.Hadrons()[part2];

    auto sqrts = (particle1.Momentum() + particle2.Momentum()).M();

    InteractionResults results;

    auto QECrossSections = Oset.QECrossSection(event, part1, part2);


    spdlog::debug("Incoming:");
    spdlog::debug("{}", particle1);
    spdlog::debug("{}", particle2);
    // Vector of pairs of outgoing states
    auto outgoing_states = in_out_states.at({particle1.ID(), particle2.ID()});
    for(const auto &state : outgoing_states) {
        spdlog::debug("Outgoing:");
        spdlog::debug("{} + {}", state.first, state.second);
        PID outgoing_pi_PID = state.first;

        auto outgoing_masses =
            ParticleInfo(outgoing_pi_PID).Mass() + ParticleInfo(state.second).Mass();
        if(sqrts < outgoing_masses) continue;

        double CS = QECrossSections[{particle1.ID(), outgoing_pi_PID}];
        if(CS == 0.) continue;
        results.push_back({{outgoing_pi_PID, state.second}, CS});
    }

    return results;
}

std::vector<std::pair<PID, PID>> OsetMesonBaryonInteraction::InitialStates() const {
    return {{PID::pion0(), PID::proton()},  {PID::pionp(), PID::proton()},
            {-PID::pionp(), PID::proton()}, {PID::pion0(), PID::neutron()},
            {PID::pionp(), PID::neutron()}, {-PID::pionp(), PID::neutron()}};
}

std::vector<Particle> OsetMesonBaryonInteraction::GenerateMomentum(const Particle &particle1,
                                                                   const Particle &particle2,
                                                                   const std::vector<PID> &out_pids,
                                                                   Random &random) const {
    // TODO implement S wave and P wave outgoing angles from Oset

    spdlog::debug("We chose pion QE scattering");
    // Boost to center of mass
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    spdlog::debug("{}: {}", particle1.ID(), particle1.Momentum());
    spdlog::debug("{}: {}", particle2.ID(), particle2.Momentum());

    FourVector p1CM = particle1.Momentum().Boost(-boostCM);
    FourVector p2CM = particle2.Momentum().Boost(-boostCM);

    FourVector pTotalCM = p1CM + p2CM;

    auto s = pTotalCM.M2();
    auto sqrts = sqrt(s);

    auto ma = ParticleInfo(out_pids[0]).Mass();
    auto mb = ParticleInfo(out_pids[1]).Mass();
    spdlog::debug("out 1 mass = {}", ma);
    spdlog::debug("out 2 mass = {}", mb);

    const double Eacms = sqrts / 2. * (1. + ma * ma / s - mb * mb / s);
    const double Ebcms = sqrts / 2. * (1. + mb * mb / s - ma * ma / s);
    auto lambda = sqrt(pow(s - ma * ma - mb * mb, 2) - 4. * ma * ma * mb * mb);
    auto pfCMS = lambda / 2. / sqrts;

    std::vector<double> rans(2);
    random.Generate(rans);

    double cos_cms = 2 * rans[0] - 1;
    double sin_cms = sqrt(1. - cos_cms * cos_cms);
    double phi_cms = 2 * M_PI * rans[1];
    double cosphi_cms = cos(phi_cms);
    double sinphi_cms = sin(phi_cms);

    FourVector paOut = FourVector(Eacms, pfCMS * sin_cms * cosphi_cms, pfCMS * sin_cms * sinphi_cms,
                                  pfCMS * cos_cms);
    FourVector pbOut = FourVector(Ebcms, -pfCMS * sin_cms * cosphi_cms,
                                  -pfCMS * sin_cms * sinphi_cms, -pfCMS * cos_cms);

    paOut = paOut.Boost(boostCM);
    pbOut = pbOut.Boost(boostCM);

    spdlog::debug("out1: {}, {}", out_pids[0], paOut.Momentum());
    spdlog::debug("out2: {}, {}", out_pids[1], pbOut.Momentum());

    return {{out_pids[0], paOut, particle1.Position()}, {out_pids[1], pbOut, particle2.Position()}};
}
