#include "Achilles/CascadeInteractions/PionAbsorption.hh"
#include "Achilles/Event.hh"
#include "Achilles/Random.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

PionAbsorption::PionAbsorption(const YAML::Node &) {
    // Maps of pion + nucleon as key
    // value is vector of pairs (absorption partner, outgoing final states)
    states[{PID::pionp(), PID::proton()}] = {{PID::neutron(), {PID::proton(), PID::proton()}}};

    states[{PID::pionp(), PID::neutron()}] = {{PID::neutron(), {PID::proton(), PID::neutron()}},
                                              {PID::neutron(), {PID::neutron(), PID::proton()}},
                                              {PID::proton(), {PID::proton(), PID::proton()}}};

    states[{-PID::pionp(), PID::proton()}] = {{PID::neutron(), {PID::neutron(), PID::neutron()}},
                                              {PID::proton(), {PID::neutron(), PID::proton()}},
                                              {PID::proton(), {PID::proton(), PID::neutron()}}};

    states[{-PID::pionp(), PID::neutron()}] = {{PID::proton(), {PID::neutron(), PID::neutron()}}};

    states[{PID::pion0(), PID::proton()}] = {{PID::proton(), {PID::proton(), PID::proton()}},
                                             {PID::neutron(), {PID::proton(), PID::neutron()}},
                                             {PID::neutron(), {PID::neutron(), PID::proton()}}};

    states[{PID::pion0(), PID::neutron()}] = {{PID::neutron(), {PID::neutron(), PID::neutron()}},
                                              {PID::proton(), {PID::neutron(), PID::proton()}},
                                              {PID::proton(), {PID::proton(), PID::neutron()}}};
}

std::vector<std::pair<PID, PID>> PionAbsorption::InitialStates() const {
    return {{PID::pion0(), PID::proton()},  {PID::pionp(), PID::proton()},
            {-PID::pionp(), PID::proton()}, {PID::pion0(), PID::neutron()},
            {PID::pionp(), PID::neutron()}, {-PID::pionp(), PID::neutron()}};
}

bool PionAbsorptionOneStep::AllowedAbsorption(Event &event, size_t part1, size_t part2) const {
    auto closest = FindClosest(event, part1, part2);
    //spdlog::debug("closest: {}, {}", closest.first, closest.second);
    //spdlog::debug("part1: {}, part2: {}", event.Hadrons()[part1].ID(), event.Hadrons()[part2].ID());

    auto modes = states.at({event.Hadrons()[part1].ID(), event.Hadrons()[part2].ID()});

    m_states.clear();
    for(const auto &mode : modes) {
        if(mode.absorption_partner == PID::proton()) {
            m_states.push_back({closest.first, mode.outgoing});
        } else {
            m_states.push_back({closest.second, mode.outgoing});
        }
    }

    return true;
}

std::pair<size_t, size_t> PionAbsorptionOneStep::FindClosest(Event &event, size_t part1,
                                                             size_t part2) const {
    // Particle 2 is the incoming nucleon
    const auto &particle2 = event.Hadrons()[part2];

    double shortest_p_distance = 10000; // fm^2
    double shortest_n_distance = 10000; // fm^2
    size_t closest_p_idx = SIZE_MAX, closest_n_idx = SIZE_MAX;

    for(std::size_t i = 0; i < event.Hadrons().size(); ++i) {
        if(i == part1 || i == part2) continue;
        if(!event.Hadrons()[i].IsBackground()) continue;
        auto distance = (event.Hadrons()[i].Position() - particle2.Position()).Magnitude2();
        //spdlog::debug("distance = {}", distance);
        if(distance < shortest_p_distance && event.Hadrons()[i].ID() == PID::proton()) {
            shortest_p_distance = distance;
            closest_p_idx = i;
        } else if(distance < shortest_n_distance && event.Hadrons()[i].ID() == PID::neutron()) {
            shortest_n_distance = distance;
            closest_n_idx = i;
        }
    }
    //spdlog::debug("p idx = {}, n idx = {}", closest_p_idx, closest_n_idx);
    return {closest_p_idx, closest_n_idx};
}

InteractionResults PionAbsorption::CrossSection(Event &event, size_t part1, size_t part2) const {
    const auto &particle2 = event.Hadrons()[part2];

    InteractionResults results;

    // Fills m_states which is a vector of <absorption partner idx, <out PID1, outPID2>>
    if(!AllowedAbsorption(event, part1, part2)) return results;

    auto oset_abs_xsec = Oset.AbsCrossSection(event, part1, part2);
    //spdlog::debug("Oset abs xsec = {}", oset_abs_xsec);
    //oset_abs_xsec = 0.0;

    // Nuclear Physics A568 (1994) 855-872 Table 1
    auto opposite_isospin_xsec = (5. / 6.) * oset_abs_xsec;
    auto same_isospin_xsec = (1. / 6.) * oset_abs_xsec;

    auto same_isospin_counter = 0;
    auto opp_isospin_counter = 0;

    for(const auto &state : m_states) {
        size_t abs_partner_idx = state.first;
        //spdlog::debug("abs idx = {}", abs_partner_idx);
        if(abs_partner_idx == SIZE_MAX) continue;
        absorption_partners.push_back(event.Hadrons()[abs_partner_idx]);
        if(event.Hadrons()[abs_partner_idx].ID() == particle2.ID()) {
            same_isospin_counter += 1;
        } else {
            opp_isospin_counter += 1;
        }
    }

    for(const auto &state : m_states) {
        size_t abs_partner_idx = state.first;
        if(abs_partner_idx == SIZE_MAX) continue;

        if(event.Hadrons()[abs_partner_idx].ID() == particle2.ID()) {
            results.push_back({state.second, same_isospin_xsec / same_isospin_counter});
            /*spdlog::debug("Absorption xsec: {} + {} + {} -> {} + {} = {}",
                          event.Hadrons()[part1].ID(), event.Hadrons()[part2].ID(),
                          event.Hadrons()[abs_partner_idx].ID(), state.second[0], state.second[1],
                          same_isospin_xsec / same_isospin_counter);*/
        } else {
            results.push_back({state.second, opposite_isospin_xsec / opp_isospin_counter});
            /*spdlog::debug("Absorption xsec: {} + {} + {} -> {} + {} = {}",
                          event.Hadrons()[part1].ID(), event.Hadrons()[part2].ID(),
                          event.Hadrons()[abs_partner_idx].ID(), state.second[0], state.second[1],
                          opposite_isospin_xsec / opp_isospin_counter); */
        }
    }

    return results;
}

std::vector<Particle> PionAbsorptionOneStep::GenerateMomentum(const Particle &particle1,
                                                              const Particle &particle2,
                                                              const std::vector<PID> &out_pids,
                                                              Random &random) const {
    Particle absorption_partner;

    auto partner_charge =
        (ParticleInfo(out_pids[0]).IntCharge() + ParticleInfo(out_pids[1]).IntCharge() -
         particle1.Info().IntCharge() - particle2.Info().IntCharge()) /
        3;

    for(const auto &partner : absorption_partners) {
        if(partner.Info().IntCharge() / 3 == partner_charge) {
            absorption_partner = partner;
            break;
        }
    }

    absorption_partner.Status() = ParticleStatus::absorption_partner;

    spdlog::debug("We chose pion absorption!");
    // Boost to center of mass
    ThreeVector boostCM =
        (particle1.Momentum() + particle2.Momentum() + absorption_partner.Momentum()).BoostVector();
    spdlog::debug("{}: {}", particle1.ID(), particle1.Momentum());
    spdlog::debug("{}: {}", particle2.ID(), particle2.Momentum());
    spdlog::debug("{}: {}", absorption_partner.ID(), absorption_partner.Momentum());

    FourVector p1CM = particle1.Momentum().Boost(-boostCM);
    FourVector p2CM = particle2.Momentum().Boost(-boostCM);
    FourVector p3CM = absorption_partner.Momentum().Boost(-boostCM);

    FourVector pTotalCM = p1CM + p2CM + p3CM;

    auto s = pTotalCM.M2();
    auto sqrts = sqrt(s);

    auto ma = ParticleInfo(out_pids[0]).Mass();
    auto mb = ParticleInfo(out_pids[1]).Mass();

    const double Eacms = sqrts / 2 * (1 + ma * ma / s - mb * mb / s);
    const double Ebcms = sqrts / 2 * (1 + mb * mb / s - ma * ma / s);
    auto lambda = sqrt(pow(s - ma * ma - mb * mb, 2) - 4 * ma * ma * mb * mb);
    auto pfCMS = lambda / 2 / sqrts;

    std::vector<double> rans(2);
    random.Generate(rans);

    double cos_cms = 2 * rans[0] - 1;
    double sin_cms = sqrt(1. - cos_cms * cos_cms);
    double phi_cms = 2 * M_PI * rans[1];
    double cosphi_cms = cos(phi_cms);
    double sinphi_cms = sin(phi_cms);

    FourVector paOut = FourVector(Eacms, pfCMS * sin_cms * cosphi_cms, pfCMS * sin_cms * sinphi_cms,
                                  pfCMS * cos_cms);
    FourVector pbOut = FourVector(Ebcms, -pfCMS * sin_cms * cosphi_cms, -pfCMS * sin_cms * sinphi_cms,
                                  -pfCMS * cos_cms);

    paOut = paOut.Boost(boostCM);
    pbOut = pbOut.Boost(boostCM);

    spdlog::debug("out1: {}, {}", out_pids[0], paOut.Momentum());
    spdlog::debug("out2: {}, {}", out_pids[1], pbOut.Momentum());

    return {{out_pids[0], paOut, particle1.Position()}, {out_pids[1], pbOut, particle2.Position()}};
}

bool PionAbsorptionTwoStep::AllowedAbsorption(Event &, size_t, size_t) const {
    return false;
}

std::vector<Particle> PionAbsorptionTwoStep::GenerateMomentum(const Particle &, const Particle &,
                                                              const std::vector<PID> &,
                                                              Random &) const {
    return {};
}
