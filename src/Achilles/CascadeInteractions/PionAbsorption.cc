#include "Achilles/CascadeInteractions/PionAbsorption.hh"
#include "Achilles/Event.hh"
#include "Achilles/Random.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

PionAbsorption::PionAbsorption(const YAML::Node &) {
    states[{PID::pionp(), PID::proton()}] = {{PID::neutron(), {PID::proton(), PID::neutron()}},
                                             {PID::neutron(), {PID::neutron(), PID::proton()}}};
}

std::vector<std::pair<PID, PID>> PionAbsorption::InitialStates() const {
    return {{PID::pion0(), PID::proton()},  {PID::pionp(), PID::proton()},
            {-PID::pionp(), PID::proton()}, {PID::pion0(), PID::neutron()},
            {PID::pionp(), PID::neutron()}, {-PID::pionp(), PID::neutron()}};
}

bool PionAbsorptionOneStep::AllowedAbsorption(Event &event, size_t part1, size_t part2) const {
    auto closest = FindClosest(event, part1, part2);
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
    const auto &particle2 = event.Hadrons()[part2];

    // Start shortest distance at well beyond nuclear diameter D ~ 10 fm = 1e-11 mm
    double shortest_p_distance = 6 / Constant::HBARC;
    double shortest_n_distance = 6 / Constant::HBARC;
    size_t closest_p_idx = SIZE_MAX, closest_n_idx = SIZE_MAX;

    // Particle 2 is the incoming nucleon
    for(std::size_t i = 0; i < event.Hadrons().size(); ++i) {
        if(i == part1 || i == part2) continue;
        if(!event.Hadrons()[i].IsBackground()) continue;
        auto distance = (event.Hadrons()[i].Position() - particle2.Position()).Magnitude2();
        if(distance < shortest_p_distance && event.Hadrons()[i].ID() == PID::proton()) {
            shortest_p_distance = distance;
            closest_p_idx = i;
        } else if(distance < shortest_n_distance && event.Hadrons()[i].ID() == PID::neutron()) {
            shortest_n_distance = distance;
            closest_n_idx = i;
        }
    }

    return {closest_p_idx, closest_n_idx};
}

InteractionResults PionAbsorption::CrossSection(Event &event, size_t part1, size_t part2) const {
    // const auto &particle1 = event.Hadrons()[part1];

    spdlog::debug("Did we get inside pion absorption cross section");
    const auto &particle2 = event.Hadrons()[part2];

    InteractionResults results;

    if(!AllowedAbsorption(event, part1, part2)) return results;

    absorption_partner.Status() = ParticleStatus::absorption_partner;

    double oset_abs_xsec = Oset_abs.AbsorptionCrossSection(event, part1, part2);

    // TODO
    // Insert logic for allowed outgoing states
    results.push_back({{particle2.ID(), absorption_partner.ID()}, oset_abs_xsec});
    return results;
}

std::vector<Particle> PionAbsorptionOneStep::GenerateMomentum(const Particle &particle1,
                                                              const Particle &particle2,
                                                              const std::vector<PID> &out_pids,
                                                              Random &random) const {
    // Boost to center of mass
    ThreeVector boostCM =
        (particle1.Momentum() + particle2.Momentum() + absorption_partner.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);

    // Generate outgoing momentum
    const double pcm = p1CM.Vec3().Magnitude();
    std::vector<double> rans(2);
    random.Generate(rans);

    double pR = pcm;
    double pTheta = acos(2 * rans[0] - 1);
    double pPhi = 2 * M_PI * rans[1];

    auto momentum = ThreeVector(ToCartesian({pR, pTheta, pPhi}));

    FourVector p1Out = FourVector(p1CM.E(), momentum[0], momentum[1], momentum[2]);
    FourVector p2Out = FourVector(p1CM.E(), -momentum[0], -momentum[1], -momentum[2]);

    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    spdlog::info("Absorption: Pin = {}, Pout = {} -> PSum = {}",
                 particle1.Momentum() + particle2.Momentum() + absorption_partner.Momentum(),
                 p1Out + p2Out,
                 particle1.Momentum() + particle2.Momentum() + absorption_partner.Momentum() -
                     p1Out - p2Out);

    // TODO
    // write momentum conservation unit test
    return {{out_pids[0], p1Out, particle1.Position()}, {out_pids[1], p2Out, particle2.Position()}};
}

bool PionAbsorptionTwoStep::AllowedAbsorption(Event &, size_t, size_t) const {
    return false;
}

std::vector<Particle> PionAbsorptionTwoStep::GenerateMomentum(const Particle &, const Particle &,
                                                              const std::vector<PID> &,
                                                              Random &) const {
    return {};
}
