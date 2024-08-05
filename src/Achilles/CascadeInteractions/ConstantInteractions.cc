#include "Achilles/CascadeInteractions/ConstantInteractions.hh"
#include "Achilles/Event.hh"
#include "Achilles/Random.hh"
#include "Achilles/Utilities.hh"
#include <iostream>

using namespace achilles;

ConstantInteraction::ConstantInteraction(const YAML::Node &node) {
    for(const auto &interactions : node["InitialStates"]) {
        auto incoming = std::make_pair<PID, PID>(interactions["Incoming"][0].as<PID>(),
                                                 interactions["Incoming"][1].as<PID>());
        AddInteraction(incoming, interactions["Outgoing"].as<InteractionResults>());
    }
}

InteractionResults ConstantInteraction::CrossSection(Event &event, size_t part1,
                                                     size_t part2) const {
    const auto &p1 = event.Hadrons()[part1];
    const auto &p2 = event.Hadrons()[part2];

    try {
        return m_interactions.at({p1.ID(), p2.ID()});
    } catch(std::out_of_range &e) {
        std::cout << p1 << " " << p2 << std::endl;
        throw;
    }
}

// TODO: Implement generic n-body final state phase space
std::vector<Particle> ConstantInteraction::GenerateMomentum(const Particle &particle1,
                                                            const Particle &particle2,
                                                            const std::vector<PID> &out_pids,
                                                            Random &random) const {
    // Boost to center of mass
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);

    // Generate outgoing momentum
    bool samePID = particle1.ID() == particle2.ID();
    const double pcm = p1CM.Vec3().Magnitude();
    std::vector<double> rans(2);
    random.Generate(rans);
    ThreeVector momentum = MakeMomentum(samePID, pcm, rans);

    FourVector p1Out = FourVector(p1CM.E(), momentum[0], momentum[1], momentum[2]);
    FourVector p2Out = FourVector(p1CM.E(), -momentum[0], -momentum[1], -momentum[2]);

    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    return {{out_pids[0], p1Out, particle1.Position()}, {out_pids[1], p2Out, particle2.Position()}};
}

ThreeVector ConstantInteraction::MakeMomentum(bool, double pcm,
                                              const std::vector<double> &rans) const {
    double pR = pcm;
    double pTheta = acos(2 * rans[0] - 1);
    double pPhi = 2 * M_PI * rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}
