#include "Achilles/CascadeInteractions/NasaInteractions.hh"
#include "Achilles/Event.hh"
#include "Achilles/Random.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

InteractionResults NasaInteraction::CrossSection(Event &event, size_t part1, size_t part2) const {
    const auto &particle1 = event.Hadrons()[part1];
    const auto &particle2 = event.Hadrons()[part2];

    bool samePID = particle1.ID() == particle2.ID();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // Generate outgoing momentum
    double s = (p1Lab + p2Lab).M2();
    double smin = pow(particle1.Mass() + particle2.Mass(), 2);
    double plab = sqrt(pow(s, 2) / smin - s);
    return {{{particle1.ID(), particle2.ID()}, CrossSectionLab(samePID, plab)}};
}

std::vector<std::pair<PID, PID>> NasaInteraction::InitialStates() const {
    return {{PID::proton(), PID::proton()},
            {PID::neutron(), PID::proton()},
            {PID::neutron(), PID::neutron()}};
}

std::vector<Particle> NasaInteraction::GenerateMomentum(const Particle &particle1,
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

ThreeVector NasaInteraction::MakeMomentum(bool, double pcm, const std::vector<double> &rans) const {
    double pR = pcm;
    double pTheta = acos(2 * rans[0] - 1);
    double pPhi = 2 * M_PI * rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}
