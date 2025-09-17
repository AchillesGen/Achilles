#include <iostream>
#include <map>

#include "spdlog/spdlog.h"

#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/NasaInteractions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

REGISTER_INTERACTION(NasaInteractions);

double NasaInteractions::CrossSection(const Particle &particle1, const Particle &particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // Generate outgoing momentum
    double s = (p1Lab + p2Lab).M2();
    double plab = sqrt(pow(s, 2) / (4.0 * pow(Constant::mN, 2)) - s);
    return CrossSectionLab(samePID, plab);
}

ThreeVector NasaInteractions::MakeMomentum(bool, const double &pcm,
                                           const std::array<double, 2> &rans) const {
    double pR = pcm;
    double pTheta = acos(2 * rans[0] - 1);
    double pPhi = 2 * M_PI * rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}
