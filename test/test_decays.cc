#include "Achilles/Particle.hh"
#include "catch2/catch.hpp"

#include "Achilles/DecayHandler.hh"
#include "catch_utils.hh"

using achilles::PID;

TEST_CASE("Parsing database", "[DecayHandler]") {
    achilles::DecayHandler decay("data/decays.txt");

    auto decays = decay.AllowedDecays(PID::deltapp());
    double br = 0;
    for(const auto &mode : decays) { br += mode.branching_ratio; }

    CHECK_THAT(br, Catch::Matchers::WithinAbs(1, 1e-8));
}

TEST_CASE("Two body decay", "[DecayHandler]") {
    achilles::DecayHandler decay("data/decays.txt");

    auto mom =
        GENERATE(take(30, randomMomentum(1000, achilles::ParticleInfo(PID::deltapp()).Mass())));
    achilles::Particle delta{PID::deltapp(), mom};

    auto outgoing = decay.Decay(delta);

    CHECK(outgoing.size() == 2);
    auto mom2 = outgoing[0].Momentum() + outgoing[1].Momentum();
    CHECK_THAT(mom2[0], Catch::Matchers::WithinAbs(mom[0], 1e-8));
    CHECK_THAT(mom2[1], Catch::Matchers::WithinAbs(mom[1], 1e-8));
    CHECK_THAT(mom2[2], Catch::Matchers::WithinAbs(mom[2], 1e-8));
    CHECK_THAT(mom2[3], Catch::Matchers::WithinAbs(mom[3], 1e-8));
}
