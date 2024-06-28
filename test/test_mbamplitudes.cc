#include "catch2/catch.hpp"

#include "Achilles/MesonBaryonAmplitudes.hh"

TEST_CASE("Initialization", "[MBAmplitudes]") {
    MBAmplitudes mbamp;
    REQUIRE(mbamp.NChargeChannels() == 16);
    REQUIRE(mbamp.NMesonBaryonChannels() == 4);
}
