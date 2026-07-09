// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch.hpp"

#include "Achilles/MesonBaryonAmplitudes.hh"

TEST_CASE("Initialization", "[MBAmplitudes]") {
    MBAmplitudes mbamp;
    REQUIRE(mbamp.NChargeChannels() == 16);
    REQUIRE(mbamp.NMesonBaryonChannels() == 4);
}
