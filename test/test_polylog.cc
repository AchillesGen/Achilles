#include "catch2/catch.hpp"

#include "nuchic/Utilities.hh"
#include "spdlog/spdlog.h"

TEST_CASE("Zeta and PolyLog", "[Utilities]") {
    SECTION("Zeta function results are as expected") {
        CHECK(nuchic::zeta(2) == Approx(M_PI*M_PI/6));
        CHECK(nuchic::zeta(4) == Approx(pow(M_PI, 4)/90));
        CHECK(nuchic::zeta(12) == Approx(691*pow(M_PI, 12)/638512875));
        CHECK(nuchic::zeta(14) == Approx(2*pow(M_PI, 14)/18243225));
    }

    SECTION("Polylog function results are as expected") {
        CHECK(nuchic::PolyLog(-1, 0.5) == Approx(2));
        CHECK(nuchic::PolyLog(0, 0.5) == Approx(1));
        CHECK(nuchic::PolyLog(2, 0.5) == Approx(0.5822405264650127));
        CHECK(nuchic::PolyLog(7, 0.5) == Approx(0.5020145633247085));
        CHECK(nuchic::PolyLog(15, 0.5) == Approx(0.5000076381652628));
        CHECK(nuchic::PolyLog(3, 0.995) == Approx(1.193896987191249));
        CHECK(nuchic::PolyLog(3, 0.5) == Approx(0.5372131936080403));

        double y = 0.5;
        CHECK(nuchic::PolyLog(3, 1.0-y) == Approx(nuchic::PolyLog(3, y)));
    }
}
