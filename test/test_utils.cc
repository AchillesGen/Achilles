#include "catch2/catch.hpp"

#include "nuchic/Utilities.hh"

TEST_CASE("Coordiante Transforms", "[Utilities]") {
    SECTION("ToCartesian") {
        std::array<double, 3> rcoord{10, 0, 0}, xcoord{0, 0, 10};

        CHECK(nuchic::ToCartesian(rcoord) == xcoord);
    }
}

TEST_CASE("Brent", "[Utilities]") {
    auto func = [](const double &x) -> double {
        return (x+1)*(x-2);
    };

    SECTION("Find Roots") {
        nuchic::Brent brent(func);

        CHECK(brent.CalcRoot(1, 3) == Approx(2));
        CHECK(brent.CalcRoot(-2, 0) == Approx(-1));
        CHECK_THROWS_AS(brent.CalcRoot(3, 4), std::domain_error);
    }
}

TEST_CASE("GridSpace Generation", "[Utilities]") {
    SECTION("Logspace") {
        auto points = nuchic::Logspace(0, 10, 11);
        std::vector<double> points2{1, 10, 1E2, 1E3, 1E4, 1E5, 1E6, 1E7, 1E8, 1E9, 1E10};
        CHECK(points == points2);
    }

    SECTION("Linspace") {
        auto points = nuchic::Linspace(0, 10, 11);
        std::vector<double> points2{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        CHECK(points == points2);
    }
}
