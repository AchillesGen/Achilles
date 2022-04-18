#include "catch2/catch.hpp"

#include "Achilles/Utilities.hh"

TEST_CASE("Bit operations", "[Utilities]") {

    SECTION("NextPermutation") {
        unsigned int inp = 0b100;
        auto perm1 = achilles::NextPermutation(inp);
        auto perm2 = achilles::NextPermutation(perm1);
        CHECK(perm1 == 0b1000);
        CHECK(perm2 == 0b10000);
    }

    SECTION("SetBit") {
        unsigned int inp = 0b100;
        CHECK(achilles::SetBit(inp, 0) == false);    
        CHECK(achilles::SetBit(inp, 1) == false);    
        CHECK(achilles::SetBit(inp, 2) == true);    
    }

    SECTION("SetBits") {
        unsigned int inp = 0b1010;
        auto result = achilles::SetBits(inp, 4);
        CHECK(result.size() == 2);
        CHECK(result[0] == 2);
        CHECK(result[1] == 8);
    }

    SECTION("IsPower2") {
        unsigned int inp = 0b100;
        CHECK(achilles::IsPower2(inp) == true);
        CHECK(achilles::IsPower2(inp+1) == false);
    }

    SECTION("Log2") {
        for(unsigned int i = 0; i < 16; ++i) {
            CHECK(achilles::Log2(1 << i) == i);
        }
    }
}

TEST_CASE("Tokenize", "[Utilities]") {
    std::string str = "This is a string";
    std::vector<std::string> tokens;
    achilles::tokenize(str, tokens);
    CHECK(tokens.size() == 4);
    CHECK(tokens[0] == "This");
    CHECK(tokens[1] == "is");
    CHECK(tokens[2] == "a");
    CHECK(tokens[3] == "string");
}

TEST_CASE("Coordiante Transforms", "[Utilities]") {
    SECTION("ToCartesian") {
        std::array<double, 3> rcoord{10, 0, 0}, xcoord{0, 0, 10};

        CHECK(achilles::ToCartesian(rcoord) == xcoord);
    }
}

TEST_CASE("Brent", "[Utilities]") {
    auto func = [](const double &x) -> double {
        return (x+1)*(x-2);
    };

    SECTION("Find Roots") {
        achilles::Brent brent(func);

        CHECK(brent.CalcRoot(1, 3) == Approx(2));
        CHECK(brent.CalcRoot(-2, 0) == Approx(-1));
        CHECK_THROWS_AS(brent.CalcRoot(3, 4), std::domain_error);
    }

    SECTION("Minimize") {
        achilles::Brent brent(func);
        CHECK(brent.Minimize(-2, 1) == Approx(0.5));
    }
}

TEST_CASE("GridSpace Generation", "[Utilities]") {
    SECTION("Logspace") {
        auto points = achilles::Logspace(0, 10, 11);
        std::vector<double> points2{1, 10, 1E2, 1E3, 1E4, 1E5, 1E6, 1E7, 1E8, 1E9, 1E10};
        CHECK(points == points2);
    }

    SECTION("Linspace") {
        auto points = achilles::Linspace(0, 10, 11);
        std::vector<double> points2{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        CHECK(points == points2);
    }
}
