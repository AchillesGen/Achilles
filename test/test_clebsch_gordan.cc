#include "catch2/catch.hpp"

#include "Achilles/ClebschGordan.hh"

using ZSpins = std::vector<std::pair<double, double>>;

void CompareAll(double j1, double j2, const ZSpins &zspins, const std::vector<double> &total_spin,
                const std::vector<double> expected) {
    size_t idx = 0;
    for(const auto &[m1, m2] : zspins) {
        for(const auto &j3 : total_spin) {
            achilles::SpinState s1{j1, m1};
            achilles::SpinState s2{j2, m2};
            achilles::SpinState s3{j3, m1 + m2};

            CHECK_THAT(achilles::ClebschGordan(s1, s2, s3),
                       Catch::Matchers::WithinAbs(expected[idx++], 1e-8));
        }
    }
}

TEST_CASE("Integer Spins", "[ClebschGordan]") {
    const double sqrt2 = sqrt(2.0);
    const double sqrt3 = sqrt(3.0);
    const double sqrt5 = sqrt(5.0);
    SECTION("j1=1, j2=1") {
        SECTION("m3=2") {
            achilles::SpinState s1{1, 1};
            achilles::SpinState s2{1, 1};
            achilles::SpinState s3{2, 1};

            double expected = 1.0;
            CHECK_THAT(achilles::ClebschGordan(s1, s2, s3),
                       Catch::Matchers::WithinAbs(expected, 1e-8));
        }

        SECTION("m3=1") {
            ZSpins zspins = {{1, 0}, {0, 1}};
            std::vector<double> total_spin = {2, 1};
            std::vector<double> expected = {1 / sqrt2, 1 / sqrt2, 1 / sqrt2, -1 / sqrt2};
            CompareAll(1, 1, zspins, total_spin, expected);
        }

        SECTION("m3=0") {
            ZSpins zspins = {{1, -1}, {0, 0}, {-1, 1}};
            std::vector<double> total_spin = {2, 1, 0};
            std::vector<double> expected = {1 / (sqrt2 * sqrt3), 1 / sqrt2,  1 / sqrt3,
                                            sqrt2 / sqrt3,       0,          -1 / sqrt3,
                                            1 / (sqrt2 * sqrt3), -1 / sqrt2, 1 / sqrt3};
            CompareAll(1, 1, zspins, total_spin, expected);
        }

        SECTION("m3=-1") {
            ZSpins zspins = {{-1, 0}, {0, -1}};
            std::vector<double> total_spin = {2, 1};
            std::vector<double> expected = {1 / sqrt2, -1 / sqrt2, 1 / sqrt2, 1 / sqrt2};
            CompareAll(1, 1, zspins, total_spin, expected);
        }
    }

    SECTION("j1=2, j1=1") {
        double j1 = 2;
        double j2 = 1;
        SECTION("m3=3") {
            achilles::SpinState s1{j1, 2};
            achilles::SpinState s2{j2, 1};
            achilles::SpinState s3{3, 3};

            double expected = 1.0;
            CHECK_THAT(achilles::ClebschGordan(s1, s2, s3),
                       Catch::Matchers::WithinAbs(expected, 1e-8));
        }

        SECTION("m3=2") {
            ZSpins zspins = {{2, 0}, {1, 1}};
            std::vector<double> total_spin = {3, 2};
            std::vector<double> expected = {1 / sqrt3, sqrt2 / sqrt3, sqrt2 / sqrt3, -1 / sqrt3};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }

        SECTION("m3=1") {
            ZSpins zspins = {{2, -1}, {1, 0}, {0, 1}};
            std::vector<double> total_spin = {3, 2, 1};
            std::vector<double> expected = {1 / (sqrt3 * sqrt5), 1 / sqrt3,
                                            sqrt3 / sqrt5,       2 * sqrt2 / (sqrt3 * sqrt5),
                                            1 / (sqrt2 * sqrt3), -sqrt3 / (sqrt2 * sqrt5),
                                            sqrt2 / sqrt5,       -1 / sqrt2,
                                            1 / (sqrt2 * sqrt5)};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }

        SECTION("m3=0") {
            ZSpins zspins = {{1, -1}, {0, 0}, {-1, 1}};
            std::vector<double> total_spin = {3, 2, 1};
            std::vector<double> expected = {1 / sqrt5,     1 / sqrt2,  sqrt3 / sqrt(10),
                                            sqrt3 / sqrt5, 0,          -sqrt2 / sqrt5,
                                            1 / sqrt5,     -1 / sqrt2, sqrt3 / sqrt(10)};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }
    }

    SECTION("j1=2, j1=2") {
        double j1 = 2;
        double j2 = 2;
        SECTION("m3=4") {
            achilles::SpinState s1{j1, 2};
            achilles::SpinState s2{j2, 2};
            achilles::SpinState s3{4, 4};

            double expected = 1.0;
            CHECK_THAT(achilles::ClebschGordan(s1, s2, s3),
                       Catch::Matchers::WithinAbs(expected, 1e-8));
        }

        SECTION("m3=3") {
            ZSpins zspins = {{2, 1}, {1, 2}};
            std::vector<double> total_spin = {4, 3};
            std::vector<double> expected = {1 / sqrt2, 1 / sqrt2, 1 / sqrt2, -1 / sqrt2};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }

        SECTION("m3=2") {
            ZSpins zspins = {{2, 0}, {1, 1}, {0, 2}};
            std::vector<double> total_spin = {4, 3, 2};
            std::vector<double> expected = {sqrt3 / sqrt(14), 1 / sqrt2,  sqrt2 / sqrt(7),
                                            2 / (sqrt(7)),    0,          -sqrt3 / sqrt(7),
                                            sqrt3 / sqrt(14), -1 / sqrt2, sqrt2 / sqrt(7)};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }

        SECTION("m3=1") {
            ZSpins zspins = {{2, -1}, {1, 0}, {0, 1}, {-1, 2}};
            std::vector<double> total_spin = {4, 3, 2, 1};
            std::vector<double> expected = {
                1 / sqrt(14),    sqrt3 / sqrt(10),  sqrt3 / sqrt(7), 1 / sqrt5,
                sqrt3 / sqrt(7), 1 / sqrt5,         -1 / sqrt(14),   -sqrt3 / sqrt(10),
                sqrt3 / sqrt(7), -1 / sqrt5,        -1 / sqrt(14),   sqrt3 / sqrt(10),
                1 / sqrt(14),    -sqrt3 / sqrt(10), sqrt3 / sqrt(7), -1 / sqrt(5)};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }

        SECTION("m3=0") {
            ZSpins zspins = {{2, -2}, {1, -1}, {0, 0}, {-1, 1}, {-2, 2}};
            std::vector<double> total_spin = {4, 3, 2, 1, 0};
            std::vector<double> expected = {1 / sqrt(70),     1 / sqrt(10),
                                            sqrt2 / sqrt(7),  sqrt2 / sqrt5,
                                            1 / sqrt5,        sqrt(8.0 / 35),
                                            sqrt2 / sqrt5,    1 / sqrt(14),
                                            -1 / sqrt(10),    -1 / sqrt5,
                                            sqrt(18.0 / 35),  0,
                                            -sqrt2 / sqrt(7), 0,
                                            1 / sqrt5,        sqrt(8.0 / 35),
                                            -sqrt2 / sqrt5,   1 / sqrt(14),
                                            1 / sqrt(10),     -1 / sqrt5,
                                            1 / sqrt(70),     -1 / sqrt(10),
                                            sqrt2 / sqrt(7),  -sqrt2 / sqrt5,
                                            1 / sqrt5};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }
    }
}

TEST_CASE("Half-Integer Spins", "[ClebschGordan]") {
    const double sqrt2 = sqrt(2.0);
    const double sqrt3 = sqrt(3.0);
    SECTION("j1=1/2, j2=1/2") {
        SECTION("m3=1") {
            achilles::SpinState s1{0.5, 0.5};
            achilles::SpinState s2{0.5, 0.5};
            achilles::SpinState s3{1, 1};

            double expected = 1.0;
            CHECK_THAT(achilles::ClebschGordan(s1, s2, s3),
                       Catch::Matchers::WithinAbs(expected, 1e-8));
        }

        SECTION("m3=-1") {
            achilles::SpinState s1{0.5, -0.5};
            achilles::SpinState s2{0.5, -0.5};
            achilles::SpinState s3{1, 1};

            double expected = 1.0;
            CHECK_THAT(achilles::ClebschGordan(s1, s2, s3),
                       Catch::Matchers::WithinAbs(expected, 1e-8));
        }

        SECTION("m3=0") {
            ZSpins zspins = {{0.5, -0.5}, {-0.5, 0.5}};
            std::vector<double> total_spin = {1, 0};
            std::vector<double> expected = {1 / sqrt2, 1 / sqrt2, 1 / sqrt2, -1 / sqrt2};
            CompareAll(0.5, 0.5, zspins, total_spin, expected);
        }
    }

    SECTION("j1=3/2, j1=1/2") {
        double j1 = 1.5;
        double j2 = 0.5;
        SECTION("m3=2") {
            achilles::SpinState s1{j1, 1.5};
            achilles::SpinState s2{j2, 0.5};
            achilles::SpinState s3{2, 2};

            double expected = 1.0;
            CHECK_THAT(achilles::ClebschGordan(s1, s2, s3),
                       Catch::Matchers::WithinAbs(expected, 1e-8));
        }

        SECTION("m3=1") {
            ZSpins zspins = {{1.5, -0.5}, {0.5, 0.5}};
            std::vector<double> total_spin = {2, 1};
            std::vector<double> expected = {0.5, sqrt3 / 2, sqrt3 / 2, -0.5};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }

        SECTION("m3=0") {
            ZSpins zspins = {{0.5, -0.5}, {-0.5, 0.5}};
            std::vector<double> total_spin = {2, 1};
            std::vector<double> expected = {1 / sqrt2, 1 / sqrt2, 1 / sqrt2, -1 / sqrt2};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }
    }

    SECTION("j1=3/2, j1=3/2") {
        double j1 = 1.5;
        double j2 = 1.5;
        SECTION("m3=3") {
            achilles::SpinState s1{j1, 1.5};
            achilles::SpinState s2{j2, 1.5};
            achilles::SpinState s3{3, 3};

            double expected = 1.0;
            CHECK_THAT(achilles::ClebschGordan(s1, s2, s3),
                       Catch::Matchers::WithinAbs(expected, 1e-8));
        }

        SECTION("m3=2") {
            ZSpins zspins = {{1.5, 0.5}, {0.5, 1.5}};
            std::vector<double> total_spin = {3, 2};
            std::vector<double> expected = {1 / sqrt2, 1 / sqrt2, 1 / sqrt2, -1 / sqrt2};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }

        SECTION("m3=1") {
            ZSpins zspins = {{1.5, -0.5}, {0.5, 0.5}, {-0.5, 1.5}};
            std::vector<double> total_spin = {3, 2, 1};
            std::vector<double> expected = {1 / sqrt(5),     1 / sqrt2,  sqrt3 / sqrt(10),
                                            sqrt3 / sqrt(5), 0,          -sqrt2 / sqrt(5),
                                            1 / sqrt(5),     -1 / sqrt2, sqrt3 / sqrt(10)};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }

        SECTION("m3=0") {
            ZSpins zspins = {{1.5, -1.5}, {0.5, -0.5}, {-0.5, 0.5}, {-1.5, 1.5}};
            std::vector<double> total_spin = {3, 2, 1, 0};
            std::vector<double> expected = {1 / sqrt(20),   0.5,  sqrt(9.0 / 20), 0.5,
                                            sqrt(9.0 / 20), 0.5,  -1 / sqrt(20),  -0.5,
                                            sqrt(9.0 / 20), -0.5, -1 / sqrt(20),  0.5,
                                            1 / sqrt(20),   -0.5, sqrt(9. / 20),  -0.5};
            CompareAll(j1, j2, zspins, total_spin, expected);
        }
    }
}
