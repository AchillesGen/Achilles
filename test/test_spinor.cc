#include "catch2/catch.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/Spinor.hh"
#include <iostream>
#include <sstream>

using achilles::SpinMatrix;

TEST_CASE("Spinors", "[Spinors]") {
    SECTION("Spinor Inner product") {
        achilles::FourVector mom{1000, 0, 0, 100};
        double mass = 100;
        achilles::FourVector mom2{mass, 0, 0, 0};
        for(int i = -1; i < 2; i += 2) {
            for(int j = -1; j < 2; j += 2) {
                auto s1 = USpinor(i, mom);
                auto s2 = UBarSpinor(j, mom);
                auto s3 = USpinor(i, mom2);
                auto s4 = UBarSpinor(j, mom2);
                if(i == j)
                    CHECK((s2 * s1).real() == Approx(2.0 * mom.M()));
                else
                    CHECK(s2 * s1 == std::complex<double>());
                if(i == j)
                    CHECK((s4 * s3).real() == Approx(2.0 * mom2.M()));
                else
                    CHECK(s4 * s3 == std::complex<double>());
            }
        }
    }
}

TEST_CASE("GammaMatrix", "[Spinors]") {
    SECTION("Identities") {
        for(size_t i = 0; i < 4; ++i) {
            for(size_t j = 0; j < 4; ++j) {
                auto expected = i == j ? 2 * SpinMatrix::Identity() : 0 * SpinMatrix::Identity();
                if(i > 0) expected = -expected;
                auto result = SpinMatrix::GammaMu(i) * SpinMatrix::GammaMu(j) +
                              SpinMatrix::GammaMu(j) * SpinMatrix::GammaMu(i);
                CHECK(result == expected);
            }
        }

        CHECK(std::complex<double>(0, 1) * SpinMatrix::Gamma_0() * SpinMatrix::Gamma_1() *
                  SpinMatrix::Gamma_2() * SpinMatrix::Gamma_3() ==
              SpinMatrix::Gamma_5());

        CHECK(SpinMatrix::PL() * SpinMatrix::PL() == SpinMatrix::PL());
        CHECK(SpinMatrix::PR() * SpinMatrix::PR() == SpinMatrix::PR());
        CHECK(SpinMatrix::PR() * SpinMatrix::PL() == SpinMatrix());
        CHECK(SpinMatrix::PL() * SpinMatrix::PR() == SpinMatrix());
    }

    SECTION("Spinor inner product (massless)") {
        achilles::FourVector mom{500, 0, 300, 400};
        {
            auto s1 = USpinor(1, mom);
            auto s2 = UBarSpinor(1, mom);
            for(size_t i = 0; i < 4; ++i) {
                CHECK((s2 * SpinMatrix::GammaMu(i) * s1).real() == Approx(2 * mom[i]));
            }
        }
        {
            auto s1 = USpinor(-1, mom);
            auto s2 = UBarSpinor(-1, mom);
            for(size_t i = 0; i < 4; ++i) {
                CHECK((s2 * SpinMatrix::GammaMu(i) * s1).real() == Approx(2 * mom[i]));
            }
        }
    }

    SECTION("Spinor outer product (massless)") {
        achilles::FourVector mom{1000, 0, 0, 1000};
        SpinMatrix result;
        for(int i = -1; i < 2; i += 2) {
            for(int j = -1; j < 2; j += 2) {
                auto s1 = USpinor(i, mom);
                auto s2 = UBarSpinor(j, mom);
                result += s1.outer(s2);
            }
        }

        auto expected = SpinMatrix::Slashed(mom);
        for(size_t i = 0; i < result.size(); ++i) {
            CHECK(result[i].real() == Approx(expected[i].real()));
            CHECK(result[i].imag() == Approx(expected[i].imag()));
        }
    }

    SECTION("Spinor outer product (massive)") {
        const double p = 1000;
        const double mass = 1000;
        const double energy = sqrt(mass * mass + p * p);
        achilles::FourVector mom{energy, 0, 0, p};
        SpinMatrix result;
        for(int i = -1; i < 2; i += 2) {
            for(int j = -1; j < 2; j += 2) {
                auto s1 = USpinor(i, mom);
                auto s2 = UBarSpinor(j, mom);
                result += s1.outer(s2);
            }
        }

        auto expected = SpinMatrix::Slashed(mom) + mass * SpinMatrix::Identity();
        for(size_t i = 0; i < result.size(); ++i) {
            CHECK(result[i].real() == Approx(expected[i].real()));
            CHECK(result[i].imag() == Approx(expected[i].imag()));
        }
    }

    SECTION("SigmaMuNu") {
        for(size_t mu = 0; mu < 4; ++mu) {
            for(size_t nu = 0; nu < 4; ++nu) {
                CHECK(SpinMatrix::SigmaMuNu(mu, nu) == -SpinMatrix::SigmaMuNu(nu, mu));
                CHECK(SpinMatrix::SigmaMuNu(mu, nu) ==
                      std::complex<double>(0, 1) / 2.0 *
                          (SpinMatrix::GammaMu(mu) * SpinMatrix::GammaMu(nu) -
                           SpinMatrix::GammaMu(nu) * SpinMatrix::GammaMu(mu)));
            }
        }
    }

    SECTION("SpinMatrix * Spinor") {
        achilles::FourVector mom{1000, 0, 0, 1000};
        for(int i = -1; i < 2; i += 2) {
            auto s1 = USpinor(i, mom);
            auto result = SpinMatrix::Identity() * s1;
            auto expected = USpinor(i, mom);
            CHECK(result[0] == expected[0]);
            CHECK(result[1] == expected[1]);
            CHECK(result[2] == expected[2]);
            CHECK(result[3] == expected[3]);
        }
    }
}
