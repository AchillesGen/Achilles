#include "catch2/catch.hpp"

#include "nuchic/Spinor.hh"
#include "nuchic/Constants.hh"
#include <sstream>
#include <iostream>

using nuchic::SpinMatrix;

TEST_CASE("Spinors", "[Spinors]") {
    SECTION("Spinor Inner product") {
        nuchic::FourVector mom{1000, 0, 0, 100};
        for(int i = -1; i < 2; i+=2) {
            for(int j = -1; j < 2; j+=2) {
                auto s1 = USpinor(i, mom);
                auto s2 = UBarSpinor(j, mom);
                if(i == j) CHECK((s2*s1).real() == Approx(2.0*mom.M()));
                else CHECK(s2*s1 == std::complex<double>());
            }
        }
    }
}

TEST_CASE("GammaMatrix", "[Spinors]") {
    SECTION("Identities") {
        for(size_t i = 0; i < 4; ++i) {
            for(size_t j = 0; j < 4; ++j) {
                auto expected = i == j ? 2*SpinMatrix::Identity() : 0*SpinMatrix::Identity();
                if(i > 0) expected = -expected;
                auto result = SpinMatrix::GammaMu(i)*SpinMatrix::GammaMu(j) 
                            + SpinMatrix::GammaMu(j)*SpinMatrix::GammaMu(i);
                CHECK(result == expected);
            }
        }

        CHECK(std::complex<double>(0, 1)*SpinMatrix::Gamma_0()*SpinMatrix::Gamma_1()*SpinMatrix::Gamma_2()*SpinMatrix::Gamma_3() == SpinMatrix::Gamma_5());

        CHECK(SpinMatrix::PL()*SpinMatrix::PL() == SpinMatrix::PL());
        CHECK(SpinMatrix::PR()*SpinMatrix::PR() == SpinMatrix::PR());
        CHECK(SpinMatrix::PR()*SpinMatrix::PL() == SpinMatrix());
        CHECK(SpinMatrix::PL()*SpinMatrix::PR() == SpinMatrix());
    }

    SECTION("Spinor outer product") {
        nuchic::FourVector mom{1000, 0, 0, 1000};
        SpinMatrix result;
        for(int i = -1; i < 2; i+=2) {
            for(int j = -1; j < 2; j+=2) {
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

    SECTION("SigmaMuNu") {
        for(size_t mu = 0; mu < 4; ++mu) {
            for(size_t nu = 0; nu < 4; ++nu) {
                CHECK(SpinMatrix::SigmaMuNu(mu, nu) == -SpinMatrix::SigmaMuNu(nu, mu));
            }
        }
    }
}
