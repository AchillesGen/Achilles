// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include <sstream>
#include <type_traits>

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "matchers.hh"

#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"
#include "spdlog/spdlog.h"

using achilles::FourPosition;
using achilles::FourVector;
using achilles::FourVectorT;
using achilles::ThreeBoost;
using achilles::ThreePosition;
using achilles::ThreeVector;
using achilles::ThreeVectorT;
namespace units = achilles::units;
using namespace achilles::units::literals;

namespace {
constexpr units::Dimensionless d(double x) {
    return units::Dimensionless{x};
}
constexpr units::Energy2 e2(double x) {
    return units::Energy2{x};
}
} // namespace

TEST_CASE("Vectors keep the legacy double[N] layout", "[Vectors]") {
    STATIC_REQUIRE(sizeof(ThreeVector) == 3 * sizeof(double));
    STATIC_REQUIRE(sizeof(FourVector) == 4 * sizeof(double));
    STATIC_REQUIRE(std::is_trivially_copyable_v<FourVector>);
    STATIC_REQUIRE(std::is_standard_layout_v<FourVector>);
}

TEST_CASE("Three Vector is constructed properly", "[Vectors]") {
    SECTION("Direct Constructors") {
        ThreeVector p1, p2(1_MeV, 2_MeV, 3_MeV), p3 = ThreeVector::FromNative({1, 2, 3});
        CHECK((p1[0] == 0_MeV && p1[1] == 0_MeV && p1[2] == 0_MeV));
        CHECK(p2 == p3);
    }

    SECTION("Move and Copy Constructors and Assignment") {
        ThreeVector p1(1_MeV, 2_MeV, 3_MeV);
        ThreeVector p2(p1);
        CHECK(p1 == p2);

        ThreeVector p3(std::move(p1));
        CHECK(p2 == p3);

        ThreeVector p4 = p3;
        ThreeVector p5 = std::move(p2);
        CHECK(p4 == p5);
    }
}

TEST_CASE("Four Vector is constructed properly", "[Vectors]") {
    SECTION("Direct Constructors") {
        ThreeVector p0(1_MeV, 2_MeV, 3_MeV);
        FourVector p1, p2(4_MeV, 1_MeV, 2_MeV, 3_MeV), p3 = FourVector::FromNative({4, 1, 2, 3});
        FourVector p4(p0, 4_MeV);
        CHECK((p1[1] == 0_MeV && p1[2] == 0_MeV && p1[3] == 0_MeV && p1[0] == 0_MeV));
        CHECK(p2 == p3);
        CHECK(p2 == p4);
    }

    SECTION("Move and Copy Constructors and Assignment") {
        FourVector p1(1_MeV, 2_MeV, 3_MeV, 4_MeV);
        FourVector p2(p1);
        CHECK(p1 == p2);

        FourVector p3(std::move(p1));
        CHECK(p2 == p3);

        FourVector p4 = p3;
        FourVector p5 = std::move(p2);
        CHECK(p4 == p5);
    }

    SECTION("Setters and Getters") {
        FourVector p1;
        p1.SetVectM(ThreeVector(1_MeV, 2_MeV, 3_MeV), 0_MeV);

        CHECK(p1.Px() == 1_MeV);
        CHECK(p1.Py() == 2_MeV);
        CHECK(p1.Pz() == 3_MeV);
        CHECK(p1.E().in(units::MeV) == sqrt(1 + 4 + 9));
        CHECK(p1.M() == 0_MeV);
    }
}

TEST_CASE("Accessors work as expected", "[Vectors]") {
    SECTION("Three Vector access") {
        ThreeVector p(1_MeV, 2_MeV, 3_MeV);
        CHECK(p[0] == 1_MeV);
        CHECK(p[1] == 2_MeV);
        CHECK(p[2] == 3_MeV);
        CHECK(p[0] == p.at(0));
        CHECK(p[1] == p.at(1));
        CHECK(p[2] == p.at(2));
        CHECK_THROWS_AS(p.at(3), std::range_error);
    }

    SECTION("Four Vector access") {
        FourVector p(1_MeV, 2_MeV, 3_MeV, 4_MeV);
        CHECK(p[0] == 1_MeV);
        CHECK(p[1] == 2_MeV);
        CHECK(p[2] == 3_MeV);
        CHECK(p[3] == 4_MeV);
        CHECK(p[0] == p.at(0));
        CHECK(p[1] == p.at(1));
        CHECK(p[2] == p.at(2));
        CHECK(p[3] == p.at(3));
        CHECK_THROWS_AS(p.at(4), std::range_error);
    }
}

TEST_CASE("Three Vector Overloaded Operators work as expected", "[Vectors]") {
    SECTION("Addition") {
        ThreeVector p1(1_MeV, 2_MeV, 3_MeV), p2(1_MeV, 2_MeV, 3_MeV), p3(2_MeV, 4_MeV, 6_MeV);
        CHECK(p1 + p2 == p3);

        p1 += p2;
        CHECK(p1 == p3);
    }

    SECTION("Subtraction and Negation") {
        ThreeVector p1(1_MeV, 2_MeV, 3_MeV), p2(1_MeV, 2_MeV, 3_MeV), p3(0_MeV, 0_MeV, 0_MeV),
            p4(-1_MeV, -2_MeV, -3_MeV);
        CHECK(p1 - p2 == p3);
        CHECK(-p1 == p4);

        p1 -= p2;
        CHECK(p1 == p3);
    }

    SECTION("Multiplication") {
        ThreeVector p1(1_MeV, 2_MeV, 3_MeV), p2(3_MeV, 6_MeV, 9_MeV);
        constexpr double scalar1 = 3;

        CHECK(p1 * p1 == e2(14));
        CHECK(p1 * p1 == p1.Dot(p1));
        CHECK(scalar1 * p1 == p1 * scalar1);
        CHECK(scalar1 * p1 == p2);

        p1 *= scalar1;
        CHECK(p1 == p2);
    }

    SECTION("Division") {
        ThreeVector p1(1_MeV, 2_MeV, 3_MeV), p2(3_MeV, 6_MeV, 9_MeV);
        constexpr double scalar = 3.0;

        CHECK(p2 / scalar == p1);

        p2 /= scalar;
        CHECK(p2 == p1);
    }
}

TEST_CASE("Three Vector Functions work as expected", "[Vectors]") {
    SECTION("Magnitude") {
        ThreeVector p1(4_MeV, 3_MeV, 2_MeV);
        constexpr double magnitude = 29;

        CHECK(p1.Magnitude2() == e2(magnitude));
        CHECK(p1.P2() == e2(magnitude));
        CHECK(p1.Magnitude().in(units::MeV) == sqrt(magnitude));
        CHECK(p1.P().in(units::MeV) == sqrt(magnitude));
    }

    SECTION("Transverse momentum and Angles") {
        ThreeVector p1(1_MeV, 2_MeV, 3_MeV), p2(0_MeV, 0_MeV, 4_MeV), p3(1_MeV, 1_MeV, 1_MeV);
        constexpr double pt2 = 5;

        CHECK(p1.Pt2() == e2(pt2));
        CHECK(p1.Pt().in(units::MeV) == sqrt(pt2));
        CHECK(p2.Theta() == 0);
        CHECK(p2.Phi() == 0);
        CHECK(p3.Phi() == M_PI_4);
    }

    SECTION("Cross Product") {
        // The cross product of two momenta has mass dimension 2.
        ThreeVector p1(3_MeV, 2_MeV, 1_MeV), p2(1_MeV, 2_MeV, 3_MeV);
        ThreeVectorT<2> p3(e2(4), e2(-8), e2(4));

        CHECK(p1.Cross(p2) == p3);
        CHECK(p1.Cross(p2) == -p2.Cross(p1));
    }

    SECTION("Unit Vector") {
        ThreeVector p1(3_MeV, 3_MeV, 3_MeV);
        ThreeBoost p2(d(1.0 / sqrt(3)), d(1.0 / sqrt(3)), d(1.0 / sqrt(3)));
        auto p1Unit = p1.Unit();

        STATIC_REQUIRE(std::is_same_v<decltype(p1Unit), ThreeBoost>);
        CHECK_THAT(p1Unit, ThreeVectorWithinRel(p2, 1e-12));
        CHECK(p1 / p1.Magnitude() == p1.Unit());
    }
}

TEST_CASE("Four Vector Overloaded Operators work as expected", "[Vectors]") {
    SECTION("Addition") {
        FourVector p1(1_MeV, 2_MeV, 3_MeV, 4_MeV), p2(4_MeV, 3_MeV, 2_MeV, 1_MeV),
            p3(5_MeV, 5_MeV, 5_MeV, 5_MeV);
        CHECK(p1 + p2 == p3);

        p1 += p2;
        CHECK(p1 == p3);
    }

    SECTION("Subtraction and Negation") {
        FourVector p1(1_MeV, 2_MeV, 3_MeV, 4_MeV), p2(1_MeV, 2_MeV, 3_MeV, 4_MeV),
            p3(0_MeV, 0_MeV, 0_MeV, 0_MeV), p4(-1_MeV, -2_MeV, -3_MeV, -4_MeV);
        CHECK(p1 - p2 == p3);
        CHECK(-p1 == p4);

        p1 -= p2;
        CHECK(p1 == p3);
    }

    SECTION("Multiplication") {
        FourVector p1(4_MeV, 1_MeV, 2_MeV, 3_MeV), p2(12_MeV, 3_MeV, 6_MeV, 9_MeV);
        constexpr double scalar1 = 3;

        CHECK(p1 * p1 == e2(16 - 9 - 4 - 1));
        CHECK(p1 * p1 == p1.Dot(p1));
        CHECK(scalar1 * p1 == p1 * scalar1);
        CHECK(scalar1 * p1 == p2);

        p1 *= scalar1;
        CHECK(p1 == p2);
    }

    SECTION("Division") {
        FourVector p1(4_MeV, 1_MeV, 2_MeV, 3_MeV), p2(12_MeV, 3_MeV, 6_MeV, 9_MeV);
        constexpr double scalar = 3.0;

        CHECK(p2 / scalar == p1);

        p2 /= scalar;
        CHECK(p2 == p1);
    }
}

TEST_CASE("Four Vector Functions work as expected", "[Vectors]") {
    SECTION("Magnitude and Mass") {
        FourVector p(4_MeV, 1_MeV, 2_MeV, 3_MeV);
        constexpr double mass = 2;

        CHECK(p.Magnitude2() == e2(mass));
        CHECK(p.M2() == e2(mass));
        CHECK(p.Magnitude().in(units::MeV) == sqrt(mass));
        CHECK(p.M().in(units::MeV) == sqrt(mass));
    }

    SECTION("Momentum and Transverse Momentum") {
        FourVector p(4_MeV, 1_MeV, 2_MeV, 3_MeV);
        constexpr double pvec2 = 14;
        constexpr double pt2 = 5;

        CHECK(p.P2() == e2(pvec2));
        CHECK(p.P().in(units::MeV) == sqrt(pvec2));
        CHECK(p.Pt2() == e2(pt2));
        CHECK(p.Pt().in(units::MeV) == sqrt(pt2));
    }

    SECTION("M2 == E^2 - P^2 and P2 == Pt2 + Pz^2") {
        FourVector p(1200_MeV, 300_MeV, 120_MeV, 700_MeV);
        const double E = p.E().in(units::MeV), P = p.P().in(units::MeV);

        CHECK_THAT(p.M2().in(units::MeV2), Catch::Matchers::WithinRel(E * E - P * P, 1e-12));
        CHECK(p.P2() == p.Pt2() + p.Pz() * p.Pz());
    }

    SECTION("Angles") {
        FourVector p1(4_MeV, 0_MeV, 0_MeV, 4_MeV), p2(4_MeV, 1_MeV, 1_MeV, 1_MeV),
            p3(4_MeV, 1_MeV, -1_MeV, 1_MeV);

        CHECK(p1.Theta() == 0);
        CHECK(p1.Phi() == 0);
        CHECK(p2.Phi() == M_PI_4);
        CHECK(p3.Phi() == 7 * M_PI_4);
    }

    SECTION("Cross Product") {
        FourVector p1(4_MeV, 3_MeV, 2_MeV, 1_MeV), p2(4_MeV, 1_MeV, 2_MeV, 3_MeV);
        FourVectorT<2> p3(e2(0), e2(4), e2(-8), e2(4));

        CHECK(p1.Cross(p2) == p3);
        CHECK(p1.Cross(p2) == -p2.Cross(p1));
    }

    SECTION("Boost") {
        FourVector p1(4_MeV, 3_MeV, 2_MeV, 1_MeV), p2(25_MeV, 12_MeV, 2_MeV, -3_MeV);
        auto beta = p2.BoostVector();
        STATIC_REQUIRE(std::is_same_v<decltype(beta), ThreeBoost>);
        auto partway = p1.Boost(beta);
        auto p3 = p1.Boost(beta).Boost(-beta);
        auto p4 =
            p1.Boost(beta).Boost(-beta.Px().native(), -beta.Py().native(), -beta.Pz().native());

        CHECK(partway != p1);
        CHECK_THAT(p3, FourVectorWithinRel(p1, 1e-12));
        CHECK_THAT(p4, FourVectorWithinRel(p1, 1e-12));

        FourVector pMass(units::Energy{sqrt(100 + 4)}, 0_MeV, 0_MeV, 2_MeV);
        beta = pMass.BoostVector();
        CHECK_THAT(pMass.Boost(-beta).M().in(units::MeV), Catch::Matchers::WithinRel(10, 1e-12));
    }

    SECTION("Invariant mass is boost invariant") {
        FourVector p(1200_MeV, 300_MeV, 120_MeV, 700_MeV);
        FourVector rest = p.Boost(-p.BoostVector());

        CHECK_THAT(rest.M().in(units::MeV), Catch::Matchers::WithinRel(p.M().in(units::MeV), 1e-9));
        CHECK_THAT(rest.P().in(units::MeV), Catch::Matchers::WithinAbs(0.0, 1e-6));
        CHECK_THAT(rest.E().in(units::MeV), Catch::Matchers::WithinRel(p.M().in(units::MeV), 1e-9));
    }

    SECTION("Rapidity") {
        FourVector p1(4_MeV, 1_MeV, 2_MeV, 3_MeV);
        constexpr double rapidity = 0.9729550745276566;

        CHECK_THAT(p1.Rapidity(), Catch::Matchers::WithinRel(rapidity, 1e-12));
    }

    SECTION("Angle between vectors") {
        FourVector p1(4_MeV, 1_MeV, 0_MeV, 3_MeV), p2(4_MeV, 3_MeV, 0_MeV, 1_MeV);
        double t1 = p1.Theta(), t2 = p2.Theta();

        CHECK_THAT(std::cos(t1 - t2), Catch::Matchers::WithinRel(p1.CosAngle(p2), 1e-12));
        CHECK_THAT(std::abs(t1 - t2), Catch::Matchers::WithinRel(p1.Angle(p2), 1e-12));
    }

    SECTION("DeltaR") {
        FourVector p1(4_MeV, 1_MeV, 2_MeV, 3_MeV), p2(4_MeV, 3_MeV, 2_MeV, 1_MeV);
        double DEta = p1.Rapidity() - p2.Rapidity();
        double DPhi = p1.Phi() - p2.Phi();
        double DR = sqrt(DEta * DEta + DPhi * DPhi);

        CHECK(p1.DeltaR(p2) == p2.DeltaR(p1));
        CHECK(p1.DeltaR(p2) == DR);
    }
}

TEST_CASE("Four-position uses the same template at D = -1", "[Vectors]") {
    FourPosition x(0.0_fm, 1.5_fm, 0.0_fm, 2.0_fm);

    CHECK_THAT(x.P().in(units::fm), Catch::Matchers::WithinRel(2.5, 1e-12));
    // Spacelike interval squared, read back as an area
    CHECK_THAT(x.M2().in(units::fm2), Catch::Matchers::WithinRel(-6.25, 1e-12));

    SECTION("native() on a position is MeV^-1, NOT fm") {
        CHECK_THAT(x.ToArray(units::fm)[3], Catch::Matchers::WithinRel(2.0, 1e-12));
        CHECK(x.Native()[3].native() < 0.02); // ~0.0101 MeV^-1
        CHECK_THAT(x.Native()[3].native(),
                   Catch::Matchers::WithinRel(2.0 / 197.32698045930246, 1e-12));
    }
}

TEST_CASE("Interop round-trips through explicit units", "[Vectors]") {
    FourVector p(1200_MeV, 300_MeV, 0_MeV, 700_MeV);
    FourVector p2 = FourVector::FromNative(p.ToArray(units::MeV));
    CHECK(p2 == p);

    ThreePosition x(0.0_fm, 1.5_fm, 2.0_fm);
    const auto arr = x.ToArray(units::fm);
    CHECK_THAT(arr[2], Catch::Matchers::WithinRel(2.0, 1e-12));
}

TEST_CASE("Rotations", "[Vectors]") {
    constexpr double eps = 1e-12;
    SECTION("Four Vectors") {
        FourVector p(4_MeV, 1_MeV, 2_MeV, 3_MeV);
        auto rotMat = p.AlignZ();
        auto result = p.Rotate(rotMat);
        FourVector expected(4_MeV, 0_MeV, 0_MeV, p.P());

        CHECK_THAT(result, FourVectorWithinRel(expected, eps));
    }

    SECTION("Three Vectors") {
        // Rotation axes are directions, i.e. dimensionless.
        ThreeBoost x(d(1), d(0), d(0)), y(d(0), d(1), d(0)), z(d(0), d(0), d(1));
        ThreeBoost expected(d(1.0 / sqrt(2.0)), d(1.0 / sqrt(2.0)), d(0));

        auto mat = x.Align(y);
        auto rotX = x.Rotate(mat);
        CHECK(rotX == y);

        rotX = x.Rotate(x.AlignZ());
        CHECK(rotX == z);

        constexpr std::array<double, 3> angles{0, 0, M_PI / 4};
        rotX = x.Rotate(angles);
        spdlog::info("rotX: {}", rotX);
        CHECK_THAT(rotX, ThreeVectorWithinRel(expected, eps));
    }
}

TEST_CASE("I/O and String", "[Vectors]") {
    SECTION("FourVector ToString and from string") {
        FourVector p(1_MeV, 2_MeV, 3_MeV, 4_MeV);
        FourVector p2;
        CHECK(p.ToString() == "FourVector(1.000000, 2.000000, 3.000000, 4.000000)");

        std::stringstream ss;
        ss << p;
        ss >> p2;

        CHECK(p == p2);
    }

    SECTION("ThreeVector ToString and from string") {
        ThreeVector p(1_MeV, 2_MeV, 3_MeV);
        ThreeVector p2;
        CHECK(p.ToString() == "ThreeVector(1.000000, 2.000000, 3.000000)");

        std::stringstream ss;
        ss << p;
        ss >> p2;

        CHECK(p == p2);
    }
}
