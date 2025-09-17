#include <sstream>

#include "catch2/catch.hpp"

#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"
#include "spdlog/spdlog.h"

TEST_CASE("Three Vector is constructed properly", "[Vectors]") {
    SECTION("Direct Constructors") {
        achilles::ThreeVector p1, p2(1, 2, 3), p3(std::array<double, 3>{1, 2, 3});
        CHECK((p1[0] == 0 && p1[1] == 0 && p1[2] == 0));
        CHECK(p2 == p3);
    }

    SECTION("Move and Copy Constructors and Assignment") {
        achilles::ThreeVector p1(1, 2, 3);
        achilles::ThreeVector p2(p1);
        CHECK(p1 == p2);

        achilles::ThreeVector p3(std::move(p1));
        CHECK(p2 == p3);

        achilles::ThreeVector p4 = p3;
        achilles::ThreeVector p5 = std::move(p2);
        CHECK(p4 == p5);
    }
}

TEST_CASE("Four Vector is constructed properly", "[Vectors]") {
    SECTION("Direct Constructors") {
        achilles::ThreeVector p0(1, 2, 3);
        achilles::FourVector p1, p2(4, 1, 2, 3), p3(std::array<double, 4>{4, 1, 2, 3});
        achilles::FourVector p4(p0, 4);
        CHECK((p1[1] == 0 && p1[2] == 0 && p1[3] == 0 && p1[0] == 0));
        CHECK(p2 == p3);
        CHECK(p2 == p4);
    }

    SECTION("Move and Copy Constructors and Assignment") {
        achilles::FourVector p1(1, 2, 3, 4);
        achilles::FourVector p2(p1);
        CHECK(p1 == p2);

        achilles::FourVector p3(std::move(p1));
        CHECK(p2 == p3);

        achilles::FourVector p4 = p3;
        achilles::FourVector p5 = std::move(p2);
        CHECK(p4 == p5);
    }

    SECTION("Setters and Getters") {
        achilles::FourVector p1;
        p1.SetVectM(achilles::ThreeVector(1, 2, 3), 0);

        CHECK(p1.Px() == 1);
        CHECK(p1.Py() == 2);
        CHECK(p1.Pz() == 3);
        CHECK(p1.E() == sqrt(1 + 4 + 9));
        CHECK(p1.M() == 0);
    }
}

TEST_CASE("Accessors work as expected", "[Vectors]") {
    SECTION("Three Vector access") {
        achilles::ThreeVector p(1, 2, 3);
        CHECK(p[0] == 1);
        CHECK(p[1] == 2);
        CHECK(p[2] == 3);
        CHECK(p[0] == p.at(0));
        CHECK(p[1] == p.at(1));
        CHECK(p[2] == p.at(2));
        CHECK_THROWS_AS(p.at(3), std::range_error);
    }

    SECTION("Four Vector access") {
        achilles::FourVector p(1, 2, 3, 4);
        CHECK(p[0] == 1);
        CHECK(p[1] == 2);
        CHECK(p[2] == 3);
        CHECK(p[3] == 4);
        CHECK(p[0] == p.at(0));
        CHECK(p[1] == p.at(1));
        CHECK(p[2] == p.at(2));
        CHECK(p[3] == p.at(3));
        CHECK_THROWS_AS(p.at(4), std::range_error);
    }
}

TEST_CASE("Three Vector Overloaded Operators work as expected", "[Vectors]") {
    SECTION("Addition") {
        achilles::ThreeVector p1(1, 2, 3), p2(1, 2, 3), p3(2, 4, 6);
        CHECK(p1 + p2 == p3);

        p1 += p2;
        CHECK(p1 == p3);
    }

    SECTION("Subtraction and Negation") {
        achilles::ThreeVector p1(1, 2, 3), p2(1, 2, 3), p3(0, 0, 0), p4(-1, -2, -3);
        CHECK(p1 - p2 == p3);
        CHECK(-p1 == p4);

        p1 -= p2;
        CHECK(p1 == p3);
    }

    SECTION("Multiplication") {
        achilles::ThreeVector p1(1, 2, 3), p2(3, 6, 9);
        constexpr double scalar1 = 3, scalar2 = 14;

        CHECK(p1 * p1 == scalar2);
        CHECK(p1 * p1 == p1.Dot(p1));
        CHECK(scalar1 * p1 == p1 * scalar1);
        CHECK(scalar1 * p1 == p2);

        p1 *= scalar1;
        CHECK(p1 == p2);
    }

    SECTION("Division") {
        achilles::ThreeVector p1(1, 2, 3), p2(3, 6, 9);
        constexpr double scalar = 3.0;

        CHECK(p2 / scalar == p1);

        p2 /= scalar;
        CHECK(p2 == p1);
    }
}

TEST_CASE("Three Vector Functions work as expected", "[Vectors]") {
    SECTION("Magnitude") {
        achilles::ThreeVector p1(4, 3, 2);
        constexpr double magnitude = 29;

        CHECK(p1.Magnitude2() == magnitude);
        CHECK(p1.P2() == magnitude);
        CHECK(p1.Magnitude() == sqrt(magnitude));
        CHECK(p1.P() == sqrt(magnitude));
    }

    SECTION("Transverse momentum and Angles") {
        achilles::ThreeVector p1(1, 2, 3), p2(0, 0, 4), p3(1, 1, 1);
        constexpr double pt2 = 5;

        CHECK(p1.Pt2() == pt2);
        CHECK(p1.Pt() == sqrt(pt2));
        CHECK(p2.Theta() == 0);
        CHECK(p2.Phi() == 0);
        CHECK(p3.Phi() == M_PI_4);
    }

    SECTION("Cross Product") {
        achilles::ThreeVector p1(3, 2, 1), p2(1, 2, 3), p3(4, -8, 4);

        CHECK(p1.Cross(p2) == p3);
        CHECK(p1.Cross(p2) == -p2.Cross(p1));
    }

    SECTION("Unit Vector") {
        achilles::ThreeVector p1(3, 3, 3), p2(1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3));
        auto p1Unit = p1.Unit();

        CHECK((p1Unit[0] == Approx(p2[0]) && p1Unit[1] == Approx(p2[1]) &&
               p1Unit[2] == Approx(p2[2])));
        CHECK(p1 / p1.Magnitude() == p1.Unit());
    }
}

TEST_CASE("Four Vector Overloaded Operators work as expected", "[Vectors]") {
    SECTION("Addition") {
        achilles::FourVector p1(1, 2, 3, 4), p2(4, 3, 2, 1), p3(5, 5, 5, 5);
        CHECK(p1 + p2 == p3);

        p1 += p2;
        CHECK(p1 == p3);
    }

    SECTION("Subtraction and Negation") {
        achilles::FourVector p1(1, 2, 3, 4), p2(1, 2, 3, 4), p3(0, 0, 0, 0), p4(-1, -2, -3, -4);
        CHECK(p1 - p2 == p3);
        CHECK(-p1 == p4);

        p1 -= p2;
        CHECK(p1 == p3);
    }

    SECTION("Multiplication") {
        achilles::FourVector p1(4, 1, 2, 3), p2(12, 3, 6, 9);
        constexpr double scalar1 = 3, scalar2 = 16 - 9 - 4 - 1;

        CHECK(p1 * p1 == scalar2);
        CHECK(p1 * p1 == p1.Dot(p1));
        CHECK(scalar1 * p1 == p1 * scalar1);
        CHECK(scalar1 * p1 == p2);

        p1 *= scalar1;
        CHECK(p1 == p2);
    }

    SECTION("Division") {
        achilles::FourVector p1(4, 1, 2, 3), p2(12, 3, 6, 9);
        constexpr double scalar = 3.0;

        CHECK(p2 / scalar == p1);

        p2 /= scalar;
        CHECK(p2 == p1);
    }
}

TEST_CASE("Four Vector Functions work as expected", "[Vectors]") {
    SECTION("Magnitude and Mass") {
        achilles::FourVector p(4, 1, 2, 3);
        constexpr double mass = 2;

        CHECK(p.Magnitude2() == mass);
        CHECK(p.M2() == mass);
        CHECK(p.Magnitude() == sqrt(mass));
        CHECK(p.M() == sqrt(mass));
    }

    SECTION("Momentum and Transverse Momentum") {
        achilles::FourVector p(4, 1, 2, 3);
        constexpr double pvec2 = 14;
        constexpr double pt2 = 5;

        CHECK(p.P2() == pvec2);
        CHECK(p.P() == sqrt(pvec2));
        CHECK(p.Pt2() == pt2);
        CHECK(p.Pt() == sqrt(pt2));
    }

    SECTION("Angles") {
        achilles::FourVector p1(4, 0, 0, 4), p2(4, 1, 1, 1), p3(4, 1, -1, 1);

        CHECK(p1.Theta() == 0);
        CHECK(p1.Phi() == 0);
        CHECK(p2.Phi() == M_PI_4);
        CHECK(p3.Phi() == 7 * M_PI_4);
    }

    SECTION("Cross Product") {
        achilles::FourVector p1(4, 3, 2, 1), p2(4, 1, 2, 3), p3(0, 4, -8, 4);

        CHECK(p1.Cross(p2) == p3);
        CHECK(p1.Cross(p2) == -p2.Cross(p1));
    }

    SECTION("Boost") {
        achilles::FourVector p1(4, 3, 2, 1), p2(25, 12, 2, -3);
        auto beta = p2.BoostVector();
        auto partway = p1.Boost(beta);
        auto p3 = p1.Boost(beta).Boost(-beta);
        auto p4 = p1.Boost(beta).Boost(-beta.Px(), -beta.Py(), -beta.Pz());

        CHECK(partway != p1);
        CHECK((p3[0] == Approx(p1[0]) && p3[1] == Approx(p1[1]) && p3[2] == Approx(p1[2]) &&
               p3[3] == Approx(p1[3])));
        CHECK((p4[0] == Approx(p1[0]) && p4[1] == Approx(p1[1]) && p4[2] == Approx(p1[2]) &&
               p4[3] == Approx(p1[3])));

        achilles::FourVector pMass(sqrt(100 + 4), 0, 0, 2);
        beta = pMass.BoostVector();
        CHECK(pMass.Boost(-beta).M() == Approx(10));
    }

    SECTION("Rapidity") {
        achilles::FourVector p1(4, 1, 2, 3);
        constexpr double rapidity = 0.9729550745276566;

        CHECK(p1.Rapidity() == Approx(rapidity));
    }

    SECTION("Angle between vectors") {
        achilles::FourVector p1(4, 1, 0, 3), p2(4, 3, 0, 1);
        double t1 = p1.Theta(), t2 = p2.Theta();

        CHECK(std::cos(t1 - t2) == Approx(p1.CosAngle(p2)));
        CHECK(std::abs(t1 - t2) == Approx(p1.Angle(p2)));
    }

    SECTION("DeltaR") {
        achilles::FourVector p1(4, 1, 2, 3), p2(4, 3, 2, 1);
        double DEta = p1.Rapidity() - p2.Rapidity();
        double DPhi = p1.Phi() - p2.Phi();
        double DR = sqrt(DEta * DEta + DPhi * DPhi);

        CHECK(p1.DeltaR(p2) == p2.DeltaR(p1));
        CHECK(p1.DeltaR(p2) == DR);
    }
}

TEST_CASE("Rotations", "[Vectors]") {
    constexpr double eps = 1e-12;
    SECTION("Four Vectors") {
        achilles::FourVector p(4, 1, 2, 3);
        auto rotMat = p.AlignZ();
        auto result = p.Rotate(rotMat);

        CHECK(result[0] == p[0]);
        CHECK(result[1] == Approx(0).margin(eps));
        CHECK(result[2] == Approx(0).margin(eps));
        CHECK(result[3] == Approx(p.P()).margin(eps));
    }

    SECTION("Three Vectors") {
        achilles::ThreeVector x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);

        auto mat = x.Align(y);
        auto rotX = x.Rotate(mat);
        CHECK(rotX == y);

        rotX = x.Rotate(x.AlignZ());
        CHECK(rotX == z);

        constexpr std::array<double, 3> angles{0, 0, M_PI / 4};
        rotX = x.Rotate(angles);
        spdlog::info("rotX: {}", rotX);
        CHECK(rotX[0] == Approx(1.0 / sqrt(2.0)).margin(eps));
        CHECK(rotX[1] == Approx(1.0 / sqrt(2.0)).margin(eps));
        CHECK(rotX[2] == Approx(0.0).margin(eps));
    }
}

TEST_CASE("I/O and String", "[Vectors]") {
    SECTION("FourVector ToString and from string") {
        achilles::FourVector p(1, 2, 3, 4);
        achilles::FourVector p2;
        CHECK(p.ToString() == "FourVector(1.000000, 2.000000, 3.000000, 4.000000)");

        std::stringstream ss;
        ss << p;
        ss >> p2;

        CHECK(p == p2);
    }

    SECTION("ThreeVector ToString and from string") {
        achilles::ThreeVector p(1, 2, 3);
        achilles::ThreeVector p2;
        CHECK(p.ToString() == "ThreeVector(1.000000, 2.000000, 3.000000)");

        std::stringstream ss;
        ss << p;
        ss >> p2;

        CHECK(p == p2);
    }
}
