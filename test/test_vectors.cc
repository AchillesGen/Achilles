#include "catch2/catch.hpp"

#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"

TEST_CASE("Three Vector is constructed properly", "[Vectors]") {
    SECTION("Direct Constructors") {
        nuchic::ThreeVector p1, p2(1, 2, 3), p3(std::array<double, 3>{1, 2, 3});
        CHECK((p1[0] == 0 && p1[1] == 0 && p1[2] == 0));
        CHECK(p2 == p3);
    }

    SECTION("Move and Copy Constructors and Assignment") {
        nuchic::ThreeVector p1(1, 2, 3);
        nuchic::ThreeVector p2(p1);
        CHECK(p1 == p2);

        nuchic::ThreeVector p3(std::move(p1));
        CHECK(p2 == p3);

        nuchic::ThreeVector p4 = p3;
        nuchic::ThreeVector p5 = std::move(p2);
        CHECK(p4 == p5);
    }
}

TEST_CASE("Four Vector is constructed properly", "[Vectors]") {
    SECTION("Direct Constructors") {
        nuchic::ThreeVector p0(1, 2, 3);
        nuchic::FourVector p1, p2(1, 2, 3, 4), p3(std::array<double, 4>{1, 2, 3, 4});
        nuchic::FourVector p4(p0, 4);
        CHECK((p1[0] == 0 && p1[1] == 0 && p1[2] == 0 && p1[3] == 0));
        CHECK(p2 == p3);
        CHECK(p2 == p4);
    }

    SECTION("Move and Copy Constructors and Assignment") {
        nuchic::FourVector p1(1, 2, 3, 4);
        nuchic::FourVector p2(p1);
        CHECK(p1 == p2);

        nuchic::FourVector p3(std::move(p1));
        CHECK(p2 == p3);

        nuchic::FourVector p4 = p3;
        nuchic::FourVector p5 = std::move(p2);
        CHECK(p4 == p5);
    }
}

TEST_CASE("Accessors work as expected", "[Vectors]") {
    SECTION("Three Vector access") {
        nuchic::ThreeVector p(1, 2, 3);
        CHECK(p[0] == 1);
        CHECK(p[1] == 2);
        CHECK(p[2] == 3);
        CHECK_THROWS_AS(p[3], std::range_error);
    }

    SECTION("Four Vector access") {
        nuchic::FourVector p(1, 2, 3, 4);
        CHECK(p[0] == 1);
        CHECK(p[1] == 2);
        CHECK(p[2] == 3);
        CHECK(p[3] == 4);
        CHECK_THROWS_AS(p[4], std::range_error);
    }
}

TEST_CASE("Three Vector Overloaded Operators work as expected", "[Vectors]") {
    SECTION("Addition") {
        nuchic::ThreeVector p1(1, 2, 3), p2(1, 2, 3), p3(2, 4, 6);
        CHECK(p1+p2 == p3);

        p1 += p2;
        CHECK(p1 == p3);
    }

    SECTION("Subtraction and Negation") {
        nuchic::ThreeVector p1(1, 2, 3), p2(1, 2, 3), p3(0, 0, 0), p4(-1, -2, -3);
        CHECK(p1-p2 == p3);
        CHECK(-p1 == p4);

        p1 -= p2;
        CHECK(p1 == p3);
    }

    SECTION("Multiplication") {
        nuchic::ThreeVector p1(1, 2, 3), p2(3, 6, 9);
        constexpr double scalar1 = 3, scalar2 = 14;

        CHECK(p1*p1 == scalar2);
        CHECK(p1*p1 == p1.Dot(p1));
        CHECK(scalar1 * p1 == p1 * scalar1);
        CHECK(scalar1 * p1 == p2);

        p1 *= scalar1;
        CHECK(p1 == p2);
    }

    SECTION("Division") {
        nuchic::ThreeVector p1(1, 2, 3), p2(3, 6, 9);
        constexpr double scalar = 3.0;

        CHECK(p2 / scalar == p1);

        p2 /= scalar;
        CHECK(p2 == p1);
    }
}

TEST_CASE("Three Vector Functions work as expected", "[Vectors]") {
    SECTION("Magnitude") {
        nuchic::ThreeVector p1(4, 3, 2);
        constexpr double magnitude = 29;

        CHECK(p1.Magnitude2() == magnitude);
        CHECK(p1.P2() == magnitude);
        CHECK(p1.Magnitude() == sqrt(magnitude));
        CHECK(p1.P() == sqrt(magnitude)); 
    }

    SECTION("Transverse momentum and Angles") {
        nuchic::ThreeVector p1(1, 2, 3), p2(0, 0, 4), p3(1, 1, 1);
        constexpr double pt2 = 5;

        CHECK(p1.Pt2() == pt2);
        CHECK(p1.Pt() == sqrt(pt2));
        CHECK(p2.Theta() == 0);
        CHECK(p2.Phi() == 0);
        CHECK(p3.Phi() == M_PI_4);
    }

    SECTION("Cross Product") {
        nuchic::ThreeVector p1(3, 2, 1), p2(1, 2, 3), p3(4, -8, 4);

        CHECK(p1.Cross(p2) == p3);
        CHECK(p1.Cross(p2) == -p2.Cross(p1));
    }

    SECTION("Unit Vector") {
        nuchic::ThreeVector p1(3, 3, 3), p2(1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3));
        auto p1Unit = p1.Unit();

        CHECK((p1Unit[0] == Approx(p2[0]) && p1Unit[1] == Approx(p2[1])
               && p1Unit[2] == Approx(p2[2])));
        CHECK(p1/p1.Magnitude() == p1.Unit());
    }
}

TEST_CASE("Four Vector Overloaded Operators work as expected", "[Vectors]") {
    SECTION("Addition") {
        nuchic::FourVector p1(1, 2, 3, 4), p2(4, 3, 2, 1), p3(5, 5, 5, 5);
        CHECK(p1+p2 == p3);

        p1 += p2;
        CHECK(p1 == p3);
    }

    SECTION("Subtraction and Negation") {
        nuchic::FourVector p1(1, 2, 3, 4), p2(1, 2, 3, 4), p3(0, 0, 0, 0), p4(-1, -2, -3, -4);
        CHECK(p1-p2 == p3);
        CHECK(-p1 == p4);

        p1 -= p2;
        CHECK(p1 == p3);
    }

    SECTION("Multiplication") {
        nuchic::FourVector p1(1, 2, 3, 4), p2(3, 6, 9, 12);
        constexpr double scalar1 = 3, scalar2 = 16-9-4-1;

        CHECK(p1*p1 == scalar2);
        CHECK(p1*p1 == p1.Dot(p1));
        CHECK(scalar1 * p1 == p1 * scalar1);
        CHECK(scalar1 * p1 == p2);

        p1 *= scalar1;
        CHECK(p1 == p2);
    }

    SECTION("Division") {
        nuchic::FourVector p1(1, 2, 3, 4), p2(3, 6, 9, 12);
        constexpr double scalar = 3.0;

        CHECK(p2 / scalar == p1);

        p2 /= scalar;
        CHECK(p2 == p1);
    }
}

TEST_CASE("Four Vector Functions work as expected", "[Vectors]") {
    SECTION("Magnitude and Mass") {
        nuchic::FourVector p(1, 2, 3, 4);
        constexpr double mass = 2;

        CHECK(p.Magnitude2() == mass);
        CHECK(p.M2() == mass);
        CHECK(p.Magnitude() == sqrt(mass));
        CHECK(p.M() == sqrt(mass)); 
    }

    SECTION("Momentum and Transverse Momentum") {
        nuchic::FourVector p(1, 2, 3, 4);
        constexpr double pvec2 = 14;
        constexpr double pt2 = 5;

        CHECK(p.P2() == pvec2);
        CHECK(p.P() == sqrt(pvec2));
        CHECK(p.Pt2() == pt2);
        CHECK(p.Pt() == sqrt(pt2));
    }

    SECTION("Angles") {
        nuchic::FourVector p1(0, 0, 4, 4), p2(1, 1, 1, 4), p3(1, -1, 1, 4);

        CHECK(p1.Theta() == 0);
        CHECK(p1.Phi() == 0);
        CHECK(p2.Phi() == M_PI_4);
        CHECK(p3.Phi() == 7*M_PI_4);
    }

    SECTION("Cross Product") {
        nuchic::FourVector p1(3, 2, 1, 4), p2(1, 2, 3, 4), p3(4, -8, 4, 0);

        CHECK(p1.Cross(p2) == p3);
        CHECK(p1.Cross(p2) == -p2.Cross(p1));
    }

    SECTION("Boost") {
        nuchic::FourVector p1(3, 2, 1, 4), p2(12, 2, -3, 25);
        auto beta = p2.BoostVector();
        auto p3 = p1.Boost(beta).Boost(-beta);

        CHECK((p3[0] == Approx(p1[0]) && p3[1] == Approx(p1[1])
               && p3[2] == Approx(p1[2]) && p3[3] == Approx(p1[3])));
    }

    SECTION("Rapidity") {
        nuchic::FourVector p1(1, 2, 3, 4);
        constexpr double rapidity = 0.9729550745276566;

        CHECK(p1.Rapidity() == Approx(rapidity));
    }
    
    SECTION("DeltaR") {
        nuchic::FourVector p1(1, 2, 3, 4), p2(3, 2, 1, 4);
        double DEta = p1.Rapidity() - p2.Rapidity();
        double DPhi = p1.Phi() - p2.Phi();
        double DR = sqrt(DEta*DEta + DPhi*DPhi);

        CHECK(p1.DeltaR(p2) == p2.DeltaR(p1));
        CHECK(p1.DeltaR(p2) == DR);
    }
}
