#include "catch2/catch.hpp"

#include "Achilles/Poincare.hh"
#include "Approx.hh"

TEST_CASE("Boost", "[Poincare]") {
    achilles::FourVector p1{1000, 0, 0, 500};
    achilles::FourVector p2{100, 0, 0, 100};
    achilles::FourVector expected = p2;

    achilles::Poincare poincare(p1, -1);
    poincare.Boost(p2);
    poincare.BoostBack(p2);
    CHECK_THAT(p2, FourVectorApprox(expected).margin(1e-16));
}

TEST_CASE("Rotation", "[Poincare]") {
    achilles::FourVector p1{1000, 0, 250, 500};
    achilles::FourVector p2{100, 0, 0, 100};
    achilles::FourVector expected = p2;

    achilles::Poincare poincare(p1, -1);
    poincare.Rotate(p2);
    poincare.RotateBack(p2);
    CHECK_THAT(p2, FourVectorApprox(expected).margin(1e-16));
}

TEST_CASE("Lambda", "[Poincare]") {
    achilles::FourVector p1{1000, 0, 250, 500};
    achilles::FourVector p2{100, 0, 0, 100};
    achilles::FourVector p3{1000, 250, 0, 500};
    achilles::FourVector expected = p2;

    achilles::Poincare poincare(p1, p3, 3);
    poincare.Lambda(p2);
    poincare.LambdaBack(p2);
    CHECK_THAT(p2, FourVectorApprox(expected).margin(1e-16));
}

TEST_CASE("Invert", "[Poincare]") {
    achilles::FourVector p1{1000, 0, 250, 500};
    achilles::FourVector p2{100, 0, 0, 100};
    achilles::FourVector expected = p2;

    achilles::Poincare poincare(p1, -1);
    poincare.Boost(p2);
    poincare.Invert();
    poincare.Boost(p2);
    CHECK_THAT(p2, FourVectorApprox(expected).margin(1e-8));
}
