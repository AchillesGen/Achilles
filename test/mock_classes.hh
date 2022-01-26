#ifndef MOCK_CLASSES
#define MOCK_CLASSES

#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"
#include "nuchic/Interactions.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Beams.hh"
#include "nuchic/Event.hh"

class MockDensity : public trompeloeil::mock_interface<nuchic::Density> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(GetConfiguration);
};

class MockNucleus : public trompeloeil::mock_interface<nuchic::Nucleus> {
    static constexpr bool trompeloeil_movable_mock = true;
    MAKE_MOCK0(Nucleons, nuchic::Particles&(), noexcept override);
    IMPLEMENT_MOCK0(GenerateConfig);
    MAKE_CONST_MOCK0(Radius, const double&(), noexcept override);
    MAKE_CONST_MOCK1(Rho, double(const double&), noexcept override);
    MAKE_CONST_MOCK0(NNucleons, size_t(), noexcept override);
};

class MockInteraction : public trompeloeil::mock_interface<nuchic::Interactions> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK2(CrossSection);
    IMPLEMENT_CONST_MOCK3(MakeMomentum);
    IMPLEMENT_CONST_MOCK0(Name);
};

class MockBeam : public trompeloeil::mock_interface<nuchic::Beam> {
    IMPLEMENT_CONST_MOCK2(Flux); 
};

class MockEvent : public trompeloeil::mock_interface<nuchic::Event> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(CurrentNucleus);
    IMPLEMENT_MOCK0(Hadrons);
};

#endif
