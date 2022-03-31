#ifndef MOCK_CLASSES
#define MOCK_CLASSES

#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"

// Includes to mock
#include "nuchic/Interactions.hh"
#include "nuchic/Potential.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Beams.hh"
#include "nuchic/Event.hh"
#include "nuchic/NuclearModel.hh"
#include "nuchic/FormFactor.hh"
#include "nuchic/PhaseSpaceBuilder.hh"
#include "plugins/Sherpa/SherpaMEs.hh"

class MockDensity : public trompeloeil::mock_interface<nuchic::Density> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(GetConfiguration);
};

class MockPotential : public trompeloeil::mock_interface<nuchic::Potential> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK2(Hamiltonian);
    IMPLEMENT_CONST_MOCK0(GetReference);
    nuchic::PotentialVals operator()(const double &p, const double &r) const override {
        return call_op(p, r);
    }
    MAKE_CONST_MOCK2(call_op, nuchic::PotentialVals(const double&, const double&));
};

class MockNucleus : public trompeloeil::mock_interface<nuchic::Nucleus> {
    static constexpr bool trompeloeil_movable_mock = true;
    MAKE_MOCK0(Nucleons, nuchic::Particles&(), noexcept override);
    IMPLEMENT_MOCK0(GenerateConfig);
    MAKE_CONST_MOCK0(Radius, const double&(), noexcept override);
    MAKE_CONST_MOCK1(Rho, double(const double&), noexcept override);
    MAKE_CONST_MOCK0(NNucleons, size_t(), noexcept override);
    MAKE_CONST_MOCK0(GetPotential, std::shared_ptr<nuchic::Potential>(), noexcept override);
    MAKE_CONST_MOCK1(
};

class MockNuclearModel : public trompeloeil::mock_interface<nuchic::NuclearModel> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK0(Mode);
    IMPLEMENT_CONST_MOCK0(PhaseSpace);
    IMPLEMENT_CONST_MOCK2(CalcCurrents);
    IMPLEMENT_CONST_MOCK1(AllowedStates);
    IMPLEMENT_CONST_MOCK0(NSpins);
    IMPLEMENT_CONST_MOCK2(FillNucleus);
};

class MockSherpaME : public trompeloeil::mock_interface<nuchic::SherpaMEs> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK3(Calc);
    IMPLEMENT_CONST_MOCK2(FormFactors);
};

class MockInteraction : public trompeloeil::mock_interface<nuchic::Interactions> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK2(CrossSection);
    IMPLEMENT_CONST_MOCK3(MakeMomentum);
    IMPLEMENT_CONST_MOCK3(FinalizeMomentum);
    IMPLEMENT_CONST_MOCK0(Name);
};

class MockBeam : public trompeloeil::mock_interface<nuchic::Beam> {
    IMPLEMENT_CONST_MOCK2(Flux); 
    IMPLEMENT_CONST_MOCK0(BeamIDs);
    IMPLEMENT_CONST_MOCK3(GenerateWeight);
    IMPLEMENT_CONST_MOCK0(NVariables);
};

class MockEvent : public trompeloeil::mock_interface<nuchic::Event> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(CurrentNucleus);
    IMPLEMENT_MOCK0(Hadrons);
    IMPLEMENT_MOCK1(InitializeLeptons);
    IMPLEMENT_MOCK1(InitializeHadrons);
    MAKE_CONST_MOCK0(Momentum, const std::vector<nuchic::FourVector>&());
    MAKE_MOCK0(Momentum, std::vector<nuchic::FourVector>&());
    IMPLEMENT_CONST_MOCK0(Particles);
    IMPLEMENT_CONST_MOCK0(Remnant);
    IMPLEMENT_CONST_MOCK0(Weight);
};

class MockFormFactor : public trompeloeil::mock_interface<nuchic::FormFactor> {
    static constexpr bool trompeloeil_movable_mock = true;
    nuchic::FormFactor::Values operator()(double Q2) const override { return call_op(Q2); }
    MAKE_CONST_MOCK1(call_op, nuchic::FormFactor::Values(double));
};

class MockFormFactorBuilder : public trompeloeil::mock_interface<nuchic::FormFactorBuilder> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK2(Vector);
    IMPLEMENT_MOCK2(AxialVector);
    IMPLEMENT_MOCK2(Coherent);
    IMPLEMENT_MOCK0(build);
};

class MockMapper : public trompeloeil::mock_interface<nuchic::Mapper<nuchic::FourVector>> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK2(GeneratePoint);
    IMPLEMENT_MOCK2(GenerateWeight);
    IMPLEMENT_CONST_MOCK0(NDims);
    IMPLEMENT_MOCK1(SetMasses);
    IMPLEMENT_CONST_MOCK0(Masses);
    IMPLEMENT_CONST_MOCK0(ToYAML);
};

class MockPSBuilder : public trompeloeil::mock_interface<nuchic::PSBuilder> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK2(Beam);
    IMPLEMENT_MOCK3(Hadron);
    IMPLEMENT_MOCK2(FinalState);
};

#endif
