#ifndef MOCK_CLASSES
#define MOCK_CLASSES

#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"

// Includes to mock
#include "Achilles/Beams.hh"
#include "Achilles/Event.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/Potential.hh"
#include "plugins/Sherpa/SherpaInterface.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#undef THROW
#include "METOOLS/Main/Spin_Structure.H"
#undef THROW
#define THROW TROMPELOEIL_THROW
#pragma GCC diagnostic pop
#endif

class MockDensity : public trompeloeil::mock_interface<achilles::Density> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(GetConfiguration);
};

class MockPotential : public trompeloeil::mock_interface<achilles::Potential> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK2(Hamiltonian);
    IMPLEMENT_CONST_MOCK0(GetReference);
    achilles::PotentialVals operator()(const double &p, const double &r) const override {
        return call_op(p, r);
    }
    MAKE_CONST_MOCK2(call_op, achilles::PotentialVals(const double &, const double &));
};

class MockNucleus : public trompeloeil::mock_interface<achilles::Nucleus> {
    static constexpr bool trompeloeil_movable_mock = true;
    MAKE_MOCK0(Nucleons, achilles::Particles &(), noexcept override);
    IMPLEMENT_MOCK0(GenerateConfig);
    MAKE_CONST_MOCK0(Radius, const double &(), noexcept override);
    MAKE_CONST_MOCK1(Rho, double(const double &), noexcept override);
    MAKE_CONST_MOCK0(NNucleons, size_t(), noexcept override);
    MAKE_CONST_MOCK0(GetPotential, std::shared_ptr<achilles::Potential>(), noexcept override);
};

class MockNuclearModel : public trompeloeil::mock_interface<achilles::NuclearModel> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK0(Mode);
    IMPLEMENT_CONST_MOCK0(PhaseSpace);
    IMPLEMENT_CONST_MOCK4(CalcCurrents);
    IMPLEMENT_CONST_MOCK1(AllowedStates);
    IMPLEMENT_CONST_MOCK0(NSpins);
    IMPLEMENT_CONST_MOCK2(FillNucleus);
};

class MockSherpaInterface : public trompeloeil::mock_interface<achilles::SherpaInterface> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK3(Calc);
    IMPLEMENT_CONST_MOCK2(FormFactors);
    IMPLEMENT_MOCK1(FillAmplitudes);
};

class MockInteraction : public trompeloeil::mock_interface<achilles::Interactions> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK2(CrossSection);
    IMPLEMENT_CONST_MOCK3(MakeMomentum);
    IMPLEMENT_CONST_MOCK3(FinalizeMomentum);
    IMPLEMENT_CONST_MOCK0(Name);
};

class MockBeam : public trompeloeil::mock_interface<achilles::Beam> {
    IMPLEMENT_CONST_MOCK3(Flux);
    IMPLEMENT_CONST_MOCK0(BeamIDs);
    IMPLEMENT_CONST_MOCK4(GenerateWeight);
    IMPLEMENT_CONST_MOCK0(NVariables);
};

class MockEvent : public trompeloeil::mock_interface<achilles::Event> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(CurrentNucleus);
    IMPLEMENT_MOCK0(Hadrons);
    IMPLEMENT_MOCK1(InitializeLeptons);
    IMPLEMENT_MOCK1(InitializeHadrons);
    MAKE_CONST_MOCK0(Momentum, const std::vector<achilles::FourVector> &());
    MAKE_MOCK0(Momentum, std::vector<achilles::FourVector> &());
    IMPLEMENT_CONST_MOCK0(Particles);
    IMPLEMENT_CONST_MOCK0(Remnant);
    MAKE_CONST_MOCK0(Weight, const double &());
    MAKE_MOCK0(Weight, double &());
};

class MockFormFactor : public trompeloeil::mock_interface<achilles::FormFactor> {
    static constexpr bool trompeloeil_movable_mock = true;
    achilles::FormFactor::Values operator()(double Q2) const override { return call_op(Q2); }
    MAKE_CONST_MOCK1(call_op, achilles::FormFactor::Values(double));
};

class MockFormFactorBuilder : public trompeloeil::mock_interface<achilles::FormFactorBuilder> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK2(Vector);
    IMPLEMENT_MOCK2(AxialVector);
    IMPLEMENT_MOCK2(Coherent);
    IMPLEMENT_MOCK0(build);
};

class MockMapper : public trompeloeil::mock_interface<achilles::Mapper<achilles::FourVector>> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK2(GeneratePoint);
    IMPLEMENT_MOCK2(GenerateWeight);
    IMPLEMENT_CONST_MOCK0(NDims);
    IMPLEMENT_MOCK1(SetMasses);
    IMPLEMENT_CONST_MOCK0(Masses);
    IMPLEMENT_CONST_MOCK0(ToYAML);
};

class MockPSBuilder : public trompeloeil::mock_interface<achilles::PSBuilder> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK3(Beam);
    IMPLEMENT_MOCK3(Hadron);
    IMPLEMENT_MOCK2(FinalState);
};

#endif
