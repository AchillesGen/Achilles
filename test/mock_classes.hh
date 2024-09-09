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
#include "Achilles/Process.hh"
#include "Achilles/Unweighter.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/SherpaInterface.hh"
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
#else

/*
#include <complex>
#include <map>
#include <vector>

namespace achilles {

struct FormFactorInfo;
namespace METOOLS {
using Spin_Amplitudes = std::vector<std::complex<double>>;
}
class SherpaInterface {
  public:
    using LeptonCurrents = std::map<int, std::vector<std::vector<std::complex<double>>>>;
    virtual ~SherpaInterface() = default;
    virtual LeptonCurrents CalcCurrent(const std::vector<int> &,
                                       const std::vector<std::array<double, 4>> &,
                                       const double &) { return {}; }
    virtual LeptonCurrents CalcDifferential(const std::vector<int> &,
                                            const std::vector<std::array<double, 4>> &,
                                            const double &) { return {}; }
    virtual std::vector<FormFactorInfo> FormFactors(int, int) const;
    virtual void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes> &amps);
};

} // namespace achilles
*/
#endif

using achilles::RegistrableBackend;
using achilles::XSecBackend;

class MockDensity : public trompeloeil::mock_interface<achilles::Density> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(GetConfiguration);
};

class MockPotential : public trompeloeil::mock_interface<achilles::Potential> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK2(Hamiltonian);
    achilles::PotentialVals operator()(double p, double r) const override { return call_op(p, r); }
    MAKE_CONST_MOCK2(call_op, achilles::PotentialVals(const double &, const double &));
};

class MockNucleus : public trompeloeil::mock_interface<achilles::Nucleus> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(GenerateConfig);
    IMPLEMENT_CONST_MOCK0(ID);
    MAKE_CONST_MOCK0(Radius, const double &(), noexcept override);
    MAKE_CONST_MOCK1(Rho, double(const double &), noexcept override);
    MAKE_CONST_MOCK0(NNucleons, size_t(), noexcept override);
    MAKE_CONST_MOCK0(NProtons, size_t(), noexcept override);
    MAKE_CONST_MOCK0(NNeutrons, size_t(), noexcept override);
    MAKE_CONST_MOCK0(GetPotential, std::shared_ptr<achilles::Potential>(), noexcept override);
};

class MockNuclearModel : public trompeloeil::mock_interface<achilles::NuclearModel> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK0(Mode);
    IMPLEMENT_CONST_MOCK1(PhaseSpace);
    IMPLEMENT_CONST_MOCK5(CalcCurrents);
    IMPLEMENT_CONST_MOCK1(AllowedStates);
    IMPLEMENT_CONST_MOCK0(NSpins);
    IMPLEMENT_CONST_MOCK4(InitialStateWeight);
    IMPLEMENT_CONST_MOCK0(GetName);
    IMPLEMENT_CONST_MOCK0(PSName);
    IMPLEMENT_CONST_MOCK0(InspireHEP);
    IMPLEMENT_CONST_MOCK0(Frame);
};

/*
 * TODO: Figure out why this is broken!!
class MockSherpaInterface : public trompeloeil::mock_interface<achilles::SherpaInterface> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK3(CalcCurrent);
    IMPLEMENT_MOCK3(CalcDifferential);
    IMPLEMENT_MOCK1(FillAmplitudes);
    IMPLEMENT_CONST_MOCK2(FormFactors);
};
*/

class MockInteraction : public trompeloeil::mock_interface<achilles::Interaction> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK0(InitialStates);
    IMPLEMENT_CONST_MOCK2(CrossSection);
    IMPLEMENT_CONST_MOCK4(GenerateMomentum);
    IMPLEMENT_CONST_MOCK0(GetName);
};

class MockBeam : public trompeloeil::mock_interface<achilles::Beam> {
    IMPLEMENT_CONST_MOCK0(NVariables);
    IMPLEMENT_CONST_MOCK3(Flux);
    IMPLEMENT_CONST_MOCK4(GenerateWeight);
    IMPLEMENT_CONST_MOCK0(BeamIDs);
    IMPLEMENT_CONST_MOCK2(EvaluateFlux);
};

class MockEvent : public trompeloeil::mock_interface<achilles::Event> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK0(Remnant);
    IMPLEMENT_MOCK0(CurrentNucleus);
    IMPLEMENT_CONST_MOCK0(CurrentNucleus);

    IMPLEMENT_CONST_MOCK0(Particles);
    IMPLEMENT_MOCK0(Hadrons);
    IMPLEMENT_CONST_MOCK0(Hadrons);
    IMPLEMENT_MOCK0(Leptons);
    IMPLEMENT_CONST_MOCK0(Leptons);

    IMPLEMENT_MOCK0(Momentum);
    IMPLEMENT_CONST_MOCK0(Momentum);
    IMPLEMENT_MOCK0(Weight);
    IMPLEMENT_CONST_MOCK0(Weight);
    IMPLEMENT_CONST_MOCK0(History);
};

class MockFormFactor : public trompeloeil::mock_interface<achilles::FormFactor> {
    static constexpr bool trompeloeil_movable_mock = true;
    achilles::FormFactor::Values operator()(double Q2) const override { return call_op(Q2); }
    MAKE_CONST_MOCK1(call_op, achilles::FormFactor::Values(double));
};

class MockFormFactorBuilder : public trompeloeil::mock_interface<achilles::FormFactorBuilder> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(Reset);
    IMPLEMENT_MOCK2(Vector);
    IMPLEMENT_MOCK2(AxialVector);
    IMPLEMENT_MOCK2(Coherent);
    IMPLEMENT_MOCK2(ResonanceVector);
    IMPLEMENT_MOCK2(ResonanceAxial);
    IMPLEMENT_MOCK2(MesonExchangeVector);
    IMPLEMENT_MOCK2(MesonExchangeAxial);
    IMPLEMENT_MOCK2(Hyperon);
    IMPLEMENT_MOCK0(build);
};

class MockMapper : public trompeloeil::mock_interface<achilles::Mapper<achilles::FourVector>> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK2(GeneratePoint);
    IMPLEMENT_MOCK2(GenerateWeight);
    IMPLEMENT_CONST_MOCK0(NDims);
    IMPLEMENT_MOCK1(SetMasses);
    IMPLEMENT_CONST_MOCK0(Masses);
    IMPLEMENT_MOCK1(SetGaugeBosonMass);
};

class MockPSBuilder : public trompeloeil::mock_interface<achilles::PSBuilder> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK2(Beam);
    IMPLEMENT_MOCK2(Hadron);
    IMPLEMENT_MOCK2(FinalState);
    IMPLEMENT_MOCK0(build);
};

class MockUnweighter : public trompeloeil::mock_interface<achilles::Unweighter> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK1(AddEvent);
    IMPLEMENT_MOCK1(AcceptEvent);
    IMPLEMENT_MOCK0(MaxValue);
    IMPLEMENT_CONST_MOCK1(SaveState);
    IMPLEMENT_MOCK1(LoadState);
};

class MockProcess : public trompeloeil::mock_interface<achilles::Process> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK0(Info);
    IMPLEMENT_MOCK1(AddWeight);
    IMPLEMENT_MOCK1(Unweight);
    IMPLEMENT_CONST_MOCK6(ExtractMomentum);
    IMPLEMENT_CONST_MOCK6(ExtractParticles);
};

class MockBackend : public trompeloeil::mock_interface<XSecBackend>,
                    RegistrableBackend<MockBackend> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK2(CrossSection);
    IMPLEMENT_MOCK1(SetOptions);
    IMPLEMENT_MOCK1(SetSherpa);
    IMPLEMENT_MOCK1(AddNuclearModel);
    IMPLEMENT_MOCK0(GetNuclearModel);
    IMPLEMENT_MOCK1(AddProcess);
    IMPLEMENT_MOCK0(Validate);
    IMPLEMENT_MOCK4(SetupChannels);

    // Required factory methods
    static std::unique_ptr<XSecBackend> Construct() { return std::move(self); }
    static std::string Name() { return "Mock"; }
    static std::unique_ptr<MockBackend> self;
    static void SetSelf(std::unique_ptr<MockBackend> backend) { self = std::move(backend); }
};

#endif
