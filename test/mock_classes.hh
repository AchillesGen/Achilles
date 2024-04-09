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
    IMPLEMENT_CONST_MOCK0(GetReference);
    achilles::PotentialVals operator()(double p, double r) const override { return call_op(p, r); }
    MAKE_CONST_MOCK2(call_op, achilles::PotentialVals(const double &, const double &));
    IMPLEMENT_CONST_MOCK2(Hamiltonian);
};

class MockNucleus : public trompeloeil::mock_interface<achilles::Nucleus> {
    static constexpr bool trompeloeil_movable_mock = true;
    MAKE_MOCK0(Nucleons, achilles::Particles &(), noexcept override);
    IMPLEMENT_MOCK0(GenerateConfig);
    IMPLEMENT_CONST_MOCK0(ID);
    MAKE_CONST_MOCK0(Radius, const double &(), noexcept override);
    MAKE_CONST_MOCK1(Rho, double(const double &), noexcept override);
    MAKE_CONST_MOCK0(NNucleons, size_t(), noexcept override);
    MAKE_CONST_MOCK0(GetPotential, std::shared_ptr<achilles::Potential>(), noexcept override);
    MAKE_CONST_MOCK0(ProtonsIDs, std::vector<size_t>(), noexcept override);
    MAKE_CONST_MOCK0(NeutronsIDs, std::vector<size_t>(), noexcept override);
};

class MockNuclearModel : public trompeloeil::mock_interface<achilles::NuclearModel> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK0(Mode);
    IMPLEMENT_CONST_MOCK1(PhaseSpace);
    IMPLEMENT_CONST_MOCK4(CalcCurrents);
    IMPLEMENT_CONST_MOCK1(AllowedStates);
    IMPLEMENT_CONST_MOCK0(NSpins);
    IMPLEMENT_CONST_MOCK3(InitialStateWeight);
    IMPLEMENT_CONST_MOCK0(GetName);
    IMPLEMENT_CONST_MOCK0(InspireHEP);
};

/*
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
    IMPLEMENT_CONST_MOCK3(Flux);
    IMPLEMENT_CONST_MOCK0(BeamIDs);
    IMPLEMENT_CONST_MOCK4(GenerateWeight);
    IMPLEMENT_CONST_MOCK0(NVariables);
    IMPLEMENT_CONST_MOCK2(EvaluateFlux);
};

class MockEvent : public trompeloeil::mock_interface<achilles::Event> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK0(CurrentNucleus);
    IMPLEMENT_CONST_MOCK0(CurrentNucleus);
    IMPLEMENT_MOCK0(Hadrons);
    IMPLEMENT_CONST_MOCK0(Hadrons);
    IMPLEMENT_MOCK0(Leptons);
    IMPLEMENT_CONST_MOCK0(Leptons);
    MAKE_CONST_MOCK0(Momentum, const std::vector<achilles::FourVector> &());
    MAKE_MOCK0(Momentum, std::vector<achilles::FourVector> &());
    IMPLEMENT_CONST_MOCK0(Particles);
    IMPLEMENT_CONST_MOCK0(Remnant);
    MAKE_CONST_MOCK0(Weight, const double &());
    MAKE_MOCK0(Weight, double &());
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
    IMPLEMENT_MOCK3(FinalState);
    IMPLEMENT_MOCK0(build);
};

class MockUnweighter : public trompeloeil::mock_interface<achilles::Unweighter> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK1(AddEvent);
    IMPLEMENT_MOCK1(AcceptEvent);
    IMPLEMENT_MOCK0(MaxValue);
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
    IMPLEMENT_MOCK1(AddProcess);
    IMPLEMENT_MOCK0(Validate);
    IMPLEMENT_MOCK4(SetupChannels);
    IMPLEMENT_MOCK0(GetNuclearModel);

    // Required factory methods
    static std::unique_ptr<XSecBackend> Construct() { return std::move(self); }
    static std::string Name() { return "Mock"; }
    static std::unique_ptr<MockBackend> self;
    static void SetSelf(std::unique_ptr<MockBackend> backend) { self = std::move(backend); }
};

#endif
