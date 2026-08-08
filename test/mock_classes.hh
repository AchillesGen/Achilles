#ifndef MOCK_CLASSES
#define MOCK_CLASSES

#include "catch2/catch.hpp"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#include "catch2/trompeloeil.hpp"
#pragma GCC diagnostic pop

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
#include "Plugins/Sherpa/SherpaInterface.hh"
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

// Nucleus is a concrete type, so tests build a real one. Its density profile is
// tabulated rather than read from a file, and GenerateConfig is steered through
// the (real) Density interface, which MockDensity implements.
inline std::shared_ptr<achilles::Nucleus> MakeNucleus(size_t Z = 6, size_t A = 12,
                                                      double radius = 10, double rho = 0.16) {
    auto nucleus = std::make_shared<achilles::Nucleus>();
    nucleus->Initialize(Z, A);
    nucleus->SetRadius(radius);
    nucleus->SetFermiMomentum(225);
    std::vector<double> radii, density;
    static constexpr size_t nknots = 20;
    for(size_t i = 0; i < nknots; ++i) {
        radii.push_back(radius * static_cast<double>(i) / (nknots - 1));
        density.push_back(rho);
    }
    nucleus->SetDensityProfile(radii, density, density);
    return nucleus;
}

// As above, but with the density profile sampled from `profile` so a test can
// pin Rho(r) to a specific function.
inline std::shared_ptr<achilles::Nucleus>
MakeNucleus(size_t Z, size_t A, double radius, const std::function<double(double)> &profile,
            size_t nknots = 400) {
    auto nucleus = std::make_shared<achilles::Nucleus>();
    nucleus->Initialize(Z, A);
    nucleus->SetRadius(radius);
    nucleus->SetFermiMomentum(225);
    std::vector<double> radii, density;
    for(size_t i = 0; i < nknots; ++i) {
        radii.push_back(radius * static_cast<double>(i) / static_cast<double>(nknots - 1));
        density.push_back(profile(radii.back()));
    }
    nucleus->SetDensityProfile(radii, density, density);
    return nucleus;
}

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
    IMPLEMENT_MOCK1(FillAmplitudes);
};
*/

class MockInteraction : public trompeloeil::mock_interface<achilles::Interaction> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_CONST_MOCK0(InitialStates);
    IMPLEMENT_CONST_MOCK3(CrossSection);
    IMPLEMENT_CONST_MOCK4(GenerateMomentum);
    IMPLEMENT_CONST_MOCK0(GetName);
};

// Beam itself is a concrete container, so mock the flux it holds and put that in
// a real Beam. The plain properties are the ones the Beam constructor queries
// and no test asserts on.
class MockFluxType : public trompeloeil::mock_interface<achilles::FluxType> {
  public:
    static constexpr bool trompeloeil_movable_mock = true;
    MockFluxType(int nvars = 1, double emin = 0, double emax = 1000)
        : m_nvars{nvars}, m_emin{emin}, m_emax{emax} {}
    int NVariables() const override { return m_nvars; }
    double MinEnergy() const override { return m_emin; }
    double MaxEnergy() const override { return m_emax; }
    std::string Type() const override { return "Mock"; }

    IMPLEMENT_CONST_MOCK2(Flux);
    IMPLEMENT_CONST_MOCK3(GenerateWeight);
    IMPLEMENT_CONST_MOCK1(EvaluateFlux);

  private:
    int m_nvars;
    double m_emin, m_emax;
};

inline std::shared_ptr<achilles::Beam> MakeBeam(achilles::PID pid,
                                                std::shared_ptr<achilles::FluxType> flux) {
    return std::make_shared<achilles::Beam>(achilles::Beam::BeamMap{{pid, std::move(flux)}});
}

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

class MockUnweighter : public trompeloeil::mock_interface<achilles::Unweighter> {
    static constexpr bool trompeloeil_movable_mock = true;
    IMPLEMENT_MOCK1(AddEvent);
    IMPLEMENT_MOCK1(AcceptEvent);
    IMPLEMENT_MOCK0(MaxValue);
    IMPLEMENT_CONST_MOCK1(SaveState);
    IMPLEMENT_MOCK1(LoadState);
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
