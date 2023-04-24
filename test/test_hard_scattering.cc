#include "catch2/catch.hpp"

#include "Achilles/HardScattering.hh"
#include "mock_classes.hh"

#ifdef ACHILLES_SHERPA_INTERFACE

void build_spinamps(std::vector<METOOLS::Spin_Amplitudes> &amps) {
    std::vector<ATOOLS::Flavour> flavs;
    flavs.emplace_back(11);
    flavs.emplace_back(11);
    amps.push_back(METOOLS::Spin_Amplitudes(flavs, std::complex<double>(0.0, 0.0)));
    for(size_t i = 0; i < 4; ++i) amps.back().Insert(0, i);
}

TEST_CASE("CrossSection", "[HardScattering]") {
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    info.m_mom_map = {{0, achilles::PID::electron().AsInt()},
                      {1, achilles::PID::proton().AsInt()},
                      {2, achilles::PID::electron().AsInt()},
                      {3, achilles::PID::proton().AsInt()}};
    std::vector<achilles::FourVector> momentum = {
        {100, 0, 0, 100}, {100, 0, 0, -100}, {100, 50, 0, 50}, {100, -50, 0, -50}};
    std::vector<achilles::NuclearModel::Currents> hCurrent(1);
    std::vector<achilles::NuclearModel::FFInfoMap> ffInfo(3);
    std::vector<int> pids = {11, 2212, 11, 2212};
    std::vector<std::array<double, 4>> mom = {
        {0.1, 0, 0, 0.1}, {0.1, 0, 0, -0.1}, {0.1, 0.05, 0, 0.05}, {0.1, -0.05, 0, -0.05}};
    achilles::SherpaInterface::LeptonCurrents lCurrent;

    // lCurrent[23] = {{10, 10, 10, 10}};
    // hCurrent[0][23] = {{10, 10, 10, 10}};
    // ffInfo[0][23] = {{achilles::FormFactorInfo::Type::F1p, 1}};
    // ffInfo[1][23] = {{achilles::FormFactorInfo::Type::F1n, 1}};
    // ffInfo[2][23] = {{achilles::FormFactorInfo::Type::FCoh, 1}};

    // MockEvent event;
    // REQUIRE_CALL(event, Momentum()).TIMES(3).LR_RETURN((momentum));

    // auto model = std::make_unique<MockNuclearModel>();
    // size_t nspins = 1;
    // REQUIRE_CALL(*model, CalcCurrents(trompeloeil::_, ffInfo)).TIMES(1).LR_RETURN((hCurrent));
    // REQUIRE_CALL(*model, NSpins()).TIMES(2).LR_RETURN((nspins));

    ATOOLS::s_kftable[kf_none] =
        new ATOOLS::Particle_Info(0, -1, 0, 0, 0, 0, -1, 0, 1, 0, "no_particle", "no_particle",
                                  "no_particle", "no_particle", 1, 1);
    ATOOLS::s_kftable[kf_e] = new ATOOLS::Particle_Info(11, 0.000511, .0, -3, 0, 1, 0, 1, 1, 0,
                                                        "e-", "e+", "e^{-}", "e^{+}");

    auto sherpa = new MockSherpaInterface;
    REQUIRE_CALL(*sherpa, Calc(pids, mom, 100)).TIMES(1).LR_RETURN((lCurrent));
    REQUIRE_CALL(*sherpa, FormFactors(achilles::PID::proton(), 23))
        .TIMES(1)
        .LR_RETURN((ffInfo[0].at(23)));
    REQUIRE_CALL(*sherpa, FormFactors(achilles::PID::neutron(), 23))
        .TIMES(1)
        .LR_RETURN((ffInfo[1].at(23)));
    REQUIRE_CALL(*sherpa, FormFactors(achilles::PID::carbon(), 23))
        .TIMES(1)
        .LR_RETURN((ffInfo[2].at(23)));
    REQUIRE_CALL(*sherpa, FillAmplitudes(trompeloeil::_)).TIMES(2).SIDE_EFFECT(build_spinamps(_1));

    achilles::HardScattering scattering;
    scattering.SetSherpa(sherpa);
    scattering.SetProcess(info);
    // scattering.SetNuclear(std::move(model));
    // auto xsec = scattering.CrossSection(event);
    // CHECK(xsec[0] == Approx(20633.001927574882));

    delete sherpa;
}

TEST_CASE("FillEvent", "[HardScattering]") {
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    std::vector<double> xsecs;
    MockEvent event;
    REQUIRE_CALL(event, InitializeLeptons(info)).TIMES(1);
    REQUIRE_CALL(event, InitializeHadrons(info)).TIMES(1);
    auto model = std::make_unique<MockNuclearModel>();
    REQUIRE_CALL(*model, FillNucleus(trompeloeil::_, xsecs)).TIMES(1).RETURN(true);

    achilles::HardScattering scattering;
    scattering.SetProcess(info);
    scattering.SetNuclear(std::move(model));
    scattering.FillEvent(event, xsecs);
}

#endif
