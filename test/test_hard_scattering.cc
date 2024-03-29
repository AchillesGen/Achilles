#include "catch2/catch.hpp"

#include "mock_classes.hh"
#include "Achilles/HardScattering.hh"

#ifdef ENABLE_BSM

TEST_CASE("CrossSection", "[HardScattering]") {
    achilles::Process_Info info; 
    info.m_ids = {achilles::PID::electron(), achilles::PID::electron()};
    info.m_states = {{{achilles::PID::proton()}, {achilles::PID::proton()}}};
    info.m_mom_map = {{0, achilles::PID::electron().AsInt()},
                      {1, achilles::PID::proton().AsInt()},
                      {2, achilles::PID::electron().AsInt()},
                      {3, achilles::PID::proton().AsInt()}};
    std::vector<achilles::FourVector> momentum = {{100, 0, 0, 100}, {100, 0, 0, -100},
                                                {100, 50, 0, 50}, {100, -50, 0, -50}};
    std::vector<achilles::NuclearModel::Currents> hCurrent(1);
    std::vector<achilles::NuclearModel::FFInfoMap> ffInfo(3);
    std::vector<int> pids = {11, 2212, 11, 2212};
    std::vector<std::array<double, 4>> mom = {{0.1, 0, 0, 0.1}, {0.1, 0, 0, -0.1},
                                              {0.1, 0.05, 0, 0.05}, {0.1, -0.05, 0, -0.05}};
    achilles::SherpaInterface::LeptonCurrents lCurrent;

    lCurrent[23] = {{10, 10, 10, 10}};
    hCurrent[0][23] = {{10, 10, 10, 10}};
    ffInfo[0][23] = {{achilles::FormFactorInfo::Type::F1p, 1}};
    ffInfo[1][23] = {{achilles::FormFactorInfo::Type::F1n, 1}};
    ffInfo[2][23] = {{achilles::FormFactorInfo::Type::FCoh, 1}};

    MockEvent event;
    REQUIRE_CALL(event, Momentum())
        .TIMES(3)
        .LR_RETURN((momentum));

    auto model = std::make_unique<MockNuclearModel>();
    size_t nspins = 1;
    REQUIRE_CALL(*model, CalcCurrents(trompeloeil::_, ffInfo))
        .TIMES(1)
        .LR_RETURN((hCurrent));
    REQUIRE_CALL(*model, NSpins())
        .TIMES(2)
        .LR_RETURN((nspins));

    trompeloeil::sequence seq;

    std::vector<METOOLS::Spin_Amplitudes> spin_amps;
    spin_amps.emplace_back(std::vector<int>{2, 2, 2, 2}, 0);
    spin_amps[0][0] = {-0.2, 0};
    auto sherpa = new MockSherpaInterface;
    REQUIRE_CALL(*sherpa, Calc(pids, mom, 100))
        .TIMES(1)
        .LR_RETURN((lCurrent));
    REQUIRE_CALL(*sherpa, FormFactors(achilles::PID::proton(), 23))
        .TIMES(1)
        .LR_RETURN((ffInfo[0].at(23)));
    REQUIRE_CALL(*sherpa, FormFactors(achilles::PID::neutron(), 23))
        .TIMES(1)
        .LR_RETURN((ffInfo[1].at(23)));
    REQUIRE_CALL(*sherpa, FormFactors(achilles::PID::carbon(), 23))
        .TIMES(1)
        .LR_RETURN((ffInfo[2].at(23)));
    REQUIRE_CALL(*sherpa, FillAmplitudes(trompeloeil::_))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_SIDE_EFFECT(_1.emplace_back(std::vector<int>{2, 2, 2, 2}, 0););
    REQUIRE_CALL(*sherpa, FillAmplitudes(spin_amps))
        .TIMES(1)
        .IN_SEQUENCE(seq);

    achilles::HardScattering scattering;
    scattering.SetSherpa(sherpa);
    scattering.SetProcess(info);
    scattering.SetNuclear(std::move(model));
    auto xsec = scattering.CrossSection(event);
    CHECK(xsec[0] == Approx(20633.001927574882));
}

TEST_CASE("FillEvent", "[HardScattering]") {
    achilles::Process_Info info; 
    info.m_ids = {achilles::PID::electron(), achilles::PID::electron()};
    info.m_states = {{{achilles::PID::proton()}, {achilles::PID::proton()}}};
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
