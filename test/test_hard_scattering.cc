#include "catch2/catch.hpp"

#include "mock_classes.hh"
#include "Achilles/HardScattering.hh"

#ifdef ENABLE_BSM

TEST_CASE("CrossSection", "[HardScattering]") {
    achilles::Process_Info info; 
    info.ids = {achilles::PID::electron(), achilles::PID::electron()};
    info.state = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    info.m_mom_map = {{0, achilles::PID::electron().AsInt()},
                      {1, achilles::PID::proton().AsInt()},
                      {2, achilles::PID::electron().AsInt()},
                      {3, achilles::PID::proton().AsInt()}};
    achilles::Process_Group group;
    group.AddProcess(info);
    std::vector<achilles::FourVector> momentum = {{100, 0, 0, 100}, {100, 0, 0, -100},
                                                {100, 50, 0, 50}, {100, -50, 0, -50}};
    std::vector<achilles::NuclearModel::Currents> hCurrent(1);
    std::vector<achilles::NuclearModel::FFInfoMap> ffInfo(1);
    std::vector<int> pids = {11, 2212, 11, 2212};
    std::vector<std::array<double, 4>> mom = {{0.1, 0, 0, 0.1}, {0.1, 0, 0, -0.1},
                                              {0.1, 0.05, 0, 0.05}, {0.1, -0.05, 0, -0.05}};
    achilles::SherpaMEs::LeptonCurrents lCurrent;

    lCurrent[23] = {{10, 10, 10, 10}};
    hCurrent[0][23] = {{10, 10, 10, 10}};
    ffInfo[0][23] = {{achilles::FormFactorInfo::Type::F1p, 1}};

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

    auto sherpa = new MockSherpaME;
    REQUIRE_CALL(*sherpa, Calc(pids, mom, 100))
        .TIMES(1)
        .LR_RETURN((lCurrent));
    REQUIRE_CALL(*sherpa, FormFactors(achilles::PID::proton(), 23))
        .TIMES(1)
        .LR_RETURN((ffInfo[0].at(23)));

    achilles::HardScattering scattering;
    scattering.SetSherpa(sherpa);
    scattering.SetProcessGroup(group);
    scattering.SetNuclear(std::move(model));
    auto xsec = scattering.CrossSection(event);
    CHECK(xsec[0] == Approx(20633.001927574882));
}

TEST_CASE("FillEvent", "[HardScattering]") {
    achilles::Process_Info info; 
    info.ids = {achilles::PID::electron(), achilles::PID::electron()};
    info.state = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    achilles::Process_Group group;
    group.AddProcess(info);
    std::vector<double> xsecs{10};
    MockEvent event; 
    REQUIRE_CALL(event, InitializeLeptons(info)).TIMES(1);
    REQUIRE_CALL(event, InitializeHadrons(info)).TIMES(1);

    achilles::HardScattering scattering;
    scattering.SetProcessGroup(group);
    scattering.FillEvent(event, xsecs);
}

#endif
