#include "Achilles/Process.hh"
#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include "Achilles/Particle.hh"
#include "Achilles/Version.hh"
#include "plugins/NuHepMC/NuHepMCWriter.hh"

TEST_CASE("Passes Validator", "[NuHepMC]") {
    // Setup dummy event
    achilles::Particle target(achilles::PID::carbon(),
                              {achilles::ParticleInfo(achilles::PID::carbon()).Mass(), 0, 0, 0});
    achilles::Particle nuc_in(achilles::PID::neutron(), {918.00661011686168, 23.585833577703873,
                                                         85.332385429710143, 52.378992899809319});
    achilles::Particle beam_part(achilles::PID::nu_muon(), {1e4, 0, 0, 1e4});
    achilles::Particle neutrino(achilles::PID::nu_muon(),
                                {5.3983343748755351e3, 0, 0, 5.3983343748755351e3});
    achilles::Particle nuc_out(
        achilles::PID::proton(),
        {9.6318094613481071e2, -1.8702102417549486e2, -4.7096228265225918e1, 1.0333739302250189e2});
    achilles::Particle lepton(achilles::PID::muon(), {5.3531600388575862e3, 2.1060685775319874e2,
                                                      1.3242861369493605e2, 5.3473759747528429e3});
    // TODO: Figure out how to handle remnant
    // achilles::Particle remnant(achilles::PID(1000060110));
    // remnant.Momentum() = target.Momentum() + neutrino.Momentum() - nuc_out.Momentum() -
    // lepton.Momentum();

    achilles::EventHistory history;
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {target}, {nuc_in}, achilles::EventHistoryNode::StatusCode::target);
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {beam_part}, {neutrino}, achilles::EventHistoryNode::StatusCode::beam);
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {nuc_in, neutrino}, {nuc_out, lepton},
                      achilles::EventHistoryNode::StatusCode::primary);

    const MockEvent event{};
    double wgt = 1.0;
    REQUIRE_CALL(event, History()).TIMES(1).LR_RETURN((history));
    REQUIRE_CALL(event, Weight()).TIMES(AT_LEAST(1)).LR_RETURN((wgt));

    // Force saving of HepMC file
    {
        auto beam = std::make_shared<MockBeam>();
        auto nucleus = std::make_shared<MockNucleus>();
        YAML::Node config = YAML::Load(R"config(
        Processes:
            - Leptons: [11, [11]]
        Unweighting:
            Name: Percentile
            percentile: 99
        Backend:
            Name: Mock
            Options:
        )config");

        auto model = std::make_unique<MockNuclearModel>();
        achilles::ProcessInfo process_info;
        process_info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
        std::vector<achilles::ProcessInfo> infos;
        infos.push_back(process_info);
        infos.back().m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
        infos.push_back(process_info);
        infos.back().m_hadronic = {{achilles::PID::neutron()}, {achilles::PID::neutron()}};
        REQUIRE_CALL(*model, AllowedStates(process_info)).TIMES(1).LR_RETURN(infos);
        REQUIRE_CALL(*model, Mode()).TIMES(AT_LEAST(1)).RETURN(achilles::NuclearMode::Quasielastic);
        REQUIRE_CALL(*model, GetName()).TIMES(AT_LEAST(3)).RETURN("QESpectral");
        REQUIRE_CALL(*model, InspireHEP()).TIMES(AT_LEAST(1)).RETURN("Isaacson:2022cwh");

        auto model_dummy = std::make_unique<MockNuclearModel>();
        auto backend = std::make_unique<MockBackend>();
        trompeloeil::sequence seq;
        REQUIRE_CALL(*backend, SetOptions(config["Backend"]["Options"])).TIMES(1).IN_SEQUENCE(seq);
        REQUIRE_CALL(*backend, AddNuclearModel(trompeloeil::ne(nullptr))).TIMES(1).IN_SEQUENCE(seq);
        REQUIRE_CALL(*backend, SetSherpa(nullptr)).TIMES(1).IN_SEQUENCE(seq);
        REQUIRE_CALL(*backend, Validate()).TIMES(1).IN_SEQUENCE(seq).RETURN(true);
        REQUIRE_CALL(*backend, AddProcess(trompeloeil::_)).TIMES(2).IN_SEQUENCE(seq);
        REQUIRE_CALL(*backend, GetNuclearModel()).TIMES(AT_LEAST(3)).LR_RETURN(model.get());
        MockBackend::SetSelf(std::move(backend));

        auto groups = achilles::ProcessGroup::ConstructGroups(config, model.get(), beam, nucleus);
        std::vector<achilles::ProcessGroup> process_groups;
        for(auto &group : groups) {
            group.second.SetupBackend(config, std::move(model_dummy), nullptr);
            process_groups.push_back(std::move(group.second));
        }
        achilles::NuHepMCWriter writer("data.hepmc", false);
        writer.WriteHeader("dummy", process_groups);

        writer.Write(event);
    }
    auto result = std::system("./bin/NuHepMCReferenceValidator data.hepmc");
    CHECK(result == 0);
    // result = std::system("rm data.hepmc");
}
