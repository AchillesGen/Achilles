#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include "Achilles/Beams.hh"
#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ProcessInfo.hh"

TEST_CASE("Initialize Event Parameters", "[Event]") {
    auto nuc = std::make_shared<MockNucleus>();
    auto beam = std::make_shared<MockBeam>();
    std::vector<double> rans{0};
    static constexpr achilles::FourVector lepton0{1000, 0, 0, 1000};
    static constexpr achilles::FourVector lepton1{313.073, 105.356, 174.207, -237.838};
    static constexpr achilles::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
    static constexpr achilles::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};

    std::vector<achilles::FourVector> moms = {hadron0, lepton0, lepton1, hadron1};
    achilles::Particles particles = {{achilles::PID::proton(), hadron0}};

    REQUIRE_CALL(*nuc, GenerateConfig())
        .TIMES(1);
    REQUIRE_CALL(*nuc, NNucleons())
        .LR_RETURN((12UL))
        .TIMES(1);
    static constexpr double vegas_wgt = 10;
    achilles::Event event(nuc, moms, vegas_wgt);

    SECTION("Nucleus and Beam set correctly") {
        CHECK(event.Momentum()[0] == hadron0);
        CHECK(event.Momentum()[1] == lepton0);
        CHECK(event.Momentum()[2] == lepton1);
        CHECK(event.Momentum()[3] == hadron1);
    }

    SECTION("Initialize Particles") {
        REQUIRE_CALL(*nuc, Nucleons())
            .LR_RETURN((particles))
            .TIMES(5);

        achilles::Process_Info info;
        info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
        info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
        info.m_mom_map = {{0, achilles::PID::proton().AsInt()},
                          {1, achilles::PID::electron().AsInt()},
                          {2, achilles::PID::electron().AsInt()},
                          {3, achilles::PID::proton().AsInt()}};
        event.InitializeLeptons(info);
        event.InitializeHadrons(info);

        auto leptons = event.Leptons();
        CHECK(leptons.size() == 2);
        CHECK(leptons[0].Momentum() == moms[1]); 
        CHECK(leptons[1].Momentum() == moms[2]); 
        CHECK(leptons[0].ID() == achilles::PID::electron());
        CHECK(leptons[1].ID() == achilles::PID::electron());
        CHECK(leptons[0].Status() == achilles::ParticleStatus::initial_state);
        CHECK(leptons[1].Status() == achilles::ParticleStatus::final_state);

        auto hadrons = event.Hadrons();
        CHECK(hadrons.size() == 2);
        CHECK(hadrons[0].Momentum() == moms[0]); 
        CHECK(hadrons[1].Momentum() == moms[3]); 
        CHECK(hadrons[0].ID() == achilles::PID::proton());
        CHECK(hadrons[1].ID() == achilles::PID::proton());
        CHECK(hadrons[0].Status() == achilles::ParticleStatus::initial_state);
        CHECK(hadrons[1].Status() == achilles::ParticleStatus::propagating);

        auto output = event.Particles();
        CHECK(output[0] == hadrons[0]);
        CHECK(output[1] == hadrons[1]);
        CHECK(output[2] == leptons[0]);
        CHECK(output[3] == leptons[1]);
    }

    SECTION("Weight is correct") {
        for(auto &me : event.MatrixElementWgts()) {
            me = 10;
        }

        event.TotalCrossSection();
        event.CalcWeight();
        CHECK(event.Weight() == 1200);
        
        // TODO: Rewrite this to match new framework
        // auto probs = event.EventProbs();
        // CHECK(probs.size() == 11);
        // size_t idx = 0;
        // static constexpr double base_prob = 0.1;
        // for(const auto &prob : probs) {
        //     CHECK(prob == Approx(base_prob*static_cast<double>(idx++)));
        // }
    }

    SECTION("Event can be finalized") {
        // Dummy carbon event
        achilles::Particles final = {{achilles::PID::proton(), hadron0, {}, achilles::ParticleStatus::initial_state},
                                     {achilles::PID::proton(), hadron1, {}, achilles::ParticleStatus::final_state},
                                     {achilles::PID::proton(), hadron0},
                                     {achilles::PID::proton(), hadron0},
                                     {achilles::PID::proton(), hadron0},
                                     {achilles::PID::proton(), hadron0},
                                     {achilles::PID::proton(), hadron0},
                                     {achilles::PID::neutron(), hadron0},
                                     {achilles::PID::neutron(), hadron0},
                                     {achilles::PID::neutron(), hadron0},
                                     {achilles::PID::neutron(), hadron0},
                                     {achilles::PID::neutron(), hadron0},
                                     {achilles::PID::neutron(), hadron0}};
        REQUIRE_CALL(*nuc, Nucleons())
            .LR_RETURN((final))
            .TIMES(AT_LEAST(13));

        event.Finalize();
        CHECK(event.Remnant().PID() == 1000050110);
        CHECK(event.Remnant().Mass() == 11*achilles::Constant::mN);
    }
}
