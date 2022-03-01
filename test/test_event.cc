#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include "nuchic/Beams.hh"
#include "nuchic/Event.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ProcessInfo.hh"

TEST_CASE("Initialize Event Parameters", "[Event]") {
    auto nuc = std::make_shared<MockNucleus>();
    auto beam = std::make_shared<MockBeam>();
    std::vector<double> rans{0};
    static constexpr nuchic::FourVector lepton0{1000, 0, 0, 1000};
    static constexpr nuchic::FourVector lepton1{313.073, 105.356, 174.207, -237.838};
    static constexpr nuchic::FourVector hadron0{65.4247, 26.8702, -30.5306, -10.9449};
    static constexpr nuchic::FourVector hadron1{1560.42, -78.4858, -204.738, 1226.89};

    std::vector<nuchic::FourVector> moms = {hadron0, lepton0, lepton1, hadron1};
    nuchic::Particles particles = {{nuchic::PID::proton(), hadron0}};

    REQUIRE_CALL(*nuc, GenerateConfig())
        .TIMES(1);
    REQUIRE_CALL(*nuc, NNucleons())
        .LR_RETURN((12UL))
        .TIMES(1);
    static constexpr double vegas_wgt = 10;
    nuchic::Event event(nuc, moms, vegas_wgt);

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

        nuchic::Process_Info info;
        info.m_ids = {nuchic::PID::electron(), nuchic::PID::electron()};
        info.m_states = {{{nuchic::PID::proton()}, {nuchic::PID::proton()}}};
        info.m_mom_map = {{0, nuchic::PID::proton().AsInt()},
                          {1, nuchic::PID::electron().AsInt()},
                          {2, nuchic::PID::electron().AsInt()},
                          {3, nuchic::PID::proton().AsInt()}};
        event.InitializeLeptons(info);
        event.InitializeHadrons(info);

        auto leptons = event.Leptons();
        CHECK(leptons.size() == 2);
        CHECK(leptons[0].Momentum() == moms[1]); 
        CHECK(leptons[1].Momentum() == moms[2]); 
        CHECK(leptons[0].ID() == nuchic::PID::electron());
        CHECK(leptons[1].ID() == nuchic::PID::electron());
        CHECK(leptons[0].Status() == nuchic::ParticleStatus::initial_state);
        CHECK(leptons[1].Status() == nuchic::ParticleStatus::final_state);

        auto hadrons = event.Hadrons();
        CHECK(hadrons.size() == 2);
        CHECK(hadrons[0].Momentum() == moms[0]); 
        CHECK(hadrons[1].Momentum() == moms[3]); 
        CHECK(hadrons[0].ID() == nuchic::PID::proton());
        CHECK(hadrons[1].ID() == nuchic::PID::proton());
        CHECK(hadrons[0].Status() == nuchic::ParticleStatus::initial_state);
        CHECK(hadrons[1].Status() == nuchic::ParticleStatus::propagating);

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
        nuchic::Particles final = {{nuchic::PID::proton(), hadron0, {}, nuchic::ParticleStatus::initial_state},
                                   {nuchic::PID::proton(), hadron1, {}, nuchic::ParticleStatus::final_state},
                                   {nuchic::PID::proton(), hadron0},
                                   {nuchic::PID::proton(), hadron0},
                                   {nuchic::PID::proton(), hadron0},
                                   {nuchic::PID::proton(), hadron0},
                                   {nuchic::PID::proton(), hadron0},
                                   {nuchic::PID::neutron(), hadron0},
                                   {nuchic::PID::neutron(), hadron0},
                                   {nuchic::PID::neutron(), hadron0},
                                   {nuchic::PID::neutron(), hadron0},
                                   {nuchic::PID::neutron(), hadron0},
                                   {nuchic::PID::neutron(), hadron0}};
        REQUIRE_CALL(*nuc, Nucleons())
            .LR_RETURN((final))
            .TIMES(AT_LEAST(13));

        event.Finalize();
        CHECK(event.Remnant().PID() == 1000050110);
        CHECK(event.Remnant().Mass() == 11*nuchic::Constant::mN);
    }
}
