#include "catch2/catch.hpp"
#include "mock_classes.hh"

#include "nuchic/Beams.hh"
#include "nuchic/Event.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"

TEST_CASE("Initialize Event Parameters", "[Event]") {
    auto nuc = std::make_shared<MockNucleus>();
    auto beam = std::make_shared<MockBeam>();
    std::vector<double> rans{0};
    static constexpr nuchic::FourVector lepton0{0, 0, 1000, 1000};
    static constexpr nuchic::FourVector lepton1{105.356, 174.207, -237.838, 313.073};
    static constexpr nuchic::FourVector hadron0{26.8702, -30.5306, -10.9449, 65.4247};
    static constexpr nuchic::FourVector hadron1{-78.4858, -204.738, 1226.89, 1560.42};

    std::vector<nuchic::FourVector> moms = {hadron0, lepton0, lepton1, hadron1};
    nuchic::Particles particles = {{nuchic::PID::proton(), hadron0}};

    REQUIRE_CALL(*nuc, GenerateConfig())
        .TIMES(1);
    nuchic::Event event(nuc, moms, 10);

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

        // TODO: Rewrite this to match new framework
        // event.MatrixElements().resize(1);
        // event.MatrixElement(0).inital_state.emplace_back(nuchic::PID::electron());
        // event.MatrixElement(0).inital_state.emplace_back(nuchic::PID::proton());
        // event.MatrixElement(0).final_state.emplace_back(nuchic::PID::electron());
        // event.MatrixElement(0).final_state.emplace_back(nuchic::PID::proton());

        // event.InitializeLeptons(0);
        // CHECK(event.Leptons().size() == 2);
        // CHECK(event.Leptons()[0].ID() == nuchic::PID::electron());
        // CHECK(event.Leptons()[0].Momentum() == lepton0);
        // CHECK(event.Leptons()[1].ID() == nuchic::PID::electron());
        // CHECK(event.Leptons()[1].Momentum() == lepton1);

        // event.InitializeHadrons({{0, 2, 3}});
        // auto hadrons = event.Hadrons();
        // CHECK(hadrons.size() == 2);
        // CHECK(hadrons[0].ID() == nuchic::PID::proton());
        // CHECK(hadrons[0].Momentum() == hadron0);
        // CHECK(hadrons[1].ID() == nuchic::PID::proton());
        // CHECK(hadrons[1].Momentum() == hadron1);
    }

    SECTION("Weight is correct") {
        event.MatrixElements().resize(10);
        for(size_t i = 0; i < 10; ++i) {
            event.MatrixElement(i).weight = 10;
            event.MatrixElement(i).inital_state.emplace_back(nuchic::PID::proton());
        }

        event.TotalCrossSection();

        CHECK(event.Weight() == 100);
        
        // TODO: Rewrite this to match new framework
        // auto probs = event.EventProbs();
        // CHECK(probs.size() == 11);
        // size_t idx = 0;
        // static constexpr double base_prob = 0.1;
        // for(const auto &prob : probs) {
        //     CHECK(prob == Approx(base_prob*static_cast<double>(idx++)));
        // }
    }
}
