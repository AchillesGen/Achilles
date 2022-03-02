#include "catch2/catch.hpp"

#include "nuchic/BeamMapper.hh"
#include "mock_classes.hh"
#include "Approx.hh"

TEST_CASE("BeamMapper", "[PhaseSpace]") {
    static constexpr auto electron = nuchic::PID::electron();
    const std::set<nuchic::PID> beam_ids{electron};
    const nuchic::FourVector beam_mom{1000, 0, 0, 1000};
    const std::vector<double> beam_rans{0.5};
    std::vector<double> new_rans(1);
    static constexpr int nvars = 1;

    SECTION("Variables is passed through") {
        auto beam = std::make_shared<MockBeam>();
        REQUIRE_CALL(*beam, NVariables())
            .TIMES(1)
            .LR_RETURN((nvars));
        nuchic::BeamMapper mapper(0, beam);
        CHECK(mapper.NDims() == nvars);
    }

    SECTION("Forward and Backward pass") {
        auto beam = std::make_shared<MockBeam>();
        REQUIRE_CALL(*beam, BeamIDs())
            .TIMES(2)
            .LR_RETURN((beam_ids));
        REQUIRE_CALL(*beam, Flux(electron, beam_rans))
            .TIMES(1)
            .LR_RETURN((beam_mom));
        REQUIRE_CALL(*beam, GenerateWeight(electron, beam_mom, trompeloeil::_))
            .LR_SIDE_EFFECT(_3[0] = 0.5)
            .TIMES(1)
            .RETURN(1.0);

        SECTION("Forward") {
            nuchic::BeamMapper mapper(0, beam);
            std::vector<nuchic::FourVector> mom(1);
            mapper.GeneratePoint(mom, beam_rans);
            double wgt = mapper.GenerateWeight(mom, new_rans);
            CHECK_THAT(beam_rans, 
                       Catch::Matchers::Approx(new_rans));
            CHECK(wgt == 1.0);
            CHECK_THAT(mom, AllFourVectorApprox({beam_mom}));
        }

        SECTION("Reverse") {
            nuchic::BeamMapper mapper(0, beam);
            double wgt = mapper.GenerateWeight({beam_mom}, new_rans);
            std::vector<nuchic::FourVector> mom(1);
            mapper.GeneratePoint(mom, new_rans);
            CHECK_THAT(beam_rans, 
                       Catch::Matchers::Approx(new_rans));
            CHECK(wgt == 1.0);
            CHECK_THAT(mom, AllFourVectorApprox({beam_mom}));
        }
    }
}
