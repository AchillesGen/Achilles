#include "catch2/catch.hpp"

#include "Achilles/BeamMapper.hh"
#include "Approx.hh"
#include "mock_classes.hh"

TEST_CASE("BeamMapper", "[PhaseSpace]") {
    static constexpr auto electron = achilles::PID::electron();
    const std::set<achilles::PID> beam_ids{electron};
    const achilles::FourVector beam_mom{1000, 0, 0, 1000};
    const std::vector<double> beam_rans{0.5};
    std::vector<double> new_rans(1);
    static constexpr int nvars = 1;

    SECTION("Variables is passed through") {
        auto beam = std::make_shared<MockBeam>();
        REQUIRE_CALL(*beam, NVariables()).TIMES(1).LR_RETURN((nvars));
        achilles::BeamMapper mapper(0, beam);
        CHECK(mapper.NDims() == nvars);
    }

    SECTION("Forward and Backward pass") {
        auto beam = std::make_shared<MockBeam>();
        REQUIRE_CALL(*beam, BeamIDs()).TIMES(2).LR_RETURN((beam_ids));
        REQUIRE_CALL(*beam, Flux(electron, beam_rans, trompeloeil::ge(0)))
            .TIMES(1)
            .LR_RETURN((beam_mom));
        REQUIRE_CALL(*beam, GenerateWeight(electron, beam_mom, trompeloeil::_, trompeloeil::ge(0)))
            .LR_SIDE_EFFECT(_3[0] = 0.5)
            .TIMES(1)
            .RETURN(1.0);

        SECTION("Forward") {
            achilles::BeamMapper mapper(0, beam);
            mapper.SetMasses({0});
            std::vector<achilles::FourVector> mom(1);
            mapper.GeneratePoint(mom, beam_rans);
            double wgt = mapper.GenerateWeight(mom, new_rans);
            CHECK_THAT(beam_rans, Catch::Matchers::Approx(new_rans));
            CHECK(wgt == 1.0);
            CHECK_THAT(mom, AllFourVectorApprox({beam_mom}));
        }

        SECTION("Reverse") {
            achilles::BeamMapper mapper(0, beam);
            mapper.SetMasses({0});
            double wgt = mapper.GenerateWeight({beam_mom}, new_rans);
            std::vector<achilles::FourVector> mom(1);
            mapper.GeneratePoint(mom, new_rans);
            CHECK_THAT(beam_rans, Catch::Matchers::Approx(new_rans));
            CHECK(wgt == 1.0);
            CHECK_THAT(mom, AllFourVectorApprox({beam_mom}));
        }
    }
}
