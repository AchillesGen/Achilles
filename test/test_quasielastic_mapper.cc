#include "catch2/catch.hpp"

#include "Achilles/QuasielasticTestMapper.hh"
#include "mock_classes.hh"
#include "Approx.hh"

TEST_CASE("QuasielasticTestMapper", "[PhaseSpace]") {
    static constexpr auto electron = achilles::PID::electron();
    const std::set<achilles::PID> beam_ids = {electron};
    const achilles::FourVector beam_mom = {1000, 0, 0, 1000};
    const std::vector<double> beam_rans{0};
    static constexpr int nvars = 1;
    auto beam = std::make_shared<MockBeam>();  
    REQUIRE_CALL(*beam, BeamIDs())
        .TIMES(2)
        .LR_RETURN((beam_ids));
    REQUIRE_CALL(*beam, Flux(electron, beam_rans, 0))
        .TIMES(1)
        .LR_RETURN((beam_mom));
    REQUIRE_CALL(*beam, GenerateWeight(electron, beam_mom, trompeloeil::_, 0))
        .TIMES(1)
        .RETURN(1.0);
    REQUIRE_CALL(*beam, NVariables())
        .TIMES(4)
        .RETURN(nvars);

    auto input_fmt = R"node(
        Run Mode: {}
        Angle: 90
        Energy: 500
    )node";
    auto mode = GENERATE(table<std::string, size_t>({
                 {"FixedAngleEnergy", 4},
                 {"FixedAngle",       5},
                 {"FullPhaseSpace",   6}}));
    auto input = fmt::format(input_fmt, std::get<0>(mode));
    auto node = YAML::Load(input);
    std::vector<double> rans(std::get<1>(mode)+nvars, 0.5);
    if(rans.size() == 0)
        throw;
    rans[0] = 0;
    achilles::QuasielasticTestMapper mapper(node, beam);

    // TODO: How to validate wgt is right
    SECTION("Forward Map") {
        std::vector<achilles::FourVector> mom(4);
        mapper.GeneratePoint(mom, rans);
        std::vector<double> rans_out(rans.size());
        mapper.GenerateWeight(mom, rans_out);
        CHECK_THAT(rans, Catch::Matchers::Approx(rans_out));
        CHECK_THAT(mom[0].E() + mom[1].E(),
                   Catch::Matchers::WithinRel(mom[2].E() + mom[3].E()));
        CHECK_THAT(mom[0].Px() + mom[1].Px(),
                   Catch::Matchers::WithinRel(mom[2].Px() + mom[3].Px()));
        CHECK_THAT(mom[0].Py() + mom[1].Py(),
                   Catch::Matchers::WithinRel(mom[2].Py() + mom[3].Py()));
        CHECK_THAT(mom[0].Pz() + mom[1].Pz(),
                   Catch::Matchers::WithinRel(mom[2].Pz() + mom[3].Pz()));
    }

    // TODO: Figure out how to handle this, since each mode needs different momenta
    // SECTION("Reverse Map") {
    //     achilles::FourVector nucleon_in{6.414e+02, 1.284e+02, -6.791e+01, 2.630e+02};
    //     achilles::FourVector nucleon_out{1.377e+03, 1.193e+02, -6.882e+01, 9.991e+02};
    //     achilles::FourVector beam_out{2.641e+02, 9.129e+00, 9.111e-01, 2.639e+02};

    //     std::vector<achilles::FourVector> mom{nucleon_in, beam_mom, nucleon_out, beam_out};
    //     double wgt = mapper.GenerateWeight(mom, rans);
    //     fmt::print("{}\n", fmt::join(rans, ", "));
    //     std::vector<achilles::FourVector> moms_out(mom.size());
    //     mapper.GeneratePoint(moms_out, rans);
    //     CHECK_THAT(mom, AllFourVectorApprox(moms_out));
    // }
}
