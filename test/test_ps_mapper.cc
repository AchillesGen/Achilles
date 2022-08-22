#include "catch2/catch.hpp"

#include "Achilles/PhaseSpaceMapper.hh"
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/FinalStateMapper.hh"
#include "mock_classes.hh"
#include "Approx.hh"

class DummyHadron : public achilles::HadronicBeamMapper,
                           achilles::RegistrablePS<achilles::HadronicBeamMapper, DummyHadron, size_t> {
    public:
        DummyHadron(size_t idx) : HadronicBeamMapper(idx, Name()) {}
        static std::string Name() { return "Dummy"; }
        static std::unique_ptr<achilles::HadronicBeamMapper> Construct(const size_t &idx) {
            return std::make_unique<DummyHadron>(idx);
        }

        void GeneratePoint(std::vector<achilles::FourVector> &p, const std::vector<double>&) override {
            p[HadronIdx()] = {achilles::Constant::mN, 0, 0, 0};
        }
        double GenerateWeight(const std::vector<achilles::FourVector>&, std::vector<double>&) override {
            return 1.0;
        }
        size_t NDims() const override { return 0; }
        YAML::Node ToYAML() const override { return YAML::Node(); }
};

class DummyFS : public achilles::FinalStateMapper,
                       achilles::RegistrablePS<achilles::FinalStateMapper, DummyFS, std::vector<double>> {
    public:
        DummyFS(const std::vector<double>&) : FinalStateMapper(2) {}
        static std::string Name() { return "Dummy"; }
        static std::unique_ptr<achilles::FinalStateMapper> Construct(const std::vector<double> &m) {
            return std::make_unique<DummyFS>(m);
        }

        void GeneratePoint(std::vector<achilles::FourVector> &p, const std::vector<double>&) override {
            p[2] = {achilles::Constant::mN, 0, 0, 0};
            p[3] = {achilles::Constant::mN, 0, 0, 0};
        }
        double GenerateWeight(const std::vector<achilles::FourVector>&, std::vector<double>&) override {
            return 1.0;
        }
        YAML::Node ToYAML() const override { return YAML::Node(); }
};

TEST_CASE("PhaseSpaceBuilder", "[PhaseSpace]") {
    auto beam = std::make_shared<MockBeam>();
    achilles::FourVector beam_mom = {achilles::Constant::mN, 0, 0, 0};
    std::vector<achilles::FourVector> expected = {beam_mom, beam_mom, beam_mom, beam_mom};
    std::vector<achilles::FourVector> output(4);
    auto mapper = achilles::PSBuilder(2, 2).Beam(beam, {achilles::Constant::mN2, 0.0}, 1)
                                         .Hadron("Dummy", {achilles::Constant::mN2, 0.0})
                                         .FinalState("Dummy", {achilles::Constant::mN2, 0.0})
                                         .build();
    std::vector<double> rans(2);
    std::vector<double> beam_rans;
    std::set<achilles::PID> beam_id{achilles::PID::electron()};
    REQUIRE_CALL(*beam, BeamIDs())
        .TIMES(1)
        .LR_RETURN((beam_id));
    REQUIRE_CALL(*beam, Flux(achilles::PID::electron(), beam_rans, achilles::Constant::mN2))
        .TIMES(1)
        .LR_RETURN((beam_mom));
    REQUIRE_CALL(*beam, NVariables())
        .TIMES(2)
        .RETURN(0);
    CHECK(mapper -> NDims() == 2);
    mapper -> GeneratePoint(output, rans);
    CHECK_THAT(output, AllFourVectorApprox(expected));
}

TEST_CASE("PhaseSpaceMapper", "[PhaseSpace]") {
    auto beam_map = std::make_shared<MockMapper>();
    auto hadron_map = std::make_shared<MockMapper>();
    auto final_state_map = std::make_unique<MockMapper>();
    std::vector<double> masses{0, achilles::Constant::mN2};
    size_t beam_dims = 0, hadron_dims = 4, final_state_dims = 2;
    static constexpr double sqrts = 200;
    achilles::FourVector beam_mom = {100, 0, 0, -100};
    achilles::FourVector hadron_mom = {100, 0, 0, 100};
    std::vector<achilles::FourVector> final_state_mom ={{sqrts/2, sqrts/2*sin(M_PI/4), 0, sqrts/2*cos(M_PI/4)},
                                                      {sqrts/2, -sqrts/2*sin(M_PI/4), 0, -sqrts/2*cos(M_PI/4)}}; 
    std::vector<achilles::FourVector> mom_out = {hadron_mom, beam_mom, final_state_mom[0], final_state_mom[1]};
    std::vector<achilles::FourVector> mom(4), mom2(4), mom3(4);
    mom2[1] = beam_mom;
    mom3[0] = hadron_mom;
    mom3[1] = beam_mom;
    std::vector<double> lbeam_rans{}, lbeam_rans_out{};
    std::vector<double> hadron_rans{0.5, 0.5, 0.5, 0.5}, hadron_rans_out(4);
    std::vector<double> final_state_rans{0.5, 0.5}, final_state_rans_out(2);
    std::vector<double> rans{0.5, 0.5, 0.5, 0.5, 0.5, 0.5}, rans_out(6);

    REQUIRE_CALL(*beam_map, NDims())
        .TIMES(2)
        .LR_RETURN((beam_dims));
    REQUIRE_CALL(*beam_map, GeneratePoint(mom, lbeam_rans))
        .TIMES(1)
        .LR_SIDE_EFFECT(mom[1] = beam_mom);
    REQUIRE_CALL(*beam_map, GenerateWeight(mom_out, lbeam_rans_out))
        .TIMES(1)
        .LR_SIDE_EFFECT(_2 = lbeam_rans)
        .RETURN(1.0);

    REQUIRE_CALL(*hadron_map, NDims())
        .TIMES(2)
        .LR_RETURN((hadron_dims));
    REQUIRE_CALL(*hadron_map, GeneratePoint(mom2, hadron_rans))
        .TIMES(1)
        .LR_SIDE_EFFECT(mom[0] = hadron_mom);
    REQUIRE_CALL(*hadron_map, GenerateWeight(mom_out, hadron_rans_out))
        .TIMES(1)
        .LR_SIDE_EFFECT(_2 = hadron_rans)
        .RETURN(1.0);

    REQUIRE_CALL(*final_state_map, NDims())
        .TIMES(2)
        .LR_RETURN((final_state_dims));
    REQUIRE_CALL(*final_state_map, GeneratePoint(mom3, final_state_rans))
        .TIMES(1)
        .LR_SIDE_EFFECT(mom[2] = final_state_mom[0])
        .LR_SIDE_EFFECT(mom[3] = final_state_mom[1]);
    REQUIRE_CALL(*final_state_map, GenerateWeight(mom_out, final_state_rans_out))
        .TIMES(1)
        .LR_SIDE_EFFECT(_2 = final_state_rans)
        .RETURN(1.0);

    achilles::PSMapper mapper(2, 2);
    mapper.SetLeptonBeam(beam_map);
    mapper.SetHadronBeam(hadron_map);
    mapper.SetFinalState(std::move(final_state_map));

    mapper.GeneratePoint(mom, rans);
    auto wgt = mapper.GenerateWeight(mom, rans_out);

    CHECK(wgt == 1.0);
    CHECK_THAT(rans, Catch::Matchers::Approx(rans_out));
    CHECK_THAT(mom, AllFourVectorApprox(mom_out));

    // TODO: Move this to test EventGen.cc
    // MockPSBuilder builder; 
    // trompeloeil::sequence seq;
    // REQUIRE_CALL(builder, Beam(beam, 0UL))
    //     .TIMES(1)
    //     .IN_SEQUENCE(seq)
    //     .LR_SIDE_EFFECT(mapper -> LeptonBeam() = beam_map)
    //     .LR_RETURN((builder));
    // REQUIRE_CALL(builder, Hadron("Dummy", masses, 1UL))
    //     .TIMES(1)
    //     .IN_SEQUENCE(seq)
    //     .LR_SIDE_EFFECT(mapper -> Hadron() = hadron_map)
    //     .LR_RETURN((builder));
    // REQUIRE_CALL(builder, FinalState("Dummy", masses))
    //     .TIMES(1)
    //     .IN_SEQUENCE(seq)
    //     .LR_SIDE_EFFECT(mapper -> SetFinalState(std::move(final_state_map)))
    //     .LR_RETURN((builder));
}
