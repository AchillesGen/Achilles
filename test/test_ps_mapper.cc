#include "catch2/catch.hpp"

#include "Achilles/FinalStateMapper.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/PhaseSpaceMapper.hh"
#include "Achilles/Event.hh"
#include "Approx.hh"
#include "mock_classes.hh"

class DummyHadron : public achilles::HadronicBeamMapper,
                    achilles::Registrable<achilles::HadronicBeamMapper, DummyHadron,
                                          const achilles::ProcessInfo &> {
  public:
    DummyHadron(const achilles::ProcessInfo &info) : HadronicBeamMapper(info) {}
    static std::string Name() { return "DummyHadron"; }
    static std::unique_ptr<achilles::HadronicBeamMapper>
    Construct(const achilles::ProcessInfo &info) {
        return std::make_unique<DummyHadron>(info);
    }

    void GeneratePoint(achilles::Event& event, const std::vector<double> &) override {
        event.addHadronIn({achilles::Constant::mN, 0, 0, 0});
    }
    double GenerateWeight(const achilles::Event&,
                          std::vector<double> &) override {
        return 1.0;
    }
    size_t NDims() const override { return 0; }
};

class DummyFS
    : public achilles::FinalStateMapper,
      achilles::Registrable<achilles::FinalStateMapper, DummyFS, std::vector<double>> {
  public:
    DummyFS(const std::vector<double> &) : FinalStateMapper(2) {}
    static std::string Name() { return "DummyFS"; }
    static std::unique_ptr<achilles::FinalStateMapper> Construct(const std::vector<double> &m) {
        return std::make_unique<DummyFS>(m);
    }

    void GeneratePoint(achilles::Event& event, const std::vector<double> &) override {
        event.addHadronOut({achilles::Constant::mN, 0, 0, 0});
        event.addHadronOut({achilles::Constant::mN, 0, 0, 0});
    }
    double GenerateWeight(const achilles::Event&,
                          std::vector<double> &) override {
        return 1.0;
    }
};

TEST_CASE("PhaseSpaceBuilder", "[PhaseSpace]") {
    achilles::ProcessInfo info({achilles::PID::electron(), {achilles::PID::electron()}});
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};

    auto beam = std::make_shared<MockBeam>();
    achilles::FourVector beam_mom = {achilles::Constant::mN, 0, 0, 0};
    std::vector<achilles::FourVector> expected = {beam_mom, beam_mom, beam_mom, beam_mom};
    std::vector<achilles::FourVector> output(4);
	achilles::Event event(output);
    auto mapper =
        achilles::PSBuilder(info).Beam(beam, 0).Hadron("DummyHadron").FinalState("DummyFS").build();
    std::vector<double> rans(2);
    std::vector<double> beam_rans;
    std::set<achilles::PID> beam_id{achilles::PID::electron()};
    REQUIRE_CALL(*beam, BeamIDs()).TIMES(1).LR_RETURN((beam_id));
    REQUIRE_CALL(*beam, Flux(achilles::PID::electron(), beam_rans, trompeloeil::ge(0)))
        .TIMES(1)
        .LR_RETURN((beam_mom));
    REQUIRE_CALL(*beam, NVariables()).TIMES(2).RETURN(0);
    CHECK(mapper->NDims() == 2);
    mapper->GeneratePoint(event, rans);
    CHECK_THAT(event.Momentum(), AllFourVectorApprox(expected));
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
    std::vector<achilles::FourVector> final_state_mom = {
        {sqrts / 2, sqrts / 2 * sin(M_PI / 4), 0, sqrts / 2 * cos(M_PI / 4)},
        {sqrts / 2, -sqrts / 2 * sin(M_PI / 4), 0, -sqrts / 2 * cos(M_PI / 4)}};
    std::vector<achilles::FourVector> mom_out = {hadron_mom, beam_mom, final_state_mom[0],
                                                 final_state_mom[1]};
    std::vector<achilles::FourVector> mom(4), mom2(4), mom3(4);
    mom2[1] = beam_mom;
    mom3[0] = hadron_mom;
    mom3[1] = beam_mom;
    std::vector<double> lbeam_rans{}, lbeam_rans_out{};
    std::vector<double> hadron_rans{0.5, 0.5, 0.5, 0.5}, hadron_rans_out(4);
    std::vector<double> final_state_rans{0.5, 0.5}, final_state_rans_out(2);
    std::vector<double> rans{0.5, 0.5, 0.5, 0.5, 0.5, 0.5}, rans_out(6);

	achilles::Event final_state_event(final_state_mom),event_out(mom_out),
					event(mom),event2(mom2),event3(mom3);

    REQUIRE_CALL(*beam_map, NDims()).TIMES(2).LR_RETURN((beam_dims));
    REQUIRE_CALL(*beam_map, GeneratePoint(event, lbeam_rans))
        .TIMES(1)
        .LR_SIDE_EFFECT(event.Momentum()[1] = beam_mom);
    REQUIRE_CALL(*beam_map, GenerateWeight(event_out, lbeam_rans_out))
        .TIMES(1)
        .LR_SIDE_EFFECT(_2 = lbeam_rans)
        .RETURN(1.0);

    REQUIRE_CALL(*hadron_map, NDims()).TIMES(2).LR_RETURN((hadron_dims));
    REQUIRE_CALL(*hadron_map, GeneratePoint(event2, hadron_rans))
        .TIMES(1)
        .LR_SIDE_EFFECT(event.Momentum()[0] = hadron_mom);
    REQUIRE_CALL(*hadron_map, GenerateWeight(event_out, hadron_rans_out))
        .TIMES(1)
        .LR_SIDE_EFFECT(_2 = hadron_rans)
        .RETURN(1.0);

    REQUIRE_CALL(*final_state_map, NDims()).TIMES(2).LR_RETURN((final_state_dims));
    REQUIRE_CALL(*final_state_map, GeneratePoint(event3, final_state_rans))
        .TIMES(1)
        .LR_SIDE_EFFECT(event.Momentum()[2] = final_state_mom[0])
        .LR_SIDE_EFFECT(event.Momentum()[3] = final_state_mom[1]);
    REQUIRE_CALL(*final_state_map, GenerateWeight(event_out, final_state_rans_out))
        .TIMES(1)
        .LR_SIDE_EFFECT(_2 = final_state_rans)
        .RETURN(1.0);

    achilles::PSMapper mapper(2, 2, 0);
    mapper.SetLeptonBeam(beam_map);
    mapper.SetHadronBeam(hadron_map);
    mapper.SetFinalState(std::move(final_state_map));

    mapper.GeneratePoint(event, rans);
    auto wgt = mapper.GenerateWeight(event, rans_out);

    CHECK(wgt == 1.0);
    CHECK_THAT(rans, Catch::Matchers::Approx(rans_out));
    CHECK_THAT(event.Momentum(), AllFourVectorApprox(event_out.Momentum()));

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
