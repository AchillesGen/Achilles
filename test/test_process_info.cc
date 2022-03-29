#include "catch2/catch.hpp"

#include "nuchic/ProcessInfo.hh"

TEST_CASE("Multiplicities", "[ProcessInfo]") {
    nuchic::Process_Info info;     
    info.m_ids = {nuchic::PID::electron(), nuchic::PID::electron()};
    info.m_states = {{{nuchic::PID::neutron()}, {nuchic::PID::neutron()}}};
    CHECK(info.Multiplicity() == 4);
}

TEST_CASE("Masses", "[ProcessInfo]") {
    nuchic::Process_Info info;     
    info.m_ids = {nuchic::PID::electron(), nuchic::PID::electron()};
    info.m_states = {{{nuchic::PID::neutron()}, {nuchic::PID::neutron()}}};
    CHECK(info.Masses() == std::vector<double>{pow(nuchic::ParticleInfo(nuchic::PID::neutron()).Mass(), 2),
                                               0});
}

TEST_CASE("IDs", "[ProcessInfo]") {
    SECTION("Single nucleon initial state, single nucleon final state") {
        nuchic::Process_Info info;     
        info.m_ids = {nuchic::PID::electron(), nuchic::PID::electron()};
        info.m_states = {{{nuchic::PID::neutron()}, {nuchic::PID::neutron()}}};
        CHECK(info.Ids() == std::vector<long>{nuchic::PID::neutron().AsInt(),
                                              nuchic::PID::electron().AsInt(),
                                              nuchic::PID::neutron().AsInt(),
                                              nuchic::PID::electron().AsInt()});
    }

    SECTION("Single nucleon initial state, multiple nucleon final state") {
        nuchic::Process_Info info;     
        info.m_ids = {nuchic::PID::electron(), nuchic::PID::electron()};
        info.m_states = {{{nuchic::PID::neutron()}, {nuchic::PID::neutron(), nuchic::PID::pion0()}}};
        CHECK(info.Ids() == std::vector<long>{nuchic::PID::neutron().AsInt(),
                                              nuchic::PID::electron().AsInt(),
                                              nuchic::PID::neutron().AsInt(),
                                              nuchic::PID::pion0().AsInt(),
                                              nuchic::PID::electron().AsInt()});
    }

    SECTION("Multiple nucleon initial state, multiple nucleon final state") {
        nuchic::Process_Info info;     
        info.m_ids = {nuchic::PID::electron(), nuchic::PID::electron()};
        info.m_states = {{{nuchic::PID::neutron(), nuchic::PID::neutron()},
                          {nuchic::PID::neutron(), nuchic::PID::neutron()}}};
        CHECK(info.Ids() == std::vector<long>{nuchic::PID::MECnn(),
                                              nuchic::PID::electron().AsInt(),
                                              nuchic::PID::neutron().AsInt(),
                                              nuchic::PID::neutron().AsInt(),
                                              nuchic::PID::electron().AsInt()});
    }
}

TEST_CASE("YAML Decoding", "[ProcessInfo]") {
    YAML::Node node;
    node["Model"] = "SM";
    node["Final States"][0] = 11;
    node["Final States"][1] = 12;
    node["Final States"][2] = 13;

    auto info = node.as<nuchic::Process_Info>();
    CHECK(info.m_ids == std::vector<nuchic::PID>{nuchic::PID::electron(),
                                                 nuchic::PID::nu_electron(),
                                                 nuchic::PID::muon()});
    CHECK(info.m_model == "SM");
}
