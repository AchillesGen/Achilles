#include "catch2/catch.hpp"

#include "Achilles/EventHistory.hh"
#include "Achilles/Constants.hh"

TEST_CASE("EventHistoryNode", "[EventHistory]") {
    achilles::EventHistoryNode node(0);
    achilles::Particle part1(achilles::PID::proton(), {achilles::Constant::mN, 0, 0, 0});
    achilles::Particle part2(achilles::PID::neutron(), {achilles::Constant::mN, 0, 0, 0});
    node.AddIncoming(part1);
    node.AddOutgoing(part2);

    SECTION("HasIncoming") {
        CHECK(node.HasIncoming(part1));
        CHECK(!node.HasIncoming(part2));
    }

    SECTION("HasOutgoing") {
        CHECK(!node.HasOutgoing(part1));
        CHECK(node.HasOutgoing(part2));
    }

    SECTION("Status Checks") {
        CHECK(node.IsCascade());
        CHECK(!node.IsPropagation());
        CHECK(node.IsInternal());
        CHECK(!node.IsExternal());
        CHECK(node.Status() == achilles::EventHistoryNode::StatusCode::cascade);
    }
}

TEST_CASE("EventHistory", "[EventHistory]") {
    SECTION("Node Access") {
        achilles::EventHistory history;
        history.AddVertex({});

        achilles::EventHistoryNode *node = history.Node(0);
        CHECK(node != nullptr);
        CHECK(node -> IsCascade());
        CHECK(node -> ParticlesIn().size() == 0);
        CHECK(node -> ParticlesOut().size() == 0);

        CHECK(history.Node(1) == nullptr);
    }

    SECTION("Find Node") {
        achilles::EventHistory history;
        history.AddVertex({});

        achilles::Particle part1(achilles::PID::proton(), {achilles::Constant::mN, 0, 0, 0});
        achilles::Particle part2(achilles::PID::neutron(), {achilles::Constant::mN, 0, 0, 0});
        history.AddParticleIn(0, part1);
        history.AddParticleOut(0, part2);

        auto *node = history.FindNodeIn(part1);
        CHECK(node != nullptr);
        CHECK(node -> ParticlesIn()[0] == part1);

        node = history.FindNodeOut(part2);
        CHECK(node != nullptr);
        CHECK(node -> ParticlesOut()[0] == part2);

        CHECK(history.FindNodeIn(part2) == nullptr);
        CHECK(history.FindNodeOut(part1) == nullptr);
    }

    achilles::Particle target(achilles::PID::carbon(),
                              {achilles::ParticleInfo(achilles::PID::carbon()).Mass(), 0, 0, 0});
    achilles::Particle nuc_in(achilles::PID::neutron(),
                              {918.00661011686168, 23.585833577703873,
                               85.332385429710143, 52.378992899809319});
    achilles::Particle beam(achilles::PID::nu_muon(),
                            {1e4, 0, 0, 1e4});
    achilles::Particle neutrino(achilles::PID::nu_muon(),
                                {5.3983343748755351e3, 0, 0, 5.3983343748755351e3});
    achilles::Particle nuc_out(achilles::PID::proton(),
                               {9.6318094613481071e2, -1.8702102417549486e2,
                                -4.7096228265225918e1, 1.0333739302250189e2});
    achilles::Particle lepton(achilles::PID::muon(),
                              {5.3531600388575862e3, 2.1060685775319874e2,
                               1.3242861369493605e2, 5.3473759747528429e3});
    // TODO: Figure out how to handle remnant
    // achilles::Particle remnant(achilles::PID(1000060110));
    // remnant.Momentum() = target.Momentum() + neutrino.Momentum() - nuc_out.Momentum() - lepton.Momentum();

    achilles::EventHistory history;
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {target}, {nuc_in}, achilles::EventHistoryNode::StatusCode::target);
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {beam}, {neutrino}, achilles::EventHistoryNode::StatusCode::beam);
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {nuc_in, neutrino}, {nuc_out, lepton}, achilles::EventHistoryNode::StatusCode::primary);

    SECTION("Example Event") {
        CHECK(history.size() == 3);
        auto children0 = history.Children(0ul);
        auto children1 = history.Children(1ul);
        CHECK_THAT(children0, Catch::Matchers::Equals(children1));

        auto parents = history.Parents(2ul);
        std::vector<achilles::EventHistoryNode*> expected_parents;
        expected_parents.push_back(history.Node(0));
        expected_parents.push_back(history.Node(1));
        CHECK_THAT(parents, Catch::Matchers::Equals(expected_parents));
    }

    SECTION("Primary Vertex") {
        CHECK(history.Primary() == history.Node(2));

        history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                          {nuc_in, neutrino}, {nuc_out, lepton}, achilles::EventHistoryNode::StatusCode::primary);
        CHECK_THROWS_WITH(history.Primary(), "EventHistory: Only one primary node is allowed!");
    }

    SECTION("Beam Vertex") {
        CHECK(history.Beam() == history.Node(1));

        history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                          {nuc_in, neutrino}, {nuc_out, lepton}, achilles::EventHistoryNode::StatusCode::beam);
        CHECK_THROWS_WITH(history.Beam(), "EventHistory: Only one beam node is allowed!");
    }

    SECTION("Target Vertex") {
        CHECK(history.Target() == history.Node(0));

        history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                          {nuc_in, neutrino}, {nuc_out, lepton}, achilles::EventHistoryNode::StatusCode::target);
        CHECK_THROWS_WITH(history.Target(), "EventHistory: Only one target node is allowed!");
    }

    SECTION("Find node from pointer") {
        CHECK(history.size() == 3);
        auto children0 = history.Children(history.Node(0));
        auto children1 = history.Children(history.Node(1));
        CHECK_THAT(children0, Catch::Matchers::Equals(children1));

        auto parents = history.Parents(history.Node(2));
        std::vector<achilles::EventHistoryNode*> expected_parents;
        expected_parents.push_back(history.Node(0));
        expected_parents.push_back(history.Node(1));
        CHECK_THAT(parents, Catch::Matchers::Equals(expected_parents));
    }
}
