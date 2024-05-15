#include "catch2/catch.hpp"
#include "catch2/trompeloeil.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/EventHistory.hh"
#include <string>

class MockVisitor : public trompeloeil::mock_interface<achilles::HistoryVisitor> {
    IMPLEMENT_MOCK1(visit);
};

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
        CHECK(node->IsCascade());
        CHECK(node->ParticlesIn().size() == 0);
        CHECK(node->ParticlesOut().size() == 0);

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
        CHECK(node->ParticlesIn()[0] == part1);

        node = history.FindNodeOut(part2);
        CHECK(node != nullptr);
        CHECK(node->ParticlesOut()[0] == part2);

        CHECK(history.FindNodeIn(part2) == nullptr);
        CHECK(history.FindNodeOut(part1) == nullptr);
    }

    achilles::Particle target(achilles::PID::carbon(),
                              {achilles::ParticleInfo(achilles::PID::carbon()).Mass(), 0, 0, 0});
    achilles::Particle nuc_in(achilles::PID::neutron(), {918.00661011686168, 23.585833577703873,
                                                         85.332385429710143, 52.378992899809319});
    achilles::Particle beam(achilles::PID::nu_muon(), {1e4, 0, 0, 1e4});
    achilles::Particle neutrino(achilles::PID::nu_muon(),
                                {5.3983343748755351e3, 0, 0, 5.3983343748755351e3});
    achilles::Particle nuc_out(
        achilles::PID::proton(),
        {9.6318094613481071e2, -1.8702102417549486e2, -4.7096228265225918e1, 1.0333739302250189e2});
    achilles::Particle lepton(achilles::PID::muon(), {5.3531600388575862e3, 2.1060685775319874e2,
                                                      1.3242861369493605e2, 5.3473759747528429e3});
    // TODO: Figure out how to handle remnant
    // achilles::Particle remnant(achilles::PID(1000060110));
    // remnant.Momentum() = target.Momentum() + neutrino.Momentum() - nuc_out.Momentum() -
    // lepton.Momentum();

    achilles::EventHistory history;
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {target}, {nuc_in}, achilles::EventHistoryNode::StatusCode::target);
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {beam}, {neutrino}, achilles::EventHistoryNode::StatusCode::beam);
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {nuc_in, neutrino}, {nuc_out, lepton},
                      achilles::EventHistoryNode::StatusCode::primary);

    SECTION("Example Event") {
        CHECK(history.size() == 3);
        auto children0 = history.Children(0ul);
        auto children1 = history.Children(1ul);
        CHECK_THAT(children0, Catch::Matchers::Equals(children1));

        auto parents = history.Parents(2ul);
        std::vector<achilles::EventHistoryNode *> expected_parents;
        expected_parents.push_back(history.Node(0));
        expected_parents.push_back(history.Node(1));
        CHECK_THAT(parents, Catch::Matchers::Equals(expected_parents));
    }

    SECTION("Insert Shower Vertex") {
        achilles::Particle neutrino2(achilles::PID::nu_muon(),
                                     {5.6983343748755351e3, 0, 0, 5.6983343748755351e3});
        auto mom = neutrino2.Momentum() - neutrino.Momentum();
        achilles::Particle zd(achilles::PID::Zboson(), mom);
        history.InsertShowerVert(neutrino.Position(), neutrino, neutrino, neutrino2, {zd});
        auto children0 = history.Children(1ul);
        CHECK(children0.size() == 1);
        CHECK(children0[0] == history.Node(3));

        auto children1 = history.Children(3ul);
        auto children2 = history.Children(0ul);
        CHECK_THAT(children1, Catch::Matchers::Equals(children2));
    }

    SECTION("Primary Vertex") {
        CHECK(history.Primary() == history.Node(2));

        history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                          {nuc_in, neutrino}, {nuc_out, lepton},
                          achilles::EventHistoryNode::StatusCode::primary);
        CHECK_THROWS_WITH(history.Primary(), "EventHistory: Only one primary node is allowed!");
    }

    SECTION("Beam Vertex") {
        CHECK(history.Beam() == history.Node(1));

        history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                          {nuc_in, neutrino}, {nuc_out, lepton},
                          achilles::EventHistoryNode::StatusCode::beam);
        CHECK_THROWS_WITH(history.Beam(), "EventHistory: Only one beam node is allowed!");
    }

    SECTION("Target Vertex") {
        CHECK(history.Target() == history.Node(0));

        history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                          {nuc_in, neutrino}, {nuc_out, lepton},
                          achilles::EventHistoryNode::StatusCode::target);
        CHECK_THROWS_WITH(history.Target(), "EventHistory: Only one target node is allowed!");
    }

    SECTION("Find node from pointer") {
        CHECK(history.size() == 3);
        auto children0 = history.Children(history.Node(0));
        auto children1 = history.Children(history.Node(1));
        CHECK_THAT(children0, Catch::Matchers::Equals(children1));

        auto parents = history.Parents(history.Node(2));
        std::vector<achilles::EventHistoryNode *> expected_parents;
        expected_parents.push_back(history.Node(0));
        expected_parents.push_back(history.Node(1));
        CHECK_THAT(parents, Catch::Matchers::Equals(expected_parents));
    }
}

TEST_CASE("EventHistory Visitor", "[EventHistory]") {
    achilles::Particle target(achilles::PID::carbon(),
                              {achilles::ParticleInfo(achilles::PID::carbon()).Mass(), 0, 0, 0});
    achilles::Particle nuc_in(achilles::PID::neutron(), {918.00661011686168, 23.585833577703873,
                                                         85.332385429710143, 52.378992899809319});
    achilles::Particle beam(achilles::PID::nu_muon(), {1e4, 0, 0, 1e4});
    achilles::Particle neutrino(achilles::PID::nu_muon(),
                                {5.3983343748755351e3, 0, 0, 5.3983343748755351e3});
    achilles::Particle nuc_out(
        achilles::PID::proton(),
        {9.6318094613481071e2, -1.8702102417549486e2, -4.7096228265225918e1, 1.0333739302250189e2});
    achilles::Particle lepton(achilles::PID::muon(), {5.3531600388575862e3, 2.1060685775319874e2,
                                                      1.3242861369493605e2, 5.3473759747528429e3});
    // TODO: Figure out how to handle remnant
    // achilles::Particle remnant(achilles::PID(1000060110));
    // remnant.Momentum() = target.Momentum() + neutrino.Momentum() - nuc_out.Momentum() -
    // lepton.Momentum();

    achilles::EventHistory history;
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {target}, {nuc_in}, achilles::EventHistoryNode::StatusCode::target);
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {beam}, {neutrino}, achilles::EventHistoryNode::StatusCode::beam);
    history.AddVertex({8.1502403531109633e-13, 3.6163822359943019e-13, 1.0579315614474801e-12},
                      {nuc_in, neutrino}, {nuc_out, lepton},
                      achilles::EventHistoryNode::StatusCode::primary);

    SECTION("Expected output") {
        achilles::PrintVisitor visitor;
        history.WalkHistory(visitor);

        std::string expected1 =
            R"exp(Node(primary, {Particle[2112, 24, FourVector(9.18006610e+02, 2.35858336e+01, 8.53323854e+01, 5.23789929e+01), ThreeVector(0.00000000e+00, 0.00000000e+00, 0.00000000e+00)], Particle[14, 24, FourVector(5.39833437e+03, 0.00000000e+00, 0.00000000e+00, 5.39833437e+03), ThreeVector(0.00000000e+00, 0.00000000e+00, 0.00000000e+00)]} -> {Particle[2212, 24, FourVector(9.63180946e+02, -1.87021024e+02, -4.70962283e+01, 1.03337393e+02), ThreeVector(0.00000000e+00, 0.00000000e+00, 0.00000000e+00)], Particle[13, 24, FourVector(5.35316004e+03, 2.10606858e+02, 1.32428614e+02, 5.34737597e+03), ThreeVector(0.00000000e+00, 0.00000000e+00, 0.00000000e+00)]}))exp";
        std::string expected2 =
            R"exp(Node(target, {Particle[1000060120, 24, FourVector(1.11880000e+04, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00), ThreeVector(0.00000000e+00, 0.00000000e+00, 0.00000000e+00)]} -> {Particle[2112, 24, FourVector(9.18006610e+02, 2.35858336e+01, 8.53323854e+01, 5.23789929e+01), ThreeVector(0.00000000e+00, 0.00000000e+00, 0.00000000e+00)]}))exp";
        std::string expected3 =
            R"exp(Node(beam, {Particle[14, 24, FourVector(1.00000000e+04, 0.00000000e+00, 0.00000000e+00, 1.00000000e+04), ThreeVector(0.00000000e+00, 0.00000000e+00, 0.00000000e+00)]} -> {Particle[14, 24, FourVector(5.39833437e+03, 0.00000000e+00, 0.00000000e+00, 5.39833437e+03), ThreeVector(0.00000000e+00, 0.00000000e+00, 0.00000000e+00)]}))exp";

        // Check for expected substrings, in case order of visiting changes in future
        CHECK_THAT(visitor.data, Catch::Matchers::Contains(expected1));
        CHECK_THAT(visitor.data, Catch::Matchers::Contains(expected2));
        CHECK_THAT(visitor.data, Catch::Matchers::Contains(expected3));
    }

    SECTION("Mock Visitor") {
        MockVisitor visitor;
        REQUIRE_CALL(visitor, visit(history.Node(0))).TIMES(1);
        REQUIRE_CALL(visitor, visit(history.Node(1))).TIMES(1);
        REQUIRE_CALL(visitor, visit(history.Node(2))).TIMES(1);
        history.WalkHistory(visitor);
    }
}
