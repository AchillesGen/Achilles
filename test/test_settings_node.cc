#include "catch2/catch.hpp"

#include "Achilles/Settings.hh"

#include <filesystem>
#include <fstream>

// The YAML-node Settings constructor (used by the geometry driver to forward
// the "Achilles" section of a run card) must behave like the file-based one,
// including resolution of nested !include tags.
TEST_CASE("Settings from a YAML node", "[Settings]") {
    SECTION("plain node") {
        auto node = YAML::Load("Main:\n  NEvents: 7\n");
        achilles::Settings settings(node);
        CHECK(settings.GetAs<int>("Main/NEvents") == 7);
    }

    SECTION("nested !include tags are resolved") {
        const auto inc =
            (std::filesystem::temp_directory_path() / "achilles_settings_include.yml").string();
        {
            std::ofstream out(inc);
            out << "Value: 42\n";
        }
        auto node = YAML::Load("Root:\n  Sub: !include \"" + inc + "\"\n");
        achilles::Settings settings(node);
        CHECK(settings.GetAs<int>("Root/Sub/Value") == 42);
    }
}
