#include "catch2/catch.hpp"

#include "Achilles/Beams.hh"
#include "Achilles/FluxLoaders.hh"
#include "Achilles/YAML/Beams.hh"
#include "yaml-cpp/yaml.h"

using achilles::PID;

TEST_CASE("Flux file loaders", "[Beams]") {
    SECTION("Achilles format is single-species and self-describing") {
        auto data = achilles::HistogramFluxLoader{}.Load(YAML::Load("flux/miniboone.dat"));
        REQUIRE(data.size() == 1);
        REQUIRE(data.count(PID(14)) == 1);
        CHECK(data.at(PID(14)).format == "Achilles");
        CHECK(data.at(PID(14)).min_energy == 0.0);
        CHECK(data.at(PID(14)).max_energy == 4450.0);
    }

    SECTION("MiniBooNE format yields all four species") {
        auto data = achilles::HistogramFluxLoader{}.Load(YAML::Load("flux/miniboone_nu.dat"));
        CHECK(data.size() == 4);
        CHECK(data.count(PID(14)) == 1);
        CHECK(data.count(PID(-14)) == 1);
        CHECK(data.count(PID(12)) == 1);
        CHECK(data.count(PID(-12)) == 1);
        CHECK(data.at(PID(14)).format == "MiniBooNE");
    }

    SECTION("T2K text format maps to numu") {
        auto data = achilles::HistogramFluxLoader{}.Load(YAML::Load("flux/T2K_nu.dat"));
        REQUIRE(data.size() == 1);
        CHECK(data.count(PID(14)) == 1);
        CHECK(data.at(PID(14)).format == "T2K");
        CHECK(data.at(PID(14)).max_energy == 30.0);
    }

    SECTION("CSV format maps columns to PIDs by name") {
        auto data = achilles::CsvFluxLoader{}.Load(YAML::Load("flux/t2k_flux.csv"));
        CHECK(data.size() == 4);
        CHECK(data.count(PID(14)) == 1);
        CHECK(data.count(PID(-14)) == 1);
        CHECK(data.count(PID(12)) == 1);
        CHECK(data.count(PID(-12)) == 1);
    }

    SECTION("LoaderForFile dispatches by file extension") {
        const auto entry = YAML::Load("flux/t2k_flux.csv");
        auto data = achilles::LoaderForFile(entry)->Load(entry);
        CHECK(data.size() == 4);
    }
}

TEST_CASE("SpeciesNameToPID", "[Beams]") {
    CHECK(achilles::SpeciesNameToPID("numu_flux") == PID(14));
    CHECK(achilles::SpeciesNameToPID("NuMuBar") == PID(-14));
    CHECK(achilles::SpeciesNameToPID("nue") == PID(12));
    CHECK(achilles::SpeciesNameToPID("nuebar_flux") == PID(-12));
    CHECK_THROWS(achilles::SpeciesNameToPID("proton"));
}

TEST_CASE("Spectrum from FluxData", "[Beams]") {
    auto data = achilles::HistogramFluxLoader{}.Load(YAML::Load("flux/dummy.dat"));
    achilles::Spectrum spectrum(data.at(PID(14)));

    CHECK(spectrum.NVariables() == 1);
    CHECK(spectrum.Format() == "Achilles");
    CHECK(spectrum.Flux({0.5}, 0) == achilles::FourVector(500, 0, 0, 500));
    std::vector<double> rans(1);
    CHECK(spectrum.GenerateWeight({500, 0, 0, 500}, rans, 0) == 1.0);
    CHECK(rans[0] == 0.5);
}

TEST_CASE("Beam from YAML", "[Beams]") {
    SECTION("Monochromatic with multiple species") {
        auto node = YAML::Load(R"beam(
Type: Monochromatic
Energy: 100
Species: [12, 13]
)beam");
        auto beam = node.as<achilles::Beam>();

        CHECK(beam.NBeams() == 2);
        CHECK(beam.BeamIDs() == std::set<PID>{PID(12), PID(13)});
        CHECK(beam.Domain().Min() == 100);
        CHECK(beam.Domain().Max() == 100);
        CHECK(beam.Flux(PID(12), {}, 0) == achilles::FourVector(100, 0, 0, 100));
        CHECK(beam.Flux(PID(13), {}, 0) == achilles::FourVector(100, 0, 0, 100));
    }

    SECTION("Spectrum aggregates single-species files over a shared support") {
        auto node = YAML::Load(R"beam(
Type: Spectrum
Files:
  - flux/minerva_numu_fhc.dat
  - flux/minerva_nue_fhc.dat
)beam");
        auto beam = node.as<achilles::Beam>();

        CHECK(beam.BeamIDs() == std::set<PID>{PID(14), PID(12)});
        CHECK(beam.Domain().Min() == 0.0);
        CHECK(beam.Domain().Max() == 100.0);
    }

    SECTION("Spectrum from a single multi-species file") {
        auto node = YAML::Load(R"beam(
Type: Spectrum
Files: [flux/miniboone_nu.dat]
)beam");
        auto beam = node.as<achilles::Beam>();

        CHECK(beam.BeamIDs() == std::set<PID>{PID(14), PID(-14), PID(12), PID(-12)});
    }

    SECTION("A multi-species file must be the only file in the beam") {
        auto node = YAML::Load(R"beam(
Type: Spectrum
Files:
  - flux/t2k_flux.csv
  - flux/minerva_numu_fhc.dat
)beam");
        CHECK_THROWS_WITH(node.as<achilles::Beam>(),
                          "Beam: a multi-species flux file must be the only file in the beam");
    }

    SECTION("Duplicate species across files are rejected") {
        auto node = YAML::Load(R"beam(
Type: Spectrum
Files:
  - flux/minerva_numu_fhc.dat
  - flux/minerva_numu_rhc.dat
)beam");
        CHECK_THROWS_WITH(node.as<achilles::Beam>(), "Multiple beams exist for PID: 14");
    }

    SECTION("Mismatched energy support is rejected") {
        auto node = YAML::Load(R"beam(
Type: Spectrum
Files:
  - flux/miniboone.dat
  - flux/minerva_nue_fhc.dat
)beam");
        CHECK_THROWS(node.as<achilles::Beam>());
    }
}
