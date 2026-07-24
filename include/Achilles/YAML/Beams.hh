#ifndef YAML_BEAMS_HH
#define YAML_BEAMS_HH

#include "Achilles/Beams.hh"
#include "Achilles/ElectronPDF.hh"
#include "Achilles/FluxLoaders.hh"

#include <map>
#include <memory>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace YAML {

// A beam has a single flux Type shared by all of its species (kinds are never
// mixed). Spectrum draws from one or more files (each file self-describes its
// species); the analytic kinds list their species explicitly.
template <> struct convert<achilles::Beam> {
    static bool decode(const Node &node, achilles::Beam &rhs) {
        achilles::Beam::BeamMap beams;
        const auto type = node["Type"].as<std::string>();

        auto add_flux = [&beams](achilles::PID pid, std::shared_ptr<achilles::FluxType> flux) {
            if(beams.find(pid) != beams.end())
                throw std::logic_error(fmt::format("Multiple beams exist for PID: {}", int(pid)));
            beams[pid] = std::move(flux);
        };

        if(type == "Spectrum") {
            std::vector<std::map<achilles::PID, achilles::FluxData>> loaded;
            for(const auto &entry : node["Files"])
                loaded.push_back(achilles::LoaderForFile(entry)->Load(entry));

            bool multi_species_file = false;
            for(const auto &file : loaded)
                if(file.size() > 1) multi_species_file = true;
            if(multi_species_file && loaded.size() > 1)
                throw std::logic_error(
                    "Beam: a multi-species flux file must be the only file in the beam");

            for(const auto &file : loaded)
                for(const auto &flux : file)
                    add_flux(flux.first, std::make_shared<achilles::Spectrum>(flux.second));
        } else if(type == "Monochromatic") {
            auto flux = std::make_shared<achilles::Monochromatic>(node["Energy"].as<double>());
            for(const auto &s : node["Species"]) add_flux(achilles::PID(s.as<int>()), flux);
        } else if(type == "FlatFlux") {
            auto flux = std::make_shared<achilles::FlatFlux>(node);
            for(const auto &s : node["Species"]) add_flux(achilles::PID(s.as<int>()), flux);
        } else if(type == "Geometry") {
            // One shared GeometryBeam serves every species the driver may inject;
            // the ray (and species/target selection) are set per event at runtime.
            auto flux = std::make_shared<achilles::GeometryBeam>(node);
            for(const auto &s : node["Species"]) add_flux(achilles::PID(s.as<int>()), flux);
        } else {
            throw std::logic_error(fmt::format("Beam: unknown flux Type '{}'", type));
        }

        if(beams.empty()) throw std::logic_error("Beam: no species defined");

        rhs = achilles::Beam(std::move(beams));
        return true;
    }
};

} // namespace YAML

#endif
