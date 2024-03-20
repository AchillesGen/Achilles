#ifndef YAML_BEAMS_HH
#define YAML_BEAMS_HH

#include "Achilles/Beams.hh"
#include "Achilles/ElectronPDF.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace YAML {

template<>
struct convert<std::shared_ptr<achilles::FluxType>> {
    // FIXME: How to encode a generic flux
    static Node encode(const std::shared_ptr<achilles::FluxType> &rhs) {
        Node node;
        node["Type"] = rhs -> Type();

        return node;
    }

    static bool decode(const Node &node, std::shared_ptr<achilles::FluxType> &rhs) {
        // TODO: Improve checks to ensure the node is a valid beam (mainly validation)
        std::string type = node["Beam Params"]["Type"].as<std::string>();
        if(type == "Monochromatic") {
            auto energy = node["Beam Params"]["Energy"].as<double>();
            rhs = std::make_shared<achilles::Monochromatic>(energy);
            return true;
        } else if(type == "Spectrum") {
            rhs = std::make_shared<achilles::Spectrum>(node);
            return true;
        } else if(type == "PDFBeam") {
            rhs = std::make_shared<achilles::PDFBeam>(node);
        }

        return false;
    }
};

template<>
struct convert<achilles::Beam> {
    // TODO: Implement encoding
    static Node encode(const achilles::Beam &rhs) {
        Node node;

        for(const auto &beam : rhs.m_beams) {
            Node subnode;
            subnode["Beam"]["PID"] = static_cast<int>(beam.first);
            subnode["Beam"]["Beam Params"] = beam.second;
            node.push_back(subnode);
        }

        return node;
    } 

    static bool decode(const Node &node, achilles::Beam &rhs) {
        achilles::Beam::BeamMap beams; 

        for(YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
            YAML::Node beamNode = *it;
            auto pid = achilles::PID(beamNode["Beam"]["PID"].as<int>());
            auto beam = beamNode["Beam"].as<std::shared_ptr<achilles::FluxType>>();
            if(beams.find(pid) != beams.end())
                throw std::logic_error(fmt::format("Multiple beams exist for PID: {}", int(pid)));
            beams[pid] = beam;
        } 

        rhs = achilles::Beam(beams);
        return true;
    }
};

}

#endif
