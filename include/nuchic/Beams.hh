#ifndef BEAMS_HH
#define BEAMS_HH

#include "nuchic/Histogram.hh"
#include <set>
#include <memory>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "nuchic/FourVector.hh"
#include "nuchic/ParticleInfo.hh"

namespace nuchic {

class FluxType {
    public:
        FluxType() = default;
        FluxType(const FluxType&) = default;
        FluxType(FluxType&&) = default;
        FluxType& operator=(const FluxType&) = default;
        FluxType& operator=(FluxType&&) = default;
        virtual ~FluxType() = default;

        virtual int NVariables() const = 0;
        virtual FourVector Flux(const std::vector<double>&) const = 0;
        virtual double GenerateWeight(const FourVector&, std::vector<double>&) const = 0;
};

class Monochromatic : public FluxType {
    public:
        Monochromatic(const double &energy) : m_energy(energy) {}

        int NVariables() const override {
            return 0;
        }

        FourVector Flux(const std::vector<double>&) const override {
            return {m_energy, 0, 0, m_energy};
        }

        double GenerateWeight(const FourVector&, std::vector<double>&) const override {
            return 1;
        }

    private:
        double m_energy;
};


// TODO: Figure out how to handle spectrums
// How to handle this correctly?
// 1. Generate the energy according to some distribution
// 2. A. Multiply maximum energy by a fraction to get neutrino energy
//    B. Calculate the probability for this energy
// The first requires being able to appropriately sample from the distribution
// (i.e. having a form for the inverse of the CDF). The second would take advantage
// of the importance sampling of Vegas. Furthermore, the second would make it easier to
// combine beams of different initial state particles in a straightforward manner
class Spectrum : public FluxType {
    public:
        enum class Type {
            Histogram,
        };

        Spectrum(const YAML::Node&);

        int NVariables() const override { return 1; }

        FourVector Flux(const std::vector<double>&) const override;

        double GenerateWeight(const FourVector&, std::vector<double>&) const override;

    private:
        Histogram m_flux{};
        double m_min_energy{}, m_max_energy{};
        double m_delta_energy{};
};

class Beam {
    public:
        using BeamMap = std::map<nuchic::PID, std::shared_ptr<FluxType>>;

        Beam(BeamMap beams);
        Beam(const Beam&) = delete;
        Beam(Beam&&) = default;
        Beam& operator=(const Beam&) = delete;
        Beam& operator=(Beam&&) = default;
        virtual ~Beam() = default;

        Beam() { n_vars = 0; }
        int NVariables() const { return  n_vars; }
        virtual FourVector Flux(const PID pid, const std::vector<double> &rans) const { 
            return m_beams.at(pid) -> Flux(rans); 
        }
        double GenerateWeight(const PID pid, const FourVector &p, std::vector<double> &rans) const { 
            return m_beams.at(pid) -> GenerateWeight(p, rans); 
        }
        size_t NBeams() const { return m_beams.size(); }
        const std::set<PID>& BeamIDs() const { return m_pids; }

        // Accessors
        std::shared_ptr<FluxType>& operator[](const PID pid) { return m_beams[pid]; }
        std::shared_ptr<FluxType> at(const PID pid) const { return m_beams.at(pid); }
        std::shared_ptr<FluxType> operator[](const PID pid) const { return m_beams.at(pid); }

        friend YAML::convert<Beam>;

    private:
        int n_vars;
        std::set<PID> m_pids;
        BeamMap m_beams;
};

}

namespace YAML {

template<>
struct convert<std::shared_ptr<nuchic::FluxType>> {
    // FIXME: How to encode a generic flux
    static Node encode(const std::shared_ptr<nuchic::FluxType> &rhs) {
        Node node;
        node["Type"] = "Monochromatic";
        node["Energy"] = rhs->Flux({}).E();

        return node;
    }

    static bool decode(const Node &node, std::shared_ptr<nuchic::FluxType> &rhs) {
        // TODO: Improve checks to ensure the node is a valid beam (mainly validation)
        if(node["Type"].as<std::string>() == "Monochromatic") {
            auto energy = node["Energy"].as<double>();
            rhs = std::make_shared<nuchic::Monochromatic>(energy);
            return true;
        } else if(node["Type"].as<std::string>() == "Spectrum") {
            rhs = std::make_shared<nuchic::Spectrum>(node);
            return true;
        }

        return false;
    }
};

template<>
struct convert<nuchic::Beam> {
    // TODO: Implement encoding
    static Node encode(const nuchic::Beam &rhs) {
        Node node;

        for(const auto &beam : rhs.m_beams) {
            Node subnode;
            subnode["Beam"]["PID"] = static_cast<int>(beam.first);
            subnode["Beam"]["Beam Params"] = beam.second;
            node.push_back(subnode);
        }

        return node;
    } 

    static bool decode(const Node &node, nuchic::Beam &rhs) {
        nuchic::Beam::BeamMap beams; 

        for(YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
            YAML::Node beamNode = *it;
            auto pid = nuchic::PID(beamNode["Beam"]["PID"].as<int>());
            auto beam = beamNode["Beam"]["Beam Params"].as<std::shared_ptr<nuchic::FluxType>>();
            if(beams.find(pid) != beams.end())
                throw std::logic_error(fmt::format("Multiple beams exist for PID: {}", int(pid)));
            beams[pid] = beam;
        } 

        rhs = nuchic::Beam(beams);
        return true;
    }
};

}

#endif
