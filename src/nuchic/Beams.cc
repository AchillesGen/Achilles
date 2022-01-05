#include "nuchic/Beams.hh"
#include "nuchic/Utilities.hh"

#include <fstream>

using nuchic::Spectrum;

Spectrum::Spectrum(const YAML::Node &node) {
    if(node["Histogram"]) {
        std::string filename = node["Histogram"].as<std::string>();
        std::ifstream hist(filename.c_str());
        std::string line;

        // Read in dummy header lines
        std::getline(hist, line);
        std::getline(hist, line);

        // Read in data
        std::vector<double> edges, heights;
        std::vector<std::string> tokens;
        while(std::getline(hist, line)) {
            tokens.clear();

            tokenize(line, tokens);
            edges.push_back(std::stod(tokens[0]));
            heights.push_back(std::stod(tokens[2]));
        }
        edges.push_back(std::stod(tokens[1]));

        // Create histogram
        m_flux = Histogram(edges, "flux");
        for(size_t i = 0; i < heights.size(); ++i) {
            m_flux.Fill((edges[i+1] - edges[i])/2, heights[i]);
        }
        m_min_energy = edges.front();
        m_max_energy = edges.back();
        m_delta_energy = m_max_energy - m_min_energy;
    } else {
        throw std::runtime_error("Spectrum: Only histogram fluxes are implemented");
    }
}

nuchic::FourVector Spectrum::Flux(const std::vector<double> &ran) const {
    double energy = ran[0]*m_delta_energy + m_min_energy;
    return {energy, 0, 0, energy};
}

double Spectrum::GenerateWeight(const FourVector &beam, std::vector<double> &ran) const {
    ran[0] = (beam.E() - m_min_energy) / m_delta_energy;
    return m_delta_energy;
}

nuchic::Beam::Beam(BeamMap beams) : m_beams{std::move(beams)} {
    n_vars = 0;
    for(const auto& beam : m_beams) {
        if(beam.second -> NVariables() > n_vars)
            n_vars = beam.second -> NVariables();
        if(m_pids.find(beam.first) != m_pids.end())
            throw std::logic_error(fmt::format("Multiple beams exist for PID: {}", int(beam.first)));
        spdlog::debug("Beam with PID: {} created.", beam.first);
        m_pids.insert(beam.first);
    }
}
