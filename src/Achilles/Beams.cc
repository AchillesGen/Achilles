#include "Achilles/Beams.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Utilities.hh"

#include <fstream>

#ifdef USE_ROOT
#include "TFile.h"
#include "TH1D.h"
#endif

using achilles::PDFBeam;
using achilles::Spectrum;
using achilles::FlatFlux;

Spectrum::Spectrum(const YAML::Node &node) {
    spdlog::debug("Loading spectrum flux");
    if(node["Histogram"]) {
        std::string filename = node["Histogram"].as<std::string>();
        std::ifstream hist(filename.c_str());
        spdlog::trace("Beam::Spectrum: Loading data from: {}", filename);
        if(hist.bad()) {
            std::string msg = fmt::format("Beam::Spectrum: Can't open file {}", filename);
            throw std::runtime_error(msg);
        }
        std::string line;

        // Read in format line
        std::getline(hist, line);
        size_t elo_token, ehi_token, flux_token;
        if(line.find("Achilles") != std::string::npos) {
            m_format = FluxFormat::Achilles;
            elo_token = 0;
            ehi_token = 1;
            flux_token = 2;
            AchillesHeader(hist);
        } else if(line.find("MiniBooNE") != std::string::npos) {
            m_energy_units = 1.0 / 1.0_GeV;
            m_format = FluxFormat::MiniBooNE;
            elo_token = 0;
            ehi_token = 1;
            if(line.find("nu mode") != std::string::npos) {
                flux_token = 2;
            } else {
                flux_token = 3;
            }
            MiniBooNEHeader(hist);
        } else if(line.find("ND280") != std::string::npos) {
            m_energy_units = 1.0 / 1.0_GeV;
            m_format = FluxFormat::T2K;
            T2KHeader(hist);
            elo_token = 1;
            ehi_token = 2;
            flux_token = 3;
        } else {
            std::string msg =
                fmt::format("Beam::Spectrum: Invalid flux format from file {}", filename);
            throw std::runtime_error(msg);
        }

        // Read in data
        std::vector<std::string> tokens;
        std::vector<double> edges, heights;
        while(std::getline(hist, line)) {
            tokens.clear();

            tokenize(line, tokens);
            edges.push_back(std::stod(tokens[elo_token]));
            heights.push_back(std::stod(tokens[flux_token]));
        }
        edges.push_back(std::stod(tokens[ehi_token]));

        // Calculate Integral
        bool use_width{};
        switch(m_units) {
        case flux_units::v_m2_POT_500MeV:
            m_energy_units = 1.0 / 1.0_GeV;
            use_width = true;
            break;
        case flux_units::v_cm2_POT_50MeV:
            use_width = true;
            break;
        case flux_units::cm2_50MeV:
            use_width = true;
            break;
        case flux_units::v_cm2_POT_MeV:
        case flux_units::v_nb_POT_MeV:
            use_width = false;
            break;
        }

        for(size_t i = 0; i < heights.size(); ++i) {
            double width = use_width ? edges[i + 1] - edges[i] : 1;
            m_flux_integral += width * heights[i];
        }
        spdlog::trace("Flux integral = {}", m_flux_integral);

        // Create interpolation function
        std::vector<double> bin_centers;
        bin_centers.push_back(edges[0]);
        heights.insert(heights.begin(), heights[0]);
        for(size_t i = 1; i < edges.size(); ++i)
            bin_centers.push_back((edges[i] + edges[i - 1]) / 2);
        bin_centers.push_back(edges.back());
        heights.push_back(heights.back());

        Interp1D interp(bin_centers, heights, InterpolationType::Polynomial);
        interp.SetPolyOrder(1);

        m_flux = [=](double x) { return interp(x); };
        m_min_energy = edges.front();
        m_max_energy = edges.back();
        m_delta_energy = m_max_energy - m_min_energy;
    } else if(node["HepData"]) {
        spdlog::debug("Reading from HepData yaml file");
        std::string filename = node["HepData"].as<std::string>();
        YAML::Node root = YAML::LoadFile(filename);
        std::vector<double> heights;
        // TODO: Handle multiple fluxes??
        for(const auto &value : root["dependent_variables"][0]["values"]) {
            heights.push_back(value["value"].as<double>());
        }

        auto energy_units = root["independent_variables"][0]["header"]["units"].as<std::string>();
        if(energy_units == "GeV") { m_energy_units = 1.0 / 1_GeV; }
        size_t nbins = root["independent_variables"][0]["values"].size();
        std::vector<double> bin_centers, low, high;
        m_min_energy = root["independent_variables"][0]["values"][0]["low"].as<double>();
        m_max_energy = root["independent_variables"][0]["values"][nbins - 1]["high"].as<double>();
        for(const auto &bin : root["independent_variables"][0]["values"]) {
            low.push_back(bin["low"].as<double>());
            high.push_back(bin["high"].as<double>());
            bin_centers.push_back((low.back() + high.back()) / 2);
        }

        // TODO: Implement check on units to use bin width or not
        bool use_width = true;
        for(size_t i = 0; i < heights.size(); ++i) {
            double width = use_width ? high[i] - low[i] : 1;
            m_flux_integral += width * heights[i];
        }
        spdlog::trace("Flux integral = {}", m_flux_integral);

        // Create interpolation function
        bin_centers.insert(bin_centers.begin(), m_min_energy);
        heights.insert(heights.begin(), heights[0]);
        bin_centers.push_back(m_max_energy);
        heights.push_back(heights.back());

        Interp1D interp(bin_centers, heights, InterpolationType::Polynomial);
        interp.SetPolyOrder(1);

        m_flux = [=](double x) { return interp(x); };
        m_delta_energy = m_max_energy - m_min_energy;
    } else if(node["ROOTHist"]) {
#ifdef USE_ROOT
        std::string filename = "flux/" + node["ROOTHist"]["File"].as<std::string>();
        spdlog::trace("Reading flux file: {}", filename);
        TFile *file = new TFile(filename.c_str());
        TH1D *hist =
            static_cast<TH1D *>(file->Get(node["ROOTHist"]["HistName"].as<std::string>().c_str()));
        bool use_width = node["ROOTHist"]["UseWidth"].as<bool>();
        double norm = node["ROOTHist"]["Norm"].as<double>();
        hist->Scale(norm);
        std::vector<double> bin_centers; //(static_cast<size_t>(hist -> GetNbinsX())+2);
        std::vector<double> heights;     //(static_cast<size_t>(hist -> GetNbinsX())+2);
        double width = use_width ? hist->GetBinWidth(1) : 1;
        bin_centers.push_back(hist->GetBinLowEdge(1));
        heights.push_back(hist->GetBinContent(1) / width);
        size_t i;
        for(i = 1; i <= static_cast<size_t>(hist->GetNbinsX()); ++i) {
            if(use_width) width = hist->GetBinWidth(static_cast<int>(i));
            double height = hist->GetBinContent(static_cast<int>(i)) / width;
            if(height == 0) break;
            bin_centers.push_back(hist->GetBinCenter(static_cast<int>(i)));
            heights.push_back(height);
        }
        bin_centers.push_back(hist->GetBinLowEdge(static_cast<int>(i)) +
                              hist->GetBinWidth(static_cast<int>(i)));
        heights.push_back(heights.back());

        Interp1D interp(bin_centers, heights, InterpolationType::Polynomial);
        interp.SetPolyOrder(1);
        m_flux_integral = interp.Integrate();
        spdlog::trace("Flux Integral = {}", m_flux_integral);

        m_flux = [=](double x) { return interp(x); };
        m_min_energy = bin_centers.front();
        m_max_energy = bin_centers.back();
        m_delta_energy = m_max_energy - m_min_energy;
        m_energy_units = 1.0 / 1.0_GeV;
        spdlog::debug("Flux energy range: [{}, {}]", m_min_energy, m_max_energy);
        spdlog::trace("Flux histogram:");
        for(size_t j = 0; j < bin_centers.size(); ++j) {
            spdlog::trace("  {} -> {}", bin_centers[j], heights[j]);
        }
#else
        throw std::runtime_error("Achilles has not been compiled with ROOT support");
#endif
    } else {
        throw std::runtime_error("Spectrum: Only histogram fluxes are implemented");
    }
}

void Spectrum::AchillesHeader(std::ifstream &hist) {
    // Parse experiment description
    std::string line;
    std::getline(hist, line);

    // Parse units
    std::getline(hist, line);
    std::vector<std::string> tokens;
    tokenize(line, tokens);
    spdlog::debug("{}, {}", tokens[0], tokens[1]);
    if(tokens[0] != "units:") throw std::runtime_error("Beam::Spectrum: Invalid file format");
    if(tokens[1] == "v/cm^2/POT/MeV") {
        m_units = flux_units::v_cm2_POT_MeV;
    } else if(tokens[1] == "v/cm^2/POT/50MeV") {
        m_units = flux_units::v_cm2_POT_50MeV;
    } else if(tokens[1] == "cm^{-2}/50MeV") {
        m_units = flux_units::cm2_50MeV;
    } else if(tokens[1] == "v/m^2/POT/500MeV") {
        m_units = flux_units::v_m2_POT_500MeV;
    } else {
        throw std::runtime_error("Beam::Spectrum: Invalid flux units");
    }

    // Read in histogram header line
    std::getline(hist, line);
}

void Spectrum::MiniBooNEHeader(std::ifstream &hist) {
    // Parse header
    std::string line;
    for(size_t i = 0; i < 3; ++i) std::getline(hist, line);

    // Get units
    std::getline(hist, line);
    if(line.find("cm^2/proton-on-target/50 MeV") != std::string::npos) {
        m_units = flux_units::v_cm2_POT_50MeV;
    } else {
        throw std::runtime_error("Beam::Spectrum: Invalid flux units");
    }

    // Read remainder of header
    for(size_t i = 0; i < 4; ++i) std::getline(hist, line);
}

void Spectrum::T2KHeader(std::ifstream &hist) {
    // Get units
    std::string line;
    std::getline(hist, line);
    if(line.find("cm^{-2}/50MeV") != std::string::npos) {
        m_units = flux_units::cm2_50MeV;
    } else {
        throw std::runtime_error("Beam::Spectrum: Invalid flux units");
    }
}

std::string Spectrum::Format() const {
    switch(m_format) {
    case FluxFormat::Achilles:
        return "Achilles";
    case FluxFormat::MiniBooNE:
        return "MiniBooNE";
    case FluxFormat::T2K:
        return "T2K";
    }
    return "Undefined";
}

achilles::FourVector Spectrum::Flux(const std::vector<double> &ran, double min_energy) const {
    min_energy = std::max(min_energy * m_energy_units, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    double energy = (ran[0] * delta_energy + min_energy) / m_energy_units;
    return {energy, 0, 0, energy};
}

double Spectrum::GenerateWeight(const FourVector &beam, std::vector<double> &ran,
                                double min_energy) const {
    min_energy = std::max(min_energy * m_energy_units, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    ran[0] = (beam.E() * m_energy_units - min_energy) / delta_energy;
    return (delta_energy * m_flux(beam.E() * m_energy_units)) / m_flux_integral;
}

double Spectrum::EvaluateFlux(const FourVector &beam) const {
    return m_flux(beam.E() * m_energy_units);
}

PDFBeam::PDFBeam(const YAML::Node &) {
    throw std::runtime_error("PDFBeam: Not implemented yet!");
}

achilles::FourVector PDFBeam::Flux(const std::vector<double> &, double) const {
    throw std::runtime_error("PDFBeam: Not implemented yet!");
    return {};
}

double PDFBeam::GenerateWeight(const FourVector &, std::vector<double> &, double) const {
    throw std::runtime_error("PDFBeam: Not implemented yet!");
    return {};
}

double PDFBeam::EvaluateFlux(const FourVector &) const {
    throw std::runtime_error("PDFBeam: Not implemented yet!");
    return {};
}

FlatFlux::FlatFlux(const YAML::Node &node) {
    m_min_energy = node["MinEnergy"].as<double>();
    m_max_energy = node["MaxEnergy"].as<double>();
}

achilles::FourVector FlatFlux::Flux(const std::vector<double> &ran, double min_energy) const {
    min_energy = std::max(min_energy, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    double energy = ran[0]*delta_energy + min_energy;
    return {energy, 0, 0, energy};
}

double FlatFlux::GenerateWeight(const FourVector &beam, std::vector<double> &ran,
                                double min_energy) const {
    min_energy = std::max(min_energy, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    ran[0] = (beam.E() - min_energy) / delta_energy;
    return delta_energy;
}

achilles::Beam::Beam(BeamMap beams) : m_beams{std::move(beams)} {
    n_vars = 0;
    for(const auto &beam : m_beams) {
        if(beam.second->NVariables() > n_vars) n_vars = beam.second->NVariables();
        if(m_pids.find(beam.first) != m_pids.end())
            throw std::logic_error(
                fmt::format("Multiple beams exist for PID: {}", int(beam.first)));
        spdlog::debug("Beam with PID: {} created.", beam.first);
        m_pids.insert(beam.first);
        if(m_beams.begin()->second->MinEnergy() != beam.second->MinEnergy())
            throw std::logic_error(
                fmt::format("Beams must have the same minimum energies. Got {} and {}",
                            m_beams.begin()->second->MinEnergy(), beam.second->MinEnergy()));
        if(m_beams.begin()->second->MaxEnergy() != beam.second->MaxEnergy())
            throw std::logic_error(
                fmt::format("Beams must have the same minimum energies. Got {} and {}",
                            m_beams.begin()->second->MaxEnergy(), beam.second->MaxEnergy()));
    }
}

double achilles::Beam::MaxEnergy() const {
    double max = 0;
    for(const auto &beam : m_beams) {
        if(beam.second->MaxEnergy() > max) max = beam.second->MaxEnergy();
    }
    return max;
}
