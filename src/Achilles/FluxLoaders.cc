#include "Achilles/FluxLoaders.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Utilities.hh"

#include "spdlog/spdlog.h"

#include <fstream>
#include <stdexcept>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#ifdef USE_ROOT
#include "TFile.h"
#include "TH1D.h"
#endif

namespace achilles {

namespace {

// Flux normalization units, used to decide whether the integral weights by bin
// width. Detected while parsing a plain-text histogram header.
enum class flux_units {
    v_m2_POT_500MeV,
    v_nb_POT_MeV,
    v_cm2_POT_MeV,
    v_cm2_POT_50MeV,
    v_cm2_POT_125MeV,
    v_cm2_POT_100MeV,
    cm2_50MeV,
};

flux_units AchillesHeader(std::ifstream &hist) {
    // Parse experiment description
    std::string line;
    std::getline(hist, line);

    // Parse units
    std::getline(hist, line);
    std::vector<std::string> tokens;
    tokenize(line, tokens);
    spdlog::debug("{}, {}", tokens[0], tokens[1]);
    if(tokens[0] != "units:") throw std::runtime_error("Beam::Spectrum: Invalid file format");
    flux_units units;
    if(tokens[1] == "v/cm^2/POT/MeV") {
        units = flux_units::v_cm2_POT_MeV;
    } else if(tokens[1] == "v/cm^2/POT/50MeV") {
        units = flux_units::v_cm2_POT_50MeV;
    } else if(tokens[1] == "cm^{-2}/50MeV") {
        units = flux_units::cm2_50MeV;
    } else if(tokens[1] == "v/m^2/POT/500MeV") {
        units = flux_units::v_m2_POT_500MeV;
    } else if(tokens[1] == "v/cm^2/POT/125MeV") {
        units = flux_units::v_cm2_POT_125MeV;
    } else if(tokens[1] == "v/cm^2/POT/100MeV") {
        units = flux_units::v_cm2_POT_100MeV;
    } else {
        throw std::runtime_error("Beam::Spectrum: Invalid flux units");
    }

    // Read in histogram header line
    std::getline(hist, line);
    return units;
}

flux_units MiniBooNEHeader(std::ifstream &hist) {
    // Parse header
    std::string line;
    for(size_t i = 0; i < 3; ++i) std::getline(hist, line);

    // Get units
    std::getline(hist, line);
    flux_units units;
    if(line.find("cm^2/proton-on-target/50 MeV") != std::string::npos) {
        units = flux_units::v_cm2_POT_50MeV;
    } else {
        throw std::runtime_error("Beam::Spectrum: Invalid flux units");
    }

    // Read remainder of header
    for(size_t i = 0; i < 4; ++i) std::getline(hist, line);
    return units;
}

flux_units T2KHeader(std::ifstream &hist) {
    // Get units
    std::string line;
    std::getline(hist, line);
    flux_units units;
    if(line.find("cm^{-2}/50MeV") != std::string::npos) {
        units = flux_units::cm2_50MeV;
    } else {
        throw std::runtime_error("Beam::Spectrum: Invalid flux units");
    }
    return units;
}

} // namespace

FluxData HistogramFluxLoader::Load(const YAML::Node &node) const {
    FluxData data;
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
    flux_units units;
    if(line.find("Achilles") != std::string::npos) {
        data.format = "Achilles";
        elo_token = 0;
        ehi_token = 1;
        flux_token = 2;
        units = AchillesHeader(hist);
    } else if(line.find("MiniBooNE") != std::string::npos) {
        data.energy_units = 1.0 / 1.0_GeV;
        data.format = "MiniBooNE";
        elo_token = 0;
        ehi_token = 1;
        if(line.find("nu mode") != std::string::npos) {
            flux_token = 2;
        } else {
            flux_token = 3;
        }
        units = MiniBooNEHeader(hist);
    } else if(line.find("ND280") != std::string::npos) {
        data.energy_units = 1.0 / 1.0_GeV;
        data.format = "T2K";
        units = T2KHeader(hist);
        elo_token = 1;
        ehi_token = 2;
        flux_token = 3;
    } else {
        std::string msg = fmt::format("Beam::Spectrum: Invalid flux format from file {}", filename);
        throw std::runtime_error(msg);
    }

    // Read in data
    std::vector<std::string> tokens;
    std::vector<double> edges, heights;
    while(std::getline(hist, line)) {
        spdlog::trace("Flux line: {}", line);
        tokens.clear();

        tokenize(line, tokens);
        spdlog::trace("elo_token = {}", tokens[elo_token]);
        edges.push_back(std::stod(tokens[elo_token]));
        spdlog::trace("Parsed edge");
        heights.push_back(std::stod(tokens[flux_token]));
        spdlog::trace("Parsed flux");
    }
    edges.push_back(std::stod(tokens[ehi_token]));

    // Calculate Integral
    bool use_width{};
    switch(units) {
    case flux_units::v_m2_POT_500MeV:
        data.energy_units = 1.0 / 1.0_GeV;
        use_width = true;
        break;
    case flux_units::v_cm2_POT_100MeV:
    case flux_units::v_cm2_POT_125MeV:
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
        data.integral += width * heights[i];
    }
    spdlog::trace("Flux integral = {}", data.integral);

    // Create interpolation points
    std::vector<double> bin_centers;
    bin_centers.push_back(edges[0]);
    heights.insert(heights.begin(), heights[0]);
    for(size_t i = 1; i < edges.size(); ++i) bin_centers.push_back((edges[i] + edges[i - 1]) / 2);
    bin_centers.push_back(edges.back());
    heights.push_back(heights.back());

    data.bin_centers = std::move(bin_centers);
    data.heights = std::move(heights);
    data.min_energy = edges.front();
    data.max_energy = edges.back();
    return data;
}

FluxData HepDataFluxLoader::Load(const YAML::Node &node) const {
    FluxData data;
    spdlog::debug("Reading from HepData yaml file");
    std::string filename = node["HepData"].as<std::string>();
    YAML::Node root = YAML::LoadFile(filename);
    std::vector<double> heights;
    // TODO: Handle multiple fluxes??
    for(const auto &value : root["dependent_variables"][0]["values"]) {
        heights.push_back(value["value"].as<double>());
    }

    auto energy_units = root["independent_variables"][0]["header"]["units"].as<std::string>();
    if(energy_units == "GeV") { data.energy_units = 1.0 / 1_GeV; }
    size_t nbins = root["independent_variables"][0]["values"].size();
    std::vector<double> bin_centers, low, high;
    data.min_energy = root["independent_variables"][0]["values"][0]["low"].as<double>();
    data.max_energy = root["independent_variables"][0]["values"][nbins - 1]["high"].as<double>();
    for(const auto &bin : root["independent_variables"][0]["values"]) {
        low.push_back(bin["low"].as<double>());
        high.push_back(bin["high"].as<double>());
        bin_centers.push_back((low.back() + high.back()) / 2);
    }

    // TODO: Implement check on units to use bin width or not
    bool use_width = true;
    for(size_t i = 0; i < heights.size(); ++i) {
        double width = use_width ? high[i] - low[i] : 1;
        data.integral += width * heights[i];
    }
    spdlog::trace("Flux integral = {}", data.integral);

    // Create interpolation points
    bin_centers.insert(bin_centers.begin(), data.min_energy);
    heights.insert(heights.begin(), heights[0]);
    bin_centers.push_back(data.max_energy);
    heights.push_back(heights.back());

    data.bin_centers = std::move(bin_centers);
    data.heights = std::move(heights);
    return data;
}

FluxData RootHistFluxLoader::Load(const YAML::Node &node) const {
#ifdef USE_ROOT
    FluxData data;
    std::string filename = "flux/" + node["ROOTHist"]["File"].as<std::string>();
    spdlog::trace("Reading flux file: {}", filename);
    TFile *file = new TFile(filename.c_str());
    TH1D *hist =
        static_cast<TH1D *>(file->Get(node["ROOTHist"]["HistName"].as<std::string>().c_str()));
    bool use_width = node["ROOTHist"]["UseWidth"].as<bool>();
    double norm = node["ROOTHist"]["Norm"].as<double>();
    hist->Scale(norm);
    std::vector<double> bin_centers;
    std::vector<double> heights;
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
    data.integral = interp.Integrate();
    spdlog::trace("Flux Integral = {}", data.integral);

    data.bin_centers = std::move(bin_centers);
    data.heights = std::move(heights);
    data.min_energy = data.bin_centers.front();
    data.max_energy = data.bin_centers.back();
    data.energy_units = 1.0 / 1.0_GeV;
    spdlog::debug("Flux energy range: [{}, {}]", data.min_energy, data.max_energy);
    return data;
#else
    (void)node;
    throw std::runtime_error("Achilles has not been compiled with ROOT support");
#endif
}

} // namespace achilles
