#include "Achilles/FluxLoaders.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Utilities.hh"

#include "spdlog/spdlog.h"

#include <algorithm>
#include <fstream>
#include <sstream>
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

bool UseWidth(flux_units units) {
    switch(units) {
    case flux_units::v_cm2_POT_MeV:
    case flux_units::v_nb_POT_MeV:
        return false;
    case flux_units::v_m2_POT_500MeV:
    case flux_units::v_cm2_POT_100MeV:
    case flux_units::v_cm2_POT_125MeV:
    case flux_units::v_cm2_POT_50MeV:
    case flux_units::cm2_50MeV:
        return true;
    }
    return false;
}

flux_units ParseUnits(const std::string &units) {
    if(units == "v/cm^2/POT/MeV") return flux_units::v_cm2_POT_MeV;
    if(units == "v/cm^2/POT/50MeV") return flux_units::v_cm2_POT_50MeV;
    if(units == "cm^{-2}/50MeV") return flux_units::cm2_50MeV;
    if(units == "v/m^2/POT/500MeV") return flux_units::v_m2_POT_500MeV;
    if(units == "v/cm^2/POT/125MeV") return flux_units::v_cm2_POT_125MeV;
    if(units == "v/cm^2/POT/100MeV") return flux_units::v_cm2_POT_100MeV;
    throw std::runtime_error("Beam::Spectrum: Invalid flux units");
}

// Achilles format header: description, "units: ...", "species: <pids...>", and a
// column-header line. Returns the units and fills the per-column species list.
flux_units AchillesHeader(std::ifstream &hist, std::vector<PID> &species) {
    std::string line;
    std::getline(hist, line); // description

    std::getline(hist, line); // units
    std::vector<std::string> tokens;
    tokenize(line, tokens);
    if(tokens.size() < 2 || tokens[0] != "units:")
        throw std::runtime_error("Beam::Spectrum: Invalid file format");
    flux_units units = ParseUnits(tokens[1]);

    std::getline(hist, line); // species
    tokens.clear();
    tokenize(line, tokens);
    if(tokens.empty() || tokens[0] != "species:")
        throw std::runtime_error("Beam::Spectrum: Achilles flux missing 'species:' header");
    for(size_t i = 1; i < tokens.size(); ++i) species.emplace_back(std::stoi(tokens[i]));
    if(species.empty())
        throw std::runtime_error("Beam::Spectrum: Achilles flux declares no species");

    std::getline(hist, line); // column header line
    return units;
}

flux_units MiniBooNEHeader(std::ifstream &hist) {
    std::string line;
    for(size_t i = 0; i < 3; ++i) std::getline(hist, line);
    std::getline(hist, line);
    flux_units units;
    if(line.find("cm^2/proton-on-target/50 MeV") != std::string::npos)
        units = flux_units::v_cm2_POT_50MeV;
    else
        throw std::runtime_error("Beam::Spectrum: Invalid flux units");
    for(size_t i = 0; i < 4; ++i) std::getline(hist, line);
    return units;
}

flux_units T2KHeader(std::ifstream &hist) {
    std::string line;
    std::getline(hist, line);
    if(line.find("cm^{-2}/50MeV") == std::string::npos)
        throw std::runtime_error("Beam::Spectrum: Invalid flux units");
    return flux_units::cm2_50MeV;
}

// Turn a set of bin edges and one species' bin heights into padded interpolation
// points plus the flux integral.
FluxData BuildFluxData(const std::vector<double> &edges, std::vector<double> heights,
                       bool use_width, double energy_units, const std::string &format) {
    FluxData data;
    data.energy_units = energy_units;
    data.format = format;

    for(size_t i = 0; i < heights.size(); ++i) {
        double width = use_width ? edges[i + 1] - edges[i] : 1;
        data.integral += width * heights[i];
    }

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

std::string ToLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

std::vector<std::string> SplitCSV(const std::string &line) {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string item;
    while(std::getline(ss, item, ',')) tokens.push_back(item);
    return tokens;
}

} // namespace

PID SpeciesNameToPID(const std::string &name) {
    std::string key = ToLower(name);
    const std::string suffix = "_flux";
    if(key.size() > suffix.size() &&
       key.compare(key.size() - suffix.size(), suffix.size(), suffix) == 0)
        key = key.substr(0, key.size() - suffix.size());

    static const std::map<std::string, int> table = {
        {"numu", 14},    {"numubar", -14}, {"nue", 12},
        {"nuebar", -12}, {"nutau", 16},    {"nutaubar", -16},
    };
    auto it = table.find(key);
    if(it == table.end())
        throw std::runtime_error(fmt::format("Beam: Unknown flux species name '{}'", name));
    return PID(it->second);
}

std::string FluxEntryPath(const YAML::Node &entry) {
    if(entry.IsScalar()) return entry.as<std::string>();
    return entry["Path"].as<std::string>();
}

std::unique_ptr<FluxFileLoader> LoaderForFile(const YAML::Node &entry) {
    const std::string path = FluxEntryPath(entry);
    const auto dot = path.find_last_of('.');
    const std::string ext = dot == std::string::npos ? "" : ToLower(path.substr(dot + 1));

    std::string name;
    if(ext == "dat")
        name = "Histogram";
    else if(ext == "csv")
        name = "Csv";
    else if(ext == "root")
        name = "ROOTHist";
    else if(ext == "yaml" || ext == "yml")
        name = "HepData";
    else
        throw std::runtime_error(
            fmt::format("Beam: No flux loader for file extension '{}' ({})", ext, path));

    return Factory<FluxFileLoader>::Initialize(name);
}

std::map<PID, FluxData> HistogramFluxLoader::Load(const YAML::Node &entry) const {
    const std::string filename = FluxEntryPath(entry);
    std::ifstream hist(filename.c_str());
    spdlog::trace("Beam::Spectrum: Loading data from: {}", filename);
    if(hist.bad())
        throw std::runtime_error(fmt::format("Beam::Spectrum: Can't open file {}", filename));

    std::string line;
    std::getline(hist, line); // format line

    std::vector<PID> species;
    size_t elo_token, ehi_token;
    std::vector<size_t> flux_tokens;
    flux_units units;
    double energy_units = 1;
    std::string format;
    if(line.find("Achilles") != std::string::npos) {
        format = "Achilles";
        elo_token = 0;
        ehi_token = 1;
        units = AchillesHeader(hist, species);
        for(size_t i = 0; i < species.size(); ++i) flux_tokens.push_back(2 + i);
    } else if(line.find("MiniBooNE") != std::string::npos) {
        energy_units = 1.0 / 1.0_GeV;
        format = "MiniBooNE";
        elo_token = 0;
        ehi_token = 1;
        species = {PID(14), PID(-14), PID(12), PID(-12)};
        flux_tokens = {2, 3, 4, 5};
        units = MiniBooNEHeader(hist);
    } else if(line.find("ND280") != std::string::npos) {
        energy_units = 1.0 / 1.0_GeV;
        format = "T2K";
        units = T2KHeader(hist);
        elo_token = 1;
        ehi_token = 2;
        species = {PID(14)};
        flux_tokens = {3};
    } else {
        throw std::runtime_error(
            fmt::format("Beam::Spectrum: Invalid flux format from file {}", filename));
    }

    std::vector<std::string> tokens;
    std::vector<double> edges;
    std::vector<std::vector<double>> heights(species.size());
    while(std::getline(hist, line)) {
        if(line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        tokens.clear();
        tokenize(line, tokens);
        edges.push_back(std::stod(tokens[elo_token]));
        for(size_t j = 0; j < species.size(); ++j)
            heights[j].push_back(std::stod(tokens[flux_tokens[j]]));
    }
    edges.push_back(std::stod(tokens[ehi_token]));

    if(units == flux_units::v_m2_POT_500MeV) energy_units = 1.0 / 1.0_GeV;
    const bool use_width = UseWidth(units);

    std::map<PID, FluxData> result;
    for(size_t j = 0; j < species.size(); ++j)
        result[species[j]] = BuildFluxData(edges, heights[j], use_width, energy_units, format);
    return result;
}

std::map<PID, FluxData> CsvFluxLoader::Load(const YAML::Node &entry) const {
    const std::string filename = FluxEntryPath(entry);
    std::ifstream file(filename.c_str());
    spdlog::trace("Beam::Spectrum: Loading CSV data from: {}", filename);
    if(file.bad())
        throw std::runtime_error(fmt::format("Beam::Spectrum: Can't open file {}", filename));

    std::string line;
    std::getline(file, line); // header row
    const auto headers = SplitCSV(line);

    size_t elo_token = 0, ehi_token = 0;
    bool have_elo = false, have_ehi = false;
    std::vector<std::pair<size_t, PID>> flux_cols;
    for(size_t i = 0; i < headers.size(); ++i) {
        const std::string name = ToLower(headers[i]);
        if(name == "low") {
            elo_token = i;
            have_elo = true;
        } else if(name == "high") {
            ehi_token = i;
            have_ehi = true;
        } else if(name == "bin") {
            continue;
        } else {
            flux_cols.emplace_back(i, SpeciesNameToPID(headers[i]));
        }
    }
    if(!have_elo || !have_ehi)
        throw std::runtime_error("Beam::Spectrum: CSV flux missing 'low'/'high' columns");

    std::vector<double> edges;
    std::vector<std::vector<double>> heights(flux_cols.size());
    std::vector<std::string> tokens;
    while(std::getline(file, line)) {
        if(line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        tokens = SplitCSV(line);
        edges.push_back(std::stod(tokens[elo_token]));
        for(size_t j = 0; j < flux_cols.size(); ++j)
            heights[j].push_back(std::stod(tokens[flux_cols[j].first]));
    }
    edges.push_back(std::stod(tokens[ehi_token]));

    const double energy_units = 1.0 / 1.0_GeV;
    std::map<PID, FluxData> result;
    for(size_t j = 0; j < flux_cols.size(); ++j)
        result[flux_cols[j].second] = BuildFluxData(edges, heights[j], true, energy_units, "CSV");
    return result;
}

std::map<PID, FluxData> HepDataFluxLoader::Load(const YAML::Node &entry) const {
    if(entry.IsScalar() || !entry["PID"])
        throw std::runtime_error("Beam::Spectrum: HepData flux entry requires an explicit 'PID'");
    const PID pid{entry["PID"].as<int>()};

    FluxData data;
    data.format = "HepData";
    spdlog::debug("Reading from HepData yaml file");
    YAML::Node root = YAML::LoadFile(FluxEntryPath(entry));
    std::vector<double> heights;
    for(const auto &value : root["dependent_variables"][0]["values"])
        heights.push_back(value["value"].as<double>());

    auto energy_units = root["independent_variables"][0]["header"]["units"].as<std::string>();
    if(energy_units == "GeV") data.energy_units = 1.0 / 1_GeV;
    size_t nbins = root["independent_variables"][0]["values"].size();
    std::vector<double> bin_centers, low, high;
    data.min_energy = root["independent_variables"][0]["values"][0]["low"].as<double>();
    data.max_energy = root["independent_variables"][0]["values"][nbins - 1]["high"].as<double>();
    for(const auto &bin : root["independent_variables"][0]["values"]) {
        low.push_back(bin["low"].as<double>());
        high.push_back(bin["high"].as<double>());
        bin_centers.push_back((low.back() + high.back()) / 2);
    }

    bool use_width = true;
    for(size_t i = 0; i < heights.size(); ++i) {
        double width = use_width ? high[i] - low[i] : 1;
        data.integral += width * heights[i];
    }

    bin_centers.insert(bin_centers.begin(), data.min_energy);
    heights.insert(heights.begin(), heights[0]);
    bin_centers.push_back(data.max_energy);
    heights.push_back(heights.back());

    data.bin_centers = std::move(bin_centers);
    data.heights = std::move(heights);

    std::map<PID, FluxData> result;
    result[pid] = std::move(data);
    return result;
}

std::map<PID, FluxData> RootHistFluxLoader::Load(const YAML::Node &entry) const {
#ifdef USE_ROOT
    if(entry.IsScalar() || !entry["PID"])
        throw std::runtime_error("Beam::Spectrum: ROOT flux entry requires an explicit 'PID'");
    const PID pid{entry["PID"].as<int>()};

    FluxData data;
    data.format = "ROOTHist";
    std::string filename = "flux/" + entry["File"].as<std::string>();
    spdlog::trace("Reading flux file: {}", filename);
    TFile *file = new TFile(filename.c_str());
    TH1D *hist = static_cast<TH1D *>(file->Get(entry["HistName"].as<std::string>().c_str()));
    bool use_width = entry["UseWidth"].as<bool>();
    double norm = entry["Norm"].as<double>();
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

    data.bin_centers = std::move(bin_centers);
    data.heights = std::move(heights);
    data.min_energy = data.bin_centers.front();
    data.max_energy = data.bin_centers.back();
    data.energy_units = 1.0 / 1.0_GeV;

    std::map<PID, FluxData> result;
    result[pid] = std::move(data);
    return result;
#else
    (void)entry;
    throw std::runtime_error("Achilles has not been compiled with ROOT support");
#endif
}

} // namespace achilles
