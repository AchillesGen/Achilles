#ifndef FLUX_LOADERS_HH
#define FLUX_LOADERS_HH

#include "Achilles/Factory.hh"
#include "Achilles/ParticleInfo.hh"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace YAML {
class Node;
}

namespace achilles {

// Histogram data for a single species, ready to be turned into an interpolated
// flux spectrum. The bin_centers/heights are already padded with the end points
// so they can be handed directly to Interp1D.
struct FluxData {
    std::vector<double> bin_centers;
    std::vector<double> heights;
    double min_energy{};
    double max_energy{};
    double integral{};
    double energy_units{1};
    std::string format{"Undefined"};
};

// Strategy interface for reading a flux file. A file may contain one or many
// species, so Load returns a map keyed by PID; single-species files yield a map
// of size one. Loaders self-register with Factory<FluxFileLoader> and are
// selected by file extension via LoaderForFile.
class FluxFileLoader {
  public:
    FluxFileLoader() = default;
    FluxFileLoader(const FluxFileLoader &) = default;
    FluxFileLoader(FluxFileLoader &&) = default;
    FluxFileLoader &operator=(const FluxFileLoader &) = default;
    FluxFileLoader &operator=(FluxFileLoader &&) = default;
    virtual ~FluxFileLoader() = default;

    virtual std::map<PID, FluxData> Load(const YAML::Node &entry) const = 0;

    // Category name for the loader Factory (see Factory<FluxFileLoader>).
    static std::string Name() { return "FluxFileLoader"; }
};

// Plain-text histogram fluxes (.dat). The concrete layout (Achilles, MiniBooNE,
// T2K) is auto-detected from the file's first line. The Achilles format is
// self-describing via a "species:" header line; the external formats use a
// hardcoded column convention.
class HistogramFluxLoader : public FluxFileLoader,
                            Registrable<FluxFileLoader, HistogramFluxLoader> {
  public:
    static std::string Name() { return "Histogram"; }
    static std::unique_ptr<FluxFileLoader> Construct() {
        return std::make_unique<HistogramFluxLoader>();
    }
    std::map<PID, FluxData> Load(const YAML::Node &entry) const override;
};

// Comma-separated fluxes (.csv) with a named header row. Flux columns are mapped
// to PIDs by name (see SpeciesNameToPID).
class CsvFluxLoader : public FluxFileLoader, Registrable<FluxFileLoader, CsvFluxLoader> {
  public:
    static std::string Name() { return "Csv"; }
    static std::unique_ptr<FluxFileLoader> Construct() { return std::make_unique<CsvFluxLoader>(); }
    std::map<PID, FluxData> Load(const YAML::Node &entry) const override;
};

// HepData-style YAML flux tables (.yaml/.yml). Single species; the PID is given
// explicitly in the Files entry.
class HepDataFluxLoader : public FluxFileLoader, Registrable<FluxFileLoader, HepDataFluxLoader> {
  public:
    static std::string Name() { return "HepData"; }
    static std::unique_ptr<FluxFileLoader> Construct() {
        return std::make_unique<HepDataFluxLoader>();
    }
    std::map<PID, FluxData> Load(const YAML::Node &entry) const override;
};

// ROOT TH1D histogram fluxes (.root, requires ROOT support). Single species; the
// PID and histogram parameters are given explicitly in the Files entry.
class RootHistFluxLoader : public FluxFileLoader, Registrable<FluxFileLoader, RootHistFluxLoader> {
  public:
    static std::string Name() { return "ROOTHist"; }
    static std::unique_ptr<FluxFileLoader> Construct() {
        return std::make_unique<RootHistFluxLoader>();
    }
    std::map<PID, FluxData> Load(const YAML::Node &entry) const override;
};

// Resolve the file path of a Files entry: either a bare scalar path or a
// { Path: ... } map (used by formats that need extra parameters, e.g. ROOT).
std::string FluxEntryPath(const YAML::Node &entry);

// Pick the loader for a Files entry based on its file extension.
std::unique_ptr<FluxFileLoader> LoaderForFile(const YAML::Node &entry);

// Map a flux-column name (e.g. "numu_flux") to a PID using a fixed table. A
// trailing "_flux" is stripped and the match is case-insensitive.
PID SpeciesNameToPID(const std::string &name);

} // namespace achilles

#endif
