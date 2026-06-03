#ifndef FLUX_LOADERS_HH
#define FLUX_LOADERS_HH

#include "Achilles/Factory.hh"

#include <memory>
#include <string>
#include <vector>

namespace YAML {
class Node;
}

namespace achilles {

// Histogram data extracted from a flux file, ready to be turned into an
// interpolated flux spectrum. The bin_centers/heights are already padded with
// the end points so they can be handed directly to Interp1D.
struct FluxData {
    std::vector<double> bin_centers;
    std::vector<double> heights;
    double min_energy{};
    double max_energy{};
    double integral{};
    double energy_units{1};
    std::string format{"Undefined"};
};

// Strategy interface for reading a flux spectrum from a particular file format.
// Loaders self-register with Factory<FluxFileLoader> under the YAML key that
// selects them (e.g. "Histogram"), so adding a new format requires no changes
// to Spectrum.
class FluxFileLoader {
  public:
    FluxFileLoader() = default;
    FluxFileLoader(const FluxFileLoader &) = default;
    FluxFileLoader(FluxFileLoader &&) = default;
    FluxFileLoader &operator=(const FluxFileLoader &) = default;
    FluxFileLoader &operator=(FluxFileLoader &&) = default;
    virtual ~FluxFileLoader() = default;

    virtual FluxData Load(const YAML::Node &node) const = 0;

    // Category name for the loader Factory (see Factory<FluxFileLoader>).
    static std::string Name() { return "FluxFileLoader"; }
};

// Plain-text histogram fluxes. The concrete column/header layout (Achilles,
// MiniBooNE, T2K) is auto-detected from the file's first line.
class HistogramFluxLoader : public FluxFileLoader,
                            Registrable<FluxFileLoader, HistogramFluxLoader> {
  public:
    static std::string Name() { return "Histogram"; }
    static std::unique_ptr<FluxFileLoader> Construct() {
        return std::make_unique<HistogramFluxLoader>();
    }
    FluxData Load(const YAML::Node &node) const override;
};

// HepData-style YAML flux tables.
class HepDataFluxLoader : public FluxFileLoader, Registrable<FluxFileLoader, HepDataFluxLoader> {
  public:
    static std::string Name() { return "HepData"; }
    static std::unique_ptr<FluxFileLoader> Construct() {
        return std::make_unique<HepDataFluxLoader>();
    }
    FluxData Load(const YAML::Node &node) const override;
};

// ROOT TH1D histogram fluxes (requires ROOT support).
class RootHistFluxLoader : public FluxFileLoader, Registrable<FluxFileLoader, RootHistFluxLoader> {
  public:
    static std::string Name() { return "ROOTHist"; }
    static std::unique_ptr<FluxFileLoader> Construct() {
        return std::make_unique<RootHistFluxLoader>();
    }
    FluxData Load(const YAML::Node &node) const override;
};

} // namespace achilles

#endif
