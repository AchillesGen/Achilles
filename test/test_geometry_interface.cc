#include "catch2/catch.hpp"

#include "Achilles/AdaptiveMap.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventGen.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Integrand.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Vegas.hh"

#include <cmath>
#include <optional>
#include <string>
#include <vector>

// A geometry driver composes its own layer-1 rejection with one Achilles trial,
// and offers two ways to normalize layer 1: the unweighter envelope
// (VertexEnvelope, single shot) or the sigma(E) spline (TotalXSec, retry until
// emit). The envelope cancels out of the events-per-POT rate, so the two only
// agree when a single trial averages to the total cross section:
//
//     E[event weight | E] * VertexEnvelope(nu, target) == TotalXSec(E, nu, target)
//
// At an injected energy the beam dimension is *not* drawn from the integrand's
// density, so its mapper density and VEGAS jacobian have to be divided back out
// (Integrand::GenerateWeight). Without that the two sides disagree by the
// beam-dimension jacobian, which biases both the rate and the spectrum.

namespace {

// Two-dimensional stub: dim 0 is supplied externally (the injected beam energy),
// dim 1 is sampled from this mapping.
class FixedDimMapper : public achilles::Mapper<double> {
  public:
    FixedDimMapper(double fixed_density, double free_density, bool has_fixed)
        : m_fixed{fixed_density}, m_free{free_density}, m_has_fixed{has_fixed} {}

    void GeneratePoint(std::vector<double> &point, const std::vector<double> &rans) override {
        point = rans;
    }
    double GenerateWeight(const std::vector<double> &point, std::vector<double> &rans) override {
        rans = point;
        return m_fixed * m_free;
    }
    size_t NDims() const override { return 2; }
    std::optional<size_t> FixedDim() const override {
        if(!m_has_fixed) return std::nullopt;
        return 0;
    }
    double FixedDensity() const override { return m_fixed; }

  private:
    double m_fixed, m_free;
    bool m_has_fixed;
};

// Grid with two bins per dimension, each split at its own edge, so the jacobian
// of the first bin is 2*edge and of the second 2*(1 - edge).
achilles::Channel<double> MakeChannel(double fixed_edge, double free_edge, double fixed_density,
                                      double free_density, bool has_fixed) {
    achilles::AdaptiveMap map(2, 2);
    map.Hist() = {0, fixed_edge, 1, 0, free_edge, 1};

    achilles::Channel<double> channel;
    channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
    channel.mapping = std::make_unique<FixedDimMapper>(fixed_density, free_density, has_fixed);
    return channel;
}

} // namespace

TEST_CASE("An externally fixed dimension leaves the multichannel density",
          "[Integrand][Geometry]") {
    // Channel grids skewed differently, so the fixed dimension's jacobian cannot
    // be factored out of the multichannel sum.
    const std::vector<double> fixed_edges{0.25, 0.75};
    const std::vector<double> free_edges{0.4, 0.8};
    const std::vector<double> fixed_densities{3.0, 7.0};
    const std::vector<double> free_densities{2.0, 5.0};
    const std::vector<double> weights{0.4, 0.6};
    const std::vector<double> point{0.1, 0.9}; // bin 0 of dim 0, bin 1 of dim 1
    // Jacobians (bin width * nbins) the point sees in each channel's grid.
    const auto fixed_jacobian = [&](size_t i) { return 2 * fixed_edges[i]; };
    const auto free_jacobian = [&](size_t i) { return 2 * (1 - free_edges[i]); };

    achilles::Integrand<double> integrand;
    for(size_t i = 0; i < 2; ++i)
        integrand.AddChannel(MakeChannel(fixed_edges[i], free_edges[i], fixed_densities[i],
                                         free_densities[i], true));

    std::vector<double> densities(2);
    const double wgt = integrand.GenerateWeight(weights, point, densities);

    // Only the sampled dimension may survive: sum_i alpha_i * rho_free_i / vw_free_i.
    double expected = 0;
    for(size_t i = 0; i < 2; ++i) expected += weights[i] * free_densities[i] / free_jacobian(i);
    CHECK(wgt == Approx(1.0 / expected));

    SECTION("and the result does not depend on that dimension's density or grid") {
        achilles::Integrand<double> other;
        for(size_t i = 0; i < 2; ++i)
            other.AddChannel(MakeChannel(1 - fixed_edges[i], free_edges[i],
                                         11.0 * static_cast<double>(i + 1), free_densities[i],
                                         true));
        std::vector<double> other_densities(2);
        CHECK(other.GenerateWeight(weights, point, other_densities) == Approx(wgt));
    }

    SECTION("while an unfixed dimension is kept") {
        achilles::Integrand<double> sampled;
        for(size_t i = 0; i < 2; ++i)
            sampled.AddChannel(MakeChannel(fixed_edges[i], free_edges[i], fixed_densities[i],
                                           free_densities[i], false));
        std::vector<double> sampled_densities(2);
        double full = 0;
        for(size_t i = 0; i < 2; ++i)
            full += weights[i] * fixed_densities[i] * free_densities[i] /
                    (fixed_jacobian(i) * free_jacobian(i));
        CHECK(sampled.GenerateWeight(weights, point, sampled_densities) == Approx(1.0 / full));
    }

    SECTION("a channel that cannot reach the fixed value drops out") {
        // A zero density there is how an out-of-range injected energy is
        // rejected; it must not turn the sum into a NaN.
        achilles::Integrand<double> blocked;
        blocked.AddChannel(
            MakeChannel(fixed_edges[0], free_edges[0], 0.0, free_densities[0], true));
        blocked.AddChannel(MakeChannel(fixed_edges[1], free_edges[1], fixed_densities[1],
                                       free_densities[1], true));
        std::vector<double> blocked_densities(2);
        const double blocked_wgt = blocked.GenerateWeight(weights, point, blocked_densities);
        CHECK(std::isfinite(blocked_wgt));
        CHECK(blocked_wgt == Approx(1.0 / (weights[1] * free_densities[1] / free_jacobian(1))));
    }
}

namespace {

// Cheapest configuration that still exercises the real path: QE spectral carbon,
// no cascade, geometry beam. Data files resolve against the test working
// directory (the build directory).
const std::string kGeometryRunCard = R"(
Main:
  NEvents: 1
  HardCuts: false
  EventCuts: false
  DoRotate: false
  RunDecays: False
  Output:
      Format: Achilles
      Name: test_geometry_interface_events
      Zipped: True

SherpaOptions: !include "data/default/SherpaOptions.yml"

Processes:
  - Leptons: [14, [13]]

Beam:
  Type: Geometry
  EnergyRange: [100, 10000]
  Species: [14]

Cascade:
  Run: False
  Interactions: !include "data/default/VirtResInteractions.yml"
  Step: 0.04
  Probability: Cylinder
  InMedium: None
  PotentialProp: False
  Algorithm: Base

NuclearModels:
- NuclearModel:
   Model: QESpectral
   ConfigFile: data/info_C12_pke.data
   SpectralP: data/Spectral_Functions/pke12p_tot.data
   SpectralN: data/Spectral_Functions/pke12n_tot.data
   FormFactorFile: "FormFactors.yml"
   Ward: None

Nuclei:
  - Nucleus: !include "data/default/12C.yml"

Options:
  Initialize:
    Seed: 12345678
    Accuracy: 1e-2
  Unweighting:
    Name: Percentile
    percentile: 99

Cache:
  Load: false
  Save: false

Backend:
  Name: Default
  Options: []
)";

} // namespace

TEST_CASE("Geometry layer-1 accessors agree", "[EventGen][Geometry]") {
    static constexpr achilles::PID nu{14}, target{1000060120};
    // Enough trials for a ~1.5% error on the mean weight (acceptance ~15%).
    static constexpr size_t ntrials = 8000;
    // Spanning the beam range matters: the bias this pins follows the VEGAS grid
    // in the beam dimension, so it is largest at the ends and passes through 1
    // in between. The threshold bins of the sigma(E) spline are the least
    // reliable, so start above them.
    const std::vector<double> energies{400, 700, 1000, 1500, 2000, 3000, 4000, 6000, 8000, 9500};

    achilles::EventGen gen(YAML::Load(kGeometryRunCard), {}, "test_geometry_interface");
    gen.Initialize();

    const double envelope = gen.VertexEnvelope(nu, target);
    REQUIRE(envelope > 0);

    std::vector<double> ratios;
    for(const double energy : energies) {
        const achilles::FourVector ray{energy, 0, 0, energy};
        double sum_weights = 0;
        for(size_t i = 0; i < ntrials; ++i) {
            gen.InjectRay(ray, nu, target);
            if(gen.GenerateSingleEvent()) sum_weights += gen.LastEvent().Weight();
        }
        const double envelope_xsec = envelope * sum_weights / static_cast<double>(ntrials);
        const double spline_xsec = gen.TotalXSec(energy, nu, target);
        REQUIRE(spline_xsec > 0);

        INFO("E = " << energy << " MeV: envelope * <w> = " << envelope_xsec
                    << " cm^2, spline = " << spline_xsec << " cm^2");
        CHECK(envelope_xsec == Approx(spline_xsec).epsilon(0.2));
        ratios.push_back(envelope_xsec / spline_xsec);
    }

    // What is left over point to point is the spline's own per-bin noise (it is
    // built from the same warm-up samples), which is uncorrelated in energy and
    // averages out. The bias this pins does not: it is a smooth function of the
    // energy that swings the individual ratios by 25% or more, so it shows up in
    // the spread rather than in the mean, which stays near 1 either way.
    // Measured here: rms 0.057 with the division, 0.129 without.
    double sum = 0, sum_sq = 0;
    for(const double ratio : ratios) {
        sum += ratio;
        sum_sq += (ratio - 1) * (ratio - 1);
    }
    const auto npoints = static_cast<double>(ratios.size());
    const double mean = sum / npoints;
    const double rms = std::sqrt(sum_sq / npoints);
    INFO("mean ratio = " << mean << ", rms deviation from 1 = " << rms);
    CHECK(mean == Approx(1.0).epsilon(0.05));
    CHECK(rms < 0.09);
}
