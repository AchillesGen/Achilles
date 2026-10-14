#include "catch2/catch.hpp"

#include "Achilles/XSecSpline.hh"

#include <cmath>
#include <random>

using achilles::GeometryXSecSpline;
using achilles::XSecSpline;

// The spline returns the cross section sigma(E) estimated as the summed MC
// weight in a bin divided by the TOTAL number of samples thrown (self-tracked as
// the Fill count):
//     sigma(E) = sum_bin / (nsamples * binwidth) = height / Fill-count.
// This is the unbiased MC integral estimator restricted to a bin. It depends on
// the weight sum and the total sample count -- not on how many samples happened
// to land in a given bin -- so it stays correct as VEGAS reshapes the sampling.

namespace {
constexpr double kEMin = 1.0;
constexpr double kEMax = 11.0;
constexpr double kDeltaE = kEMax - kEMin; // 10
constexpr size_t kNBins = 100;
const double kBinWidth = kDeltaE / static_cast<double>(kNBins);

double Center(size_t i) {
    return kEMin + (static_cast<double>(i) + 0.5) * kBinWidth;
}
} // namespace

TEST_CASE("XSecSpline estimates sigma(E) as weight-sum over total samples",
          "[XSecSpline][Geometry]") {
    SECTION("uniform sampling of a constant cross section") {
        // Throw the same number of samples into every bin with weight w.
        // height = K*w/binwidth, nsamples = K*nbins, so sigma = w/DeltaE.
        XSecSpline spline("const", kEMin, kEMax, kNBins);
        const size_t k = 40;
        const double w = 6.0;
        for(size_t i = 0; i < kNBins; ++i)
            for(size_t j = 0; j < k; ++j) spline.Fill(Center(i), w);
        spline.InitializeSpline();
        REQUIRE(spline.Interpolate(5.0) == Approx(w / kDeltaE).margin(1e-9));
        REQUIRE(spline.Interpolate(9.0) == Approx(w / kDeltaE).margin(1e-9));
    }

    SECTION("energy-dependent cross section") {
        // Choose the per-bin weight so the recovered sigma(E) = E:
        // sigma(c) = w_c / DeltaE  =>  w_c = c * DeltaE.
        XSecSpline spline("linear", kEMin, kEMax, kNBins);
        const size_t k = 40;
        for(size_t i = 0; i < kNBins; ++i)
            for(size_t j = 0; j < k; ++j) spline.Fill(Center(i), Center(i) * kDeltaE);
        spline.InitializeSpline();
        REQUIRE(spline.Interpolate(5.0) == Approx(5.0).margin(1e-6));
        REQUIRE(spline.Interpolate(8.0) == Approx(8.0).margin(1e-6));
        // Below the first / above the last bin center the value is linearly
        // extrapolated along the (unit) slope, so it still tracks sigma(E) = E.
        REQUIRE(spline.Interpolate(kEMin) == Approx(kEMin).margin(1e-6));
        REQUIRE(spline.Interpolate(kEMax) == Approx(kEMax).margin(1e-6));
    }

    SECTION("zero-weight samples count toward N and lower sigma") {
        // Same weight sum, but extra zero-weight samples double the sample count
        // (Fill count), so sigma halves. A per-bin-mean estimator would ignore
        // these -- this pins the sum/N behaviour and that zeros are counted.
        XSecSpline spline("zeros", kEMin, kEMax, kNBins);
        for(size_t i = 0; i < kNBins; ++i) {
            for(size_t j = 0; j < 40; ++j) spline.Fill(Center(i), 6.0);
            for(size_t j = 0; j < 40; ++j) spline.Fill(Center(i), 0.0);
        }
        spline.InitializeSpline();
        REQUIRE(spline.Interpolate(5.0) == Approx(0.3).margin(1e-9));
    }
}

TEST_CASE("GeometryXSecSpline sums matching (nu, element) splines", "[XSecSpline][Geometry]") {
    auto make = [](double w) {
        XSecSpline s("m", kEMin, kEMax, kNBins);
        for(size_t i = 0; i < kNBins; ++i)
            for(size_t j = 0; j < 40; ++j) s.Fill(Center(i), w);
        return s;
    };

    GeometryXSecSpline geom;
    geom.AddSpline(14, 6, make(6.0)); // sigma = 0.6
    geom.AddSpline(14, 8, make(2.0)); // sigma = 0.2
    geom.InitializeSplines();

    SECTION("single element") {
        REQUIRE(geom.TotalXSec(5.0, 14, {6}) == Approx(0.6).margin(1e-9));
        REQUIRE(geom.TotalXSec(5.0, 14, {8}) == Approx(0.2).margin(1e-9));
    }

    SECTION("material made of multiple elements sums them") {
        REQUIRE(geom.TotalXSec(5.0, 14, {6, 8}) == Approx(0.8).margin(1e-9));
    }

    SECTION("missing species contributes nothing") {
        REQUIRE(geom.TotalXSec(5.0, -14, {6, 8}) == Approx(0.0).margin(1e-9));
    }
}

TEST_CASE("XSecSpline trades unsupported resolution for precision", "[XSecSpline][Geometry]") {
    // sigma(E) = E, sampled with a large per-bin scatter. The fill granularity is
    // far finer than the statistics support, so the knots must be merged until a
    // bin is known well enough for the interpolant to be worth anything.
    auto fill_noisy = [](XSecSpline &spline, size_t nbins, size_t per_bin, double spread) {
        std::mt19937 rng(1234);
        std::uniform_real_distribution<double> jitter(1 - spread, 1 + spread);
        const double bw = kDeltaE / static_cast<double>(nbins);
        for(size_t i = 0; i < nbins; ++i) {
            const double center = kEMin + (static_cast<double>(i) + 0.5) * bw;
            for(size_t j = 0; j < per_bin; ++j) spline.Fill(center, center * kDeltaE * jitter(rng));
        }
    };

    SECTION("noisy knots are merged down to the target precision") {
        XSecSpline noisy("noisy", kEMin, kEMax, kNBins);
        fill_noisy(noisy, kNBins, 1000, 0.9);
        noisy.InitializeSpline();
        // A single fill bin here is known to ~3.6%, well past the 2% target, so
        // merging must have run and must have reached the target.
        REQUIRE(noisy.BinError() > 0.0);
        REQUIRE(noisy.BinError() < 0.02);
    }

    SECTION("merging suppresses the scatter without biasing sigma(E)") {
        XSecSpline noisy("noisy", kEMin, kEMax, kNBins);
        fill_noisy(noisy, kNBins, 1000, 0.9);
        noisy.InitializeSpline();
        double worst = 0, mean = 0;
        size_t n = 0;
        for(double e = 2.0; e <= 10.0; e += 0.05) {
            const double rel = (noisy.Interpolate(e) - e) / e;
            worst = std::max(worst, std::abs(rel));
            mean += rel;
            ++n;
        }
        // The interpolant must not wander far beyond the precision of the knots
        // it was built from, and must stay unbiased on average.
        CHECK(worst < 4 * noisy.BinError());
        CHECK(std::abs(mean / static_cast<double>(n)) < 0.01);
    }

    SECTION("plenty of statistics keeps the full resolution") {
        // Noiseless fills: every bin is known exactly, so nothing is merged away
        // and a feature one bin wide survives.
        XSecSpline sharp("sharp", kEMin, kEMax, kNBins);
        const size_t k = 5000;
        for(size_t i = 0; i < kNBins; ++i) {
            const double sigma = i == 50 ? 10.0 : 1.0;
            for(size_t j = 0; j < k; ++j) sharp.Fill(Center(i), sigma * kDeltaE);
        }
        sharp.InitializeSpline();
        CHECK(sharp.BinError() == Approx(1.0 / std::sqrt(static_cast<double>(k))).margin(1e-6));
        CHECK(sharp.Interpolate(Center(50)) == Approx(10.0).margin(1e-6));
        CHECK(sharp.Interpolate(Center(48)) == Approx(1.0).margin(1e-6));
    }
}

TEST_CASE("XSecSpline clamps unsupported and unphysical values", "[XSecSpline][Geometry]") {
    SECTION("energies outside the calibrated range give zero cross section") {
        // The generator rejects every trial outside the beam EnergyRange, so
        // the consistent layer-1 value out there is exactly zero, not a linear
        // extrapolation (which can go negative below threshold).
        XSecSpline spline("range", kEMin, kEMax, kNBins);
        for(size_t i = 0; i < kNBins; ++i)
            for(size_t j = 0; j < 40; ++j) spline.Fill(Center(i), Center(i) * kDeltaE);
        spline.InitializeSpline();
        REQUIRE(spline.Interpolate(kEMin - 0.5) == 0.0);
        REQUIRE(spline.Interpolate(kEMax + 0.5) == 0.0);
    }

    SECTION("interpolated sigma is never negative") {
        // A sharp threshold (many empty bins, then a step) makes a natural
        // cubic spline ring below zero near the transition; sigma_tot must
        // clamp at zero instead.
        XSecSpline spline("step", kEMin, kEMax, kNBins);
        for(size_t i = 0; i < kNBins; ++i) {
            const double w = i >= 90 ? 100.0 : 0.0;
            for(size_t j = 0; j < 40; ++j) spline.Fill(Center(i), w);
        }
        spline.InitializeSpline();
        for(double e = kEMin; e <= kEMax; e += kDeltaE / 1000) {
            REQUIRE(spline.Interpolate(e) >= 0.0);
        }
    }
}
