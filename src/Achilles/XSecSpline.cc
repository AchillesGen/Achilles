#include "Achilles/XSecSpline.hh"
#include "Achilles/Histogram.hh"
#include "fmt/ranges.h"
#include "spdlog/spdlog.h"
#include <algorithm>
#include <stdexcept>

using achilles::GeometryXSecSpline;
using achilles::XSecSpline;

XSecSpline::XSecSpline(const std::string &hash, double emin, double emax, size_t nbins)
    : m_emin{emin}, m_emax{emax} {
    m_xsec = Histogram(nbins, emin, emax, hash);
}

void XSecSpline::Fill(double x, double wgt) {
    m_xsec.Fill(x, wgt);
}

// Relative MC error of a bin holding sum(w) and sum(w^2): sqrt(sum w^2)/sum w.
// Heights and sumw2 carry the same 1/binwidth, and the 1/N that turns the bin
// into sigma(E) divides out, so neither appears here.
static double BinRelativeError(double sumw, double sumw2) {
    return sumw > 0 ? std::sqrt(sumw2) / sumw : 0.0;
}

size_t XSecSpline::ChooseMerge(const std::vector<double> &heights,
                               const std::vector<double> &sumw2) {
    // The fill granularity is far finer than the statistics can support: with a
    // few 1e5 warm-up samples each of 200 bins carries ~5% MC error, and a cubic
    // through knots that noisy wanders by that much between them. Merge bins
    // until the typical bin is known to kTargetBinError, which trades resolution
    // sigma(E) does not have for precision it does. The floor keeps enough knots
    // to follow the turn-on; below ~20 the merge starts biasing sigma(E) itself.
    const size_t nbins = heights.size();
    for(size_t merge = 1; merge <= nbins; ++merge) {
        if(nbins % merge != 0) continue;
        if(nbins / merge < kMinBins) break;
        std::vector<double> rel;
        for(size_t i = 0; i < nbins; i += merge) {
            double sw = 0, sw2 = 0;
            for(size_t j = i; j < i + merge; ++j) {
                sw += heights[j];
                sw2 += sumw2[j];
            }
            if(sw > 0) rel.push_back(BinRelativeError(sw, sw2));
        }
        if(rel.empty()) continue;
        std::nth_element(rel.begin(), rel.begin() + static_cast<long>(rel.size() / 2), rel.end());
        if(rel[rel.size() / 2] <= kTargetBinError) return merge;
    }
    return std::max<size_t>(1, nbins / kMinBins);
}

void XSecSpline::InitializeSpline() {
    // Histogram heights are already sum(wgt)/binwidth, so dividing by the total
    // number of samples (the self-tracked Fill count) gives
    // sigma(E) = sum_bin / (nsamples * binwidth).
    const auto raw_heights = m_xsec.Heights();
    const auto raw_centers = m_xsec.Centers();
    const auto raw_sumw2 = m_xsec.SumW2();
    const auto nsamples = m_xsec.Entries();
    const double norm = nsamples > 0 ? 1.0 / static_cast<double>(nsamples) : 0.0;

    const size_t merge = ChooseMerge(raw_heights, raw_sumw2);
    const size_t nbins = raw_heights.size();
    // Group starts. A merge factor that does not divide the fill granularity
    // would leave a short group at the top; the last group instead runs to the
    // end, since dropping it would silently shorten the supported energy range.
    std::vector<size_t> starts;
    for(size_t i = 0; i + merge <= nbins; i += merge) starts.push_back(i);
    if(starts.empty()) starts.push_back(0);

    std::vector<double> heights, centers;
    double rel_sum = 0;
    size_t rel_n = 0;
    for(size_t g = 0; g < starts.size(); ++g) {
        const size_t begin = starts[g];
        const size_t end = g + 1 < starts.size() ? starts[g + 1] : nbins;
        double sw = 0, sw2 = 0, cen = 0;
        for(size_t j = begin; j < end; ++j) {
            sw += raw_heights[j];
            sw2 += raw_sumw2[j];
            cen += raw_centers[j];
        }
        const auto width = static_cast<double>(end - begin);
        // Heights are densities, so a merged bin's density is the mean of the
        // merged ones; its centre is likewise their mean.
        heights.push_back(sw * norm / width);
        centers.push_back(cen / width);
        if(sw > 0) {
            rel_sum += BinRelativeError(sw, sw2);
            ++rel_n;
        }
    }
    m_bin_error = rel_n > 0 ? rel_sum / static_cast<double>(rel_n) : 0;

    // Each bin's sigma estimate belongs at its center, not its edge; interpolate
    // over centers so the spline is not shifted by half a bin. Allow linear
    // extrapolation for energies in the outer half-bins (within the beam range).
    spdlog::debug("XSecSpline: sigma(E) heights [{}]", fmt::join(heights, ", "));
    spdlog::debug("XSecSpline: bin centers [{}]", fmt::join(centers, ", "));
    spdlog::info("XSecSpline {}: {} samples over {} bins (merged {}x), mean bin error {:.2f}%",
                 m_xsec.GetName(), nsamples, heights.size(), merge, 100 * m_bin_error);
    m_interp = Interp1D(centers, heights);
    // Shape-preserving rather than natural cubic: it cannot overshoot between
    // noisy knots, so sigma(E) stays non-negative and does not ring against the
    // hard zero below threshold.
    m_interp.MonotoneSpline();
    m_interp.AllowExtrapolation(true);
    m_initialized = true;
}

double XSecSpline::Interpolate(double energy) const {
    if(!m_initialized)
        throw std::runtime_error(
            "XSecSpline: Ensure the spline is properly created before being used");
    // Outside the calibrated beam range the generator rejects every trial, so
    // the consistent layer-1 cross section is zero (not an extrapolation).
    if(energy < m_emin || energy > m_emax) return 0;
    // The monotone interpolant cannot undershoot between knots, but the linear
    // extrapolation across the outer half-bins still can; a total cross section
    // cannot be negative, so clamp at zero.
    const double value = m_interp(energy);
    if(value < 0) {
        spdlog::trace("XSecSpline: clamping negative interpolated sigma(E = {}) = {} to 0", energy,
                      value);
        return 0;
    }
    return value;
}

bool XSecSpline::SaveState(std::ostream &os) const {
    os << m_initialized << " " << m_emin << " " << m_emax << " ";
    if(m_initialized) {
        m_interp.SaveState(os);
    } else {
        m_xsec.SaveState(os);
    }
    return true;
}

bool XSecSpline::LoadState(std::istream &is) {
    is >> m_initialized >> m_emin >> m_emax;
    if(m_initialized) {
        m_interp.LoadState(is);
        m_interp.AllowExtrapolation(true);
    } else {
        m_xsec.LoadState(is);
    }
    return true;
}

void GeometryXSecSpline::AddSpline(int nu_pid, int elm_pid, const XSecSpline &spline) {
    if(m_xsecs.find({nu_pid, elm_pid}) != m_xsecs.end()) {
        spdlog::warn("GeometryXSecSpline: Overwriting previous spline for nu {} on elm {}", nu_pid,
                     elm_pid);
    }
    m_xsecs[{nu_pid, elm_pid}] = spline;
}

XSecSpline &GeometryXSecSpline::GetSpline(int nu_pid, int elm_pid) {
    if(m_xsecs.find({nu_pid, elm_pid}) == m_xsecs.end()) {
        auto msg =
            fmt::format("GeometryXSecSpline: Invalid spline with neutrino id {} and nucleus id {}",
                        nu_pid, elm_pid);
        throw std::runtime_error(msg);
    }

    return m_xsecs[{nu_pid, elm_pid}];
}

void GeometryXSecSpline::FillSample(double energy, const std::map<int, double> &nu_weights) {
    for(auto &[nuelm, spline] : m_xsecs) {
        auto it = nu_weights.find(nuelm.first);
        spline.Fill(energy, it != nu_weights.end() ? it->second : 0.0);
    }
}

void GeometryXSecSpline::InitializeSplines() {
    for(auto &[nuelm, spline] : m_xsecs) spline.InitializeSpline();
}

double GeometryXSecSpline::TotalXSec(double energy, int nu_pid, std::vector<int> mat_pids) const {
    double total = 0;
    for(const auto elm_pid : mat_pids) {
        if(m_xsecs.find({nu_pid, elm_pid}) == m_xsecs.end()) {
            // Normal during aggregation: a material element simply has no process
            // for this neutrino species, so it contributes nothing to the total.
            spdlog::trace(
                "GeometryXSecSpline: No cross section available for nu {} on elm {}. Returning 0",
                nu_pid, elm_pid);
            continue;
        }
        total += m_xsecs.at({nu_pid, elm_pid}).Interpolate(energy);
    }
    return total;
}

bool GeometryXSecSpline::SaveState(std::ostream &os) const {
    os << m_xsecs.size() << " ";
    for(const auto &[nuelm, spline] : m_xsecs) {
        os << nuelm.first << " " << nuelm.second << " ";
        spline.SaveState(os);
    }

    return true;
}

bool GeometryXSecSpline::LoadState(std::istream &is) {
    size_t nelms;
    is >> nelms;
    for(size_t i = 0; i < nelms; ++i) {
        size_t nu = 0, elm = 0;
        XSecSpline spline;
        is >> nu >> elm;
        spline.LoadState(is);
        m_xsecs[{nu, elm}] = spline;
    }

    return true;
}
