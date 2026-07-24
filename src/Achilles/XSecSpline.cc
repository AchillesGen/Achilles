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

void XSecSpline::InitializeSpline() {
    // Histogram heights are already sum(wgt)/binwidth, so dividing by the total
    // number of samples (the self-tracked Fill count) gives
    // sigma(E) = sum_bin / (nsamples * binwidth).
    auto heights = m_xsec.Heights();
    const auto nsamples = m_xsec.Entries();
    const double norm = nsamples > 0 ? 1.0 / static_cast<double>(nsamples) : 0.0;
    for(auto &h : heights) h *= norm;
    // Each bin's sigma estimate belongs at its center, not its edge; interpolate
    // over centers so the spline is not shifted by half a bin. Allow linear
    // extrapolation for energies in the outer half-bins (within the beam range).
    auto centers = m_xsec.Centers();
    spdlog::debug("XSecSpline: sigma(E) heights [{}]", fmt::join(heights, ", "));
    spdlog::debug("XSecSpline: bin centers [{}]", fmt::join(centers, ", "));
    m_interp = Interp1D(centers, heights);
    m_interp.CubicSpline();
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
    // Natural cubic splines undershoot between noisy knots and when linearly
    // extrapolating in the outer half-bins; a total cross section cannot be
    // negative, so clamp at zero.
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
