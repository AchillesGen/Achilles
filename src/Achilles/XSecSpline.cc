#include "Achilles/XSecSpline.hh"
#include "Achilles/Histogram.hh"
#include "fmt/ranges.h"
#include "spdlog/spdlog.h"
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
    auto heights = m_xsec.Heights();
    // Double the first bin to ensure that heights and edges are the same length
    // This should be ok as long as the number of bins is sufficiently large
    // TODO: Determine a better way to handle this
    heights.insert(heights.begin(), heights[0]);
    auto edges = m_xsec.Edges();
    spdlog::info("XSecSpline: Hist heights [{}]", fmt::join(heights, ", "));
    spdlog::info("XSecSpline: Hist edges [{}]", fmt::join(edges, ", "));
    m_interp = Interp1D(edges, heights);
    m_interp.CubicSpline();
    m_initialized = true;
}

double XSecSpline::Interpolate(double energy) const {
    if(!m_initialized)
        throw std::runtime_error(
            "XSecSpline: Ensure the spline is properly created before being used");
    return m_interp(energy);
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

double GeometryXSecSpline::TotalXSec(double energy, int nu_pid, std::vector<int> mat_pids) const {
    double total = 0;
    for(const auto elm_pid : mat_pids) {
        if(m_xsecs.find({nu_pid, elm_pid}) == m_xsecs.end()) {
            spdlog::warn(
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
