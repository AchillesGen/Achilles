#include "Achilles/XSecSpline.hh"
#include "Achilles/Histogram.hh"
#include "spdlog/spdlog.h"

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

void GeometryXSecSpline::AddSpline(size_t nu_pid, size_t elm_pid, const XSecSpline &spline) {
    if(m_xsecs.find({nu_pid, elm_pid}) != m_xsecs.end()) {
        spdlog::warn("GeometryXSecSpline: Overwriting previous spline for nu {} on elm {}", nu_pid,
                     elm_pid);
    }
    m_xsecs[{nu_pid, elm_pid}] = spline;
}

double GeometryXSecSpline::TotalXSec(double energy, size_t nu_pid,
                                     std::vector<size_t> mat_pids) const {
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
