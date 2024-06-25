#include <algorithm>
#include <cmath>
#include <iomanip>

#include "Achilles/AdaptiveMap.hh"
#include "Achilles/Random.hh"
#include "Achilles/Utilities.hh"
#include "fmt/ranges.h"
#include "spdlog/spdlog.h"

using achilles::AdaptiveMap;

bool AdaptiveMap::Deserialize(std::istream &in) {
    in >> m_bins >> m_dims;
    m_hist.resize((m_bins + 1) * m_dims);

    for(auto &x : m_hist) in >> x;

    return true;
}

bool AdaptiveMap::Serialize(std::ostream &out) const {
    out << m_bins << ' ' << m_dims;

    for(const auto &x : m_hist) {
        out << ' ' << std::scientific
            << std::setprecision(std::numeric_limits<double>::max_digits10 - 1) << x;
    }

    return true;
}

size_t AdaptiveMap::FindBin(size_t dim, double x) const {
    const auto edges = Edges(dim);
    auto it = std::lower_bound(edges.begin(), edges.end(), x);
    return static_cast<size_t>(std::distance(edges.begin(), it)) - 1;
}

double AdaptiveMap::operator()(std::vector<double> &rans) {
    double jacobian = 1.0;
    for(std::size_t i = 0; i < m_dims; ++i) {
        const auto position = rans[i] * static_cast<double>(m_bins);
        const auto index = static_cast<size_t>(position);
        const auto loc = position - static_cast<double>(index);
        const double left = lower_edge(i, index);
        const double size = width(i, index);

        // Calculate inverse CDF
        rans[i] = left + loc * size;

        jacobian *= size * static_cast<double>(m_bins);
    }

    return jacobian;
}

double AdaptiveMap::GenerateWeight(const std::vector<double> &rans) const {
    double jacobian = 1.0;
    for(std::size_t i = 0; i < m_dims; ++i) {
        const auto index = FindBin(i, rans[i]);
        jacobian *= width(i, index) * static_cast<double>(m_bins);
    }
    return jacobian;
}

void AdaptiveMap::Adapt(const double &alpha, const std::vector<double> &data) {
    std::vector<double> tmp(m_bins);
    std::vector<double> new_hist(m_hist.size());
    spdlog::trace("Starting Histogram: [{}]", fmt::join(m_hist.begin(), m_hist.end(), ", "));

    for(size_t i = 0; i < m_dims; ++i) {
        // Load data into tmp
        tmp.assign(data.begin() + static_cast<int>(i * m_bins),
                   data.begin() + static_cast<int>((i + 1) * m_bins));

        // Smooth data by averaging over neighbors
        double previous = tmp[0];
        double current = tmp[1];
        tmp[0] = (previous + current) / 2;
        double norm = tmp[0];

        for(size_t bin = 1; bin < m_bins - 1; ++bin) {
            const double sum = previous + current;
            previous = current;
            current = tmp[bin + 1];
            tmp[bin] = (sum + current) / 3;
            norm += tmp[bin];
        }
        tmp.back() = (previous + current) / 2;
        norm += tmp.back();

        // If norm is zero, then there is no data to adjust
        if(norm == 0) continue;

        // Compute the importance factor
        double average_bin = 0;
        for(size_t bin = 0; bin < m_bins; ++bin) {
            if(tmp[bin] != 0) {
                const double r = tmp[bin] / norm;
                const double fac = pow((r - 1.0) / log(r), alpha);
                average_bin += fac;
                tmp[bin] = fac;
            }
        }
        average_bin /= static_cast<double>(m_bins);

        double cbin = 0;
        size_t ibin = 0;

        // Adjust boundaries
        for(size_t nbin = 1; nbin < m_bins; ++nbin) {
            for(; cbin < average_bin; ++ibin) cbin += tmp[ibin];

            const double prev = lower_edge(i, ibin - 1);
            const double curr = lower_edge(i, ibin);

            cbin -= average_bin;
            const double delta = (curr - prev) * cbin;
            new_hist[i * (m_bins + 1) + nbin] = curr - delta / tmp[ibin - 1];
        }
        new_hist[i * (m_bins + 1) + m_bins] = 1.0;
    }

    m_hist = new_hist;
    spdlog::trace("Updated Histogram: [{}]", fmt::join(m_hist.begin(), m_hist.end(), ", "));
}

void AdaptiveMap::Split(achilles::AdaptiveMapSplit split) {
    size_t nsplit{};
    if(split == AdaptiveMapSplit::half) {
        nsplit = 2;
    } else if(split == AdaptiveMapSplit::third) {
        nsplit = 3;
    } else if(split == AdaptiveMapSplit::quarter) {
        nsplit = 4;
    } else {
        throw std::logic_error("Invalid histogram splitting");
    }

    // Create temporary histogram to store new histogram
    std::vector<double> hist(m_dims * (nsplit * m_bins + 1));

    // Split the histogram into new binnings
    for(size_t d = 0; d < m_dims; ++d) {
        for(size_t bin = 0; bin < m_bins; ++bin) {
            const double curr = lower_edge(d, bin);
            const double size = width(d, bin);

            for(size_t i = 0; i < nsplit; ++i) {
                const size_t idx = d * (nsplit * m_bins + 1) + bin * nsplit + i;
                hist[idx] = curr + static_cast<double>(i) * size / static_cast<double>(nsplit);
            }
        }
        // Add the endpoint
        hist[(d + 1) * (nsplit * m_bins) + d] = 1.0;
    }

    // Store the new histogram information
    spdlog::trace("Old Hist: [{}]", fmt::join(m_hist.begin(), m_hist.end(), ", "));
    spdlog::trace("New Hist: [{}]", fmt::join(hist.begin(), hist.end(), ", "));
    m_bins = nsplit * m_bins;
    m_hist = hist;
}
