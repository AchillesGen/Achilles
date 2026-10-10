// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Unweighter.hh"
#include "Achilles/Random.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>

using achilles::ExcessUnweighter;
using achilles::PercentileUnweighter;
using achilles::SortedWeightUnweighter;
using achilles::TailFractionUnweighter;

void SortedWeightUnweighter::AddEvent(double weight) {
    // Do not change cap after pilot run
    if(m_frozen) return;
    const double w = std::abs(weight);
    m_weights.push_back(w);
    if(w > m_max_weight) m_max_weight = w;
    m_dirty = true;
}

void SortedWeightUnweighter::EnsureCap() {
    if(!m_dirty) return;
    m_dirty = false;

    // A process that never fired keeps a zero cap, or ProcessGroup would keep sampling it.
    if(!(m_max_weight > 0.0)) {
        m_cap = 0;
        return;
    }

    std::sort(m_weights.begin(), m_weights.end());
    double total = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);

    m_cap = ComputeCap(m_weights, total);
    // A cap on a zero entry would reject everything.
    if(!(m_cap > 0.0)) m_cap = m_max_weight;
}

double SortedWeightUnweighter::MaxValue() {
    EnsureCap();
    return m_cap;
}

double SortedWeightUnweighter::AcceptEvent(double weight) {
    const double max_wgt = MaxValue();
    m_frozen = true;

    const double prob = std::abs(weight) / max_wgt;
    m_total++;

    if(prob < achilles::Random::Instance().Uniform(0.0, 1.0)) { return 0; }
    m_accepted++;
    const auto sign = weight < 0 ? -1 : 1;
    return prob > 1.0 ? sign * prob : sign * 1.0;
}

void SortedWeightUnweighter::SaveState(std::ostream &os) const {
    const auto default_precision{os.precision()};
    os << std::setprecision(std::numeric_limits<double>::max_digits10 + 1);
    os << RuleTag() << " ";
    os << m_param << " " << m_max_weight << " ";
    os << m_weights.size() << " ";
    for(const auto &w : m_weights) os << w << " ";
    os << m_accepted << " " << m_total << " ";
    os << std::setprecision(static_cast<int>(default_precision));
}

void SortedWeightUnweighter::LoadState(std::istream &is) {
    std::string tag;
    is >> tag;
    if(tag != RuleTag()) {
        auto msg = fmt::format("Unweighter::LoadState: State was written by {}, but loaded into {}",
                               tag, RuleTag());
        throw std::runtime_error(msg);
    }
    is >> m_param >> m_max_weight;
    size_t size;
    is >> size;
    m_weights.resize(size);
    for(auto &w : m_weights) is >> w;
    is >> m_accepted >> m_total;
    m_dirty = true;
    m_frozen = false;
}

PercentileUnweighter::PercentileUnweighter(const YAML::Node &node)
    : SortedWeightUnweighter{node["percentile"].as<double>() / 100} {}

double PercentileUnweighter::ComputeCap(const std::vector<double> &w, double) const {
    size_t idx = static_cast<size_t>(static_cast<double>(w.size()) * m_param);
    if(idx >= w.size()) idx = w.size() - 1;
    return w[idx];
}

std::unique_ptr<achilles::Unweighter> PercentileUnweighter::Construct(const YAML::Node &node) {
    return std::make_unique<PercentileUnweighter>(node);
}

ExcessUnweighter::ExcessUnweighter(const YAML::Node &node)
    : SortedWeightUnweighter{node["epsilon"].as<double>()} {}

double ExcessUnweighter::ComputeCap(const std::vector<double> &w, double total) const {
    const size_t N = w.size();
    const double target = m_param * total;

    std::vector<double> suffix_sum(N + 1, 0);
    for(size_t i = N; i-- > 0;) suffix_sum[i] = suffix_sum[i + 1] + w[i];

    // First w[i] whose excess sum_{j>=i} (w_j - w_i) is within the target; the cap
    // then sits where the excess of the weights above it equals the target.
    size_t i = 0;
    while(i < N && (suffix_sum[i] - static_cast<double>(N - i) * w[i]) > target) ++i;
    return std::clamp((suffix_sum[i] - target) / static_cast<double>(N - i), w.front(), w.back());
}

std::unique_ptr<achilles::Unweighter> ExcessUnweighter::Construct(const YAML::Node &node) {
    return std::make_unique<ExcessUnweighter>(node);
}

TailFractionUnweighter::TailFractionUnweighter(const YAML::Node &node)
    : SortedWeightUnweighter{node["epsilon"].as<double>()} {}

double TailFractionUnweighter::ComputeCap(const std::vector<double> &w, double total) const {
    // The largest weight whose own tail (itself and everything above) exceeds eps*total.
    const double target = m_param * total;
    double tail = 0;
    for(size_t i = w.size(); i-- > 0;) {
        if(tail + w[i] > target) return w[i];
        tail += w[i];
    }
    return w.front();
}

std::unique_ptr<achilles::Unweighter> TailFractionUnweighter::Construct(const YAML::Node &node) {
    return std::make_unique<TailFractionUnweighter>(node);
}
