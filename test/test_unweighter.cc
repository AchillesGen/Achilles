// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_template_test_macros.hpp"

#include "Achilles/Unweighter.hh"

#include <algorithm>
#include <numeric>
#include <vector>

using achilles::ExcessUnweighter;
using achilles::PercentileUnweighter;
using achilles::TailFractionUnweighter;

namespace {

// A long-tailed sample: the tail is the whole point of the cap rules.
std::vector<double> TailedWeights(size_t n) {
    std::vector<double> weights(n);
    for(size_t i = 0; i < n; ++i) { weights[i] = 1e-4 * static_cast<double>(i * i + 1); }
    return weights;
}

double Total(const std::vector<double> &weights) {
    return std::accumulate(weights.begin(), weights.end(), 0.0);
}

} // namespace

TEST_CASE("Percentile cap is the empirical percentile", "[Unweighter]") {
    PercentileUnweighter unweighter(YAML::Load("percentile: 90"));
    auto weights = TailedWeights(1000);
    for(const auto &w : weights) unweighter.AddEvent(w);

    std::sort(weights.begin(), weights.end());
    CHECK(unweighter.MaxValue() == weights[900]);
}

TEST_CASE("Excess cap bounds the weight above it", "[Unweighter]") {
    const double eps = 0.01;
    ExcessUnweighter unweighter(YAML::Load("epsilon: 0.01"));
    const auto weights = TailedWeights(1000);
    for(const auto &w : weights) unweighter.AddEvent(w);

    // The cap solves the bound exactly, so the check needs room for the two sums
    // being accumulated in a different order.
    const double target = eps * Total(weights);
    const double cap = unweighter.MaxValue();
    double excess = 0;
    for(const auto &w : weights) excess += std::max(w - cap, 0.0);
    CHECK(excess <= target * (1 + 1e-12));
    // ...and it is the smallest such cap: a lower one breaks the bound.
    double excess_lower = 0;
    for(const auto &w : weights) excess_lower += std::max(w - 0.95 * cap, 0.0);
    CHECK(excess_lower > target);
}

TEST_CASE("TailFraction cap bounds the weight in the tail", "[Unweighter]") {
    const double eps = 0.01;
    TailFractionUnweighter unweighter(YAML::Load("epsilon: 0.01"));
    const auto weights = TailedWeights(1000);
    for(const auto &w : weights) unweighter.AddEvent(w);

    const double cap = unweighter.MaxValue();
    double tail = 0;
    for(const auto &w : weights) {
        if(w > cap) tail += w;
    }
    CHECK(tail <= eps * Total(weights));
}

// A process with no allowed states is offered nothing but zeros. Its cap has to stay
// at zero: ProcessGroup turns the max weight into a selection probability, so any
// non-zero fallback makes the generator spend nearly every trial on a channel that
// can only ever return zero.
TEMPLATE_TEST_CASE("A process that never fires keeps a zero cap", "[Unweighter]",
                   PercentileUnweighter, ExcessUnweighter, TailFractionUnweighter) {
    TestType unweighter(YAML::Load("percentile: 99\nepsilon: 0.01"));
    for(size_t i = 0; i < 100; ++i) unweighter.AddEvent(0);
    CHECK(unweighter.MaxValue() == 0);

    TestType never_offered(YAML::Load("percentile: 99\nepsilon: 0.01"));
    CHECK(never_offered.MaxValue() == 0);
}

TEST_CASE("Unweighting keeps the excess above the cap", "[Unweighter]") {
    PercentileUnweighter unweighter(YAML::Load("percentile: 90"));
    for(const auto &w : TailedWeights(1000)) unweighter.AddEvent(w);
    const double cap = unweighter.MaxValue();

    // An overweight event is always accepted, carrying weight |w|/cap > 1, which is
    // what keeps the capped sample unbiased.
    const double overweight = 4 * cap;
    CHECK(unweighter.AcceptEvent(overweight) == overweight / cap);
    CHECK(unweighter.AcceptEvent(-overweight) == -overweight / cap);

    // Freezing on the first accept: a later event must not move the cap.
    unweighter.AddEvent(1e6);
    CHECK(unweighter.MaxValue() == cap);
}
