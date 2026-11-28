// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/MultiChannel.hh"

bool achilles::operator==(const MultiChannelParams &lhs, const MultiChannelParams &rhs) {
    return lhs.ncalls == rhs.ncalls && lhs.niterations == rhs.niterations && lhs.rtol == rhs.rtol &&
           lhs.nrefine == rhs.nrefine && lhs.beta == rhs.beta && lhs.min_alpha == rhs.min_alpha &&
           lhs.iteration == rhs.iteration;
}

bool achilles::operator==(const MultiChannelSummary &lhs, const MultiChannelSummary &rhs) {
    return lhs.best_weights == rhs.best_weights && lhs.results == rhs.results &&
           lhs.sum_results == rhs.sum_results;
}

void achilles::SaveState(std::ostream &os, const MultiChannelParams &params) {
    os << params.ncalls << " " << params.niterations << " " << params.rtol << " " << params.nrefine
       << " " << params.beta << " " << params.min_alpha << " " << params.iteration;
}

void achilles::LoadState(std::istream &is, MultiChannelParams &params) {
    is >> params.ncalls >> params.niterations >> params.rtol >> params.nrefine >> params.beta >>
        params.min_alpha >> params.iteration;
}

void achilles::SaveState(std::ostream &os, const MultiChannelSummary &summary) {
    summary.sum_results.SaveState(os);
    os << " " << summary.results.size() << " ";
    for(const auto &res : summary.results) {
        res.SaveState(os);
        os << " ";
    }
    os << summary.best_weights.size() << " ";
    for(const auto &wgt : summary.best_weights) { os << wgt << " "; }
}

void achilles::LoadState(std::istream &is, MultiChannelSummary &summary) {
    summary.sum_results.LoadState(is);
    size_t nresults;
    is >> nresults;
    for(size_t i = 0; i < nresults; ++i) {
        StatsData res;
        res.LoadState(is);
        summary.results.push_back(res);
    }

    size_t nweights;
    is >> nweights;
    for(size_t i = 0; i < nweights; ++i) {
        double wgt;
        is >> wgt;
        summary.best_weights.push_back(wgt);
    }
}

achilles::MultiChannel::MultiChannel(size_t dims, size_t nchannels, MultiChannelParams params_,
                                     UnitsEnum unit)
    : ndims{std::move(dims)}, params{std::move(params_)}, m_units(unit) {
    for(size_t i = 0; i < nchannels; ++i) {
        channel_weights.push_back(1.0 / static_cast<double>(nchannels));
    }
}

void achilles::MultiChannel::Adapt(const std::vector<double> &train) {
    std::vector<double> new_weights(channel_weights.size());

    spdlog::debug("MultiChannel::Adapt:");
    double sum_wgts = 0;
    for(size_t i = 0; i < new_weights.size(); ++i) {
        new_weights[i] = channel_weights[i] * pow(train[i], params.beta);
        sum_wgts += new_weights[i];
    }

    if(sum_wgts == 0) return;
    double new_sum = 0;
    for(auto &wgt : new_weights) {
        if(wgt == 0) continue;
        wgt /= sum_wgts;
        wgt = std::max(wgt, params.min_alpha);
        new_sum += wgt;
    }

    size_t idx = 0;
    for(auto &wgt : new_weights) {
        wgt /= new_sum;
        spdlog::debug("  Channel {}: {}", idx++, wgt);
    }

    channel_weights = new_weights;
}

void achilles::MultiChannel::MaxDifference(const std::vector<double> &train) {
    double max = 0;

    for(size_t i = 0; i < train.size() - 1; ++i) {
        const double wi = train[i];
        for(size_t j = i + 1; j < train.size(); ++j) {
            const double wj = train[j];

            max = std::max(max, std::abs(wi - wj));
        }
    }

    if(max < min_diff) {
        best_weights = channel_weights;
        min_diff = max;
    }
}

achilles::MultiChannelSummary achilles::MultiChannel::Summary() {
    summary.best_weights = best_weights;

    double scale = UnitScale(m_units);
    std::string unit = ToString(m_units);

    std::cout << "Final integral = "
              << fmt::format("{:^8.5e} +/- {:^8.5e} {} ( {:^4.2g} %)",
                             summary.Result().Mean() * scale, summary.Result().Error() * scale,
                             unit, summary.Result().Error() / summary.Result().Mean() * 100)
              << std::endl;
    std::cout << "Channel weights:\n";
    for(size_t i = 0; i < best_weights.size(); ++i) {
        std::cout << "  alpha(" << i << ") = " << best_weights[i] << "\n";
    }
    return summary;
}

bool achilles::MultiChannel::NeedsOptimization(double rel_err) const {
    return (rel_err > params.rtol) || summary.results.size() < params.niterations;
}

achilles::OptimizationEstimate achilles::MultiChannel::EstimateRemaining(double rel_err) const {
    // Assume the error falls as 1/ncalls while the grid adapts, or the recent rate if faster
    constexpr size_t max_iterations = 100;
    constexpr size_t window = 3;
    constexpr double fastest_decay = 0.3;

    const bool calls_grow = params.rescale_factor > 1 && params.nrefine > 0;
    double decay =
        calls_grow ? std::pow(params.rescale_factor, -1 / static_cast<double>(params.nrefine)) : 1;
    const size_t nresults_done = summary.results.size();
    if(nresults_done >= 2) {
        const size_t span = std::min(window, nresults_done - 1);
        const double first = summary.results[nresults_done - 1 - span].RelError();
        const double last = summary.results.back().RelError();
        if(first > 0 && last > 0 && std::isfinite(first) && std::isfinite(last)) {
            const double observed = std::pow(last / first, 1 / static_cast<double>(span));
            decay = std::min(decay, std::max(observed, fastest_decay));
        }
    }

    // No usable measurement: assume convergence after the minimum iterations
    const bool known = rel_err < lim::max() && std::isfinite(rel_err);
    double err = (known && decay < 1) ? rel_err : 0;

    OptimizationEstimate estimate;
    size_t nresults = nresults_done;
    size_t iteration = params.iteration;
    size_t ncalls = params.ncalls;
    while((err > params.rtol || nresults < params.niterations) &&
          estimate.iterations < max_iterations) {
        estimate.calls += ncalls;
        ++estimate.iterations;
        ++nresults;
        err *= decay;
        if(++iteration == params.nrefine) {
            iteration = 0;
            ncalls = static_cast<size_t>(static_cast<double>(ncalls) * params.rescale_factor);
        }
    }
    return estimate;
}

std::string achilles::MultiChannel::IterationLine(size_t ncalls, double seconds, double elapsed,
                                                  double left) const {
    // The cumulative result is not shown, but a nan there still means the run is broken
    std::array<double, 4> values = {summary.results.back().Mean(), summary.results.back().Error(),
                                    summary.Result().Mean(), summary.Result().Error()};
    for(double item : values) {
        if(std::isnan(item)) {
            spdlog::error("Unexpected nan encountered in MultiChannel. Aborting.");
            throw std::runtime_error("Encountered nan in MultiChannel.");
        }
    }

    using console::Color;
    using console::Colorize;
    const double scale = UnitScale(m_units);
    const std::string unit = ToString(m_units);
    const std::string pm = console::GetStyle().unicode ? "±" : "+/-";
    const auto result = Colorize(fmt::format("{:.5e} {}", values[0] * scale, unit), Color::blue);
    const auto percent =
        Colorize(fmt::format("{:5.2f}%", std::abs(values[1] / values[0]) * 100), Color::red);
    const auto calls = console::FormatCount(static_cast<double>(ncalls));
    const std::string remaining = left > 0 ? "~" + console::FormatDuration(left) + " left" : "";
    std::string line = fmt::format(
        "{:>3}  {} {} ( {} = {} )  {:>6} calls  {:>6}  ( {} elapsed{} )", summary.results.size(),
        result, pm, Colorize(fmt::format("{:.2e} {}", values[1] * scale, unit), Color::red),
        percent, calls, console::FormatDuration(seconds), console::FormatDuration(elapsed),
        remaining.empty() ? "" : " / " + remaining);
    // Width checked as if time were left, so the last iteration matches the others
    const size_t width = console::RawWidth();
    const size_t longest = console::DisplayWidth(line) + (remaining.empty() ? 16 : 0);
    if(width != 0 && longest > width)
        line = fmt::format(
            "{:>3}  {} {} ( {} = {} )  {:>6} calls  ( {}{} )", summary.results.size(), result, pm,
            Colorize(fmt::format("{:.2e} {}", values[1] * scale, unit), Color::red), percent, calls,
            console::FormatDuration(elapsed), remaining.empty() ? "" : " / " + remaining);
    return line;
}

void achilles::MultiChannel::SaveState(std::ostream &os) const {
    const auto default_precision{os.precision()};
    os << std::setprecision(std::numeric_limits<double>::max_digits10 + 1);
    os << ndims << " ";
    achilles::SaveState(os, params);
    os << " ";
    achilles::SaveState(os, summary);
    os << " " << channel_weights.size() << " ";
    for(const auto &wgt : channel_weights) { os << wgt << " "; }
    for(const auto &wgt : best_weights) { os << wgt << " "; }
    os << min_diff << " ";
    os << std::setprecision(static_cast<int>(default_precision));
}

void achilles::MultiChannel::LoadState(std::istream &is) {
    is >> ndims;
    achilles::LoadState(is, params);
    achilles::LoadState(is, summary);
    size_t nweights;
    is >> nweights;
    channel_weights.resize(nweights);
    best_weights.resize(nweights);
    for(size_t i = 0; i < nweights; ++i) {
        double wgt;
        is >> wgt;
        channel_weights[i] = wgt;
    }
    for(size_t i = 0; i < nweights; ++i) {
        double wgt;
        is >> wgt;
        best_weights[i] = wgt;
    }
    is >> min_diff;
}

bool achilles::MultiChannel::operator==(const MultiChannel &rhs) const {
    return ndims == rhs.ndims && params == rhs.params && channel_weights == rhs.channel_weights &&
           best_weights == rhs.best_weights && min_diff == rhs.min_diff && summary == rhs.summary;
}
