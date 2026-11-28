// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MULTICHANNEL_HH
#define MULTICHANNEL_HH

#include "Achilles/Console.hh"
#include "Achilles/Integrand.hh"
#include "Achilles/Progress.hh"
#include "Achilles/Vegas.hh"

#include <functional>

namespace achilles {

struct MultiChannelSummary {
    std::vector<StatsData> results;
    std::vector<double> best_weights;
    StatsData sum_results;

    StatsData Result() const { return sum_results; }
    StatsData LastResult() const { return results.back(); }
};

struct MultiChannelParams {
    size_t ncalls{ncalls_default}, niterations{nint_default};
    double rtol{rtol_default};
    size_t nrefine{nrefine_default};
    double beta{beta_default}, min_alpha{min_alpha_default};
    double rescale_factor{rescale_factor_default};
    size_t iteration{};

    static constexpr size_t ncalls_default{10000}, nint_default{7};
    static constexpr double rtol_default{1e-2};
    static constexpr size_t nrefine_default{1};
    static constexpr double beta_default{1.00}, min_alpha_default{1e-5};
    static constexpr double rescale_factor_default{1.44};
    static constexpr size_t nparams = 8;
};

struct OptimizationEstimate {
    size_t calls{}, iterations{};
};

bool operator==(const MultiChannelParams &lhs, const MultiChannelParams &rhs);
bool operator==(const MultiChannelSummary &lhs, const MultiChannelSummary &rhs);

void SaveState(std::ostream &os, const MultiChannelParams &params);
void LoadState(std::istream &is, MultiChannelParams &params);
void SaveState(std::ostream &os, const MultiChannelSummary &summary);
void LoadState(std::istream &is, MultiChannelSummary &summary);

enum class UnitsEnum { mub, nb, pb, em38cm2, LAST };

inline std::string ToString(UnitsEnum unit) {
    switch(unit) {
    case UnitsEnum::mub:
        return "mub";
    case UnitsEnum::nb:
        return "nb";
    case UnitsEnum::pb:
        return "pb";
    case UnitsEnum::em38cm2:
        return "e-38cm2";
    default:
        return "nb";
    }
}

inline double UnitScale(UnitsEnum unit) {
    switch(unit) {
    case UnitsEnum::mub:
        return 1e-3;
    case UnitsEnum::nb:
        return 1;
    case UnitsEnum::pb:
        return 1e3;
    case UnitsEnum::em38cm2:
        return 1e5;
    default:
        return 1.0;
    }
}

inline UnitsEnum ToUnitEnum(std::string unit) {
    if(unit == "mub")
        return UnitsEnum::mub;
    else if(unit == "nb")
        return UnitsEnum::nb;
    else if(unit == "pb")
        return UnitsEnum::pb;
    else if(unit == "e-38cm2")
        return UnitsEnum::em38cm2;
    else {
        spdlog::warn("Unit {} not recognized. Available display units include:", unit);
        for(int i = 0; i < static_cast<int>(UnitsEnum::LAST); ++i) {
            UnitsEnum unit_avail = static_cast<UnitsEnum>(i);
            spdlog::warn("{}", ToString(unit_avail));
        }
        spdlog::warn("Defaulting back to nb");
        return UnitsEnum::nb;
    }
}

class MultiChannel {
  public:
    MultiChannel() = default;
    MultiChannel(size_t, size_t, MultiChannelParams, UnitsEnum unit = UnitsEnum::nb);

    // Utilities
    size_t Dimensions() const { return ndims; }
    size_t NChannels() const { return channel_weights.size(); }
    MultiChannelParams Parameters() const { return params; }
    MultiChannelParams &Parameters() { return params; }

    // Optimization and event generation
    template <typename T> void operator()(Integrand<T> &);
    template <typename T> void Optimize(Integrand<T> &);
    template <typename T> std::vector<T> GeneratePoint(Integrand<T> &);
    template <typename T> double GenerateWeight(Integrand<T> &, const std::vector<T> &);

    // Getting results
    bool HasResults() const { return summary.results.size() != 0; }
    void UpdateSummary(bool update) { update_summary = update; }
    MultiChannelSummary Summary();
    StatsData LastResult() const { return summary.LastResult(); }
    StatsData Result() const { return summary.Result(); }
    size_t NIterations() const { return summary.results.size(); }
    bool NeedsOptimization(double) const;
    // rel_err is lim::max() before the first iteration
    OptimizationEstimate EstimateRemaining(double rel_err) const;

    void SetProgressCallback(std::function<void(size_t)> callback) {
        m_progress = std::move(callback);
    }
    static constexpr size_t progress_stride = 1024;

    // Cache Results
    void SaveState(std::ostream &) const;
    void LoadState(std::istream &);
    bool operator==(const MultiChannel &) const;

    // YAML interface
    friend YAML::convert<achilles::MultiChannel>;

  private:
    void Adapt(const std::vector<double> &);
    void TrainChannels();
    template <typename T> void RefineChannels(Integrand<T> &func) {
        params.iteration = 0;
        params.ncalls =
            static_cast<size_t>(static_cast<double>(params.ncalls) * params.rescale_factor);
        for(auto &channel : func.Channels()) {
            if(channel.integrator.Grid().Bins() < 100) channel.integrator.Refine();
        }
    }
    std::string IterationLine(size_t ncalls, double seconds, double elapsed, double left) const;
    void MaxDifference(const std::vector<double> &);

    bool update_summary{true};
    size_t ndims{};
    MultiChannelParams params{};
    std::vector<double> channel_weights, best_weights;
    double min_diff{lim::infinity()};
    MultiChannelSummary summary;
    UnitsEnum m_units = UnitsEnum::nb;
    std::function<void(size_t)> m_progress;
};

template <typename T> void achilles::MultiChannel::operator()(Integrand<T> &func) {
    size_t nchannels = channel_weights.size();
    std::vector<double> rans(ndims);
    std::vector<T> point(ndims);
    std::vector<double> densities(nchannels);
    std::vector<double> train_data(nchannels);

    StatsData results;
    func.InitializeTrain();

    for(size_t i = 0; i < params.ncalls; ++i) {
        if(m_progress && i % progress_stride == 0) m_progress(i);

        // Generate needed random numbers
        Random::Instance().Generate(rans);

        // Select a channel
        size_t ichannel = Random::Instance().SelectIndex(channel_weights);

        // Map the point based on the channel
        func.GeneratePoint(ichannel, rans, point);

        // Evaluate the function at this point
        double wgt = func.GenerateWeight(channel_weights, point, densities);
        double val = wgt == 0 ? 0 : func(point, wgt);
        double val2 = val * val;
        func.AddTrainData(ichannel, val2);
        results += val;

        if(val2 != 0) {
            for(size_t j = 0; j < nchannels; ++j) { train_data[j] += densities[j] * val2 * wgt; }
        }
    }

    Adapt(train_data);
    func.Train();
    MaxDifference(train_data);
    if(update_summary) {
        summary.results.push_back(results);
        summary.sum_results += results;
    }
}

template <typename T> std::vector<T> achilles::MultiChannel::GeneratePoint(Integrand<T> &func) {
    std::vector<double> rans(ndims);
    std::vector<T> point(ndims);

    // Generate needed random numbers
    Random::Instance().Generate(rans);

    // Select a channel
    size_t ichannel = Random::Instance().SelectIndex(channel_weights);

    // Map the point based on the channel
    func.GeneratePoint(ichannel, rans, point);
    return point;
}

template <typename T>
double achilles::MultiChannel::GenerateWeight(Integrand<T> &func, const std::vector<T> &point) {
    size_t nchannels = channel_weights.size();
    std::vector<double> densities(nchannels);
    return func.GenerateWeight(channel_weights, point, densities);
}

template <typename T> void achilles::MultiChannel::Optimize(Integrand<T> &func) {
    // Do not leave a callback to the local bar behind if an iteration throws
    struct ProgressScope {
        MultiChannel &integrator;
        ~ProgressScope() { integrator.m_progress = nullptr; }
    } scope{*this};

    double rel_err = lim::max();

    auto estimate = EstimateRemaining(rel_err);
    ProgressBar bar("Optimizing", static_cast<double>(estimate.calls));
    bar.SetTotal(static_cast<double>(estimate.calls), true);
    // The iteration lines already show the progress
    bar.SetLogLines(false);
    double calls_done = 0;
    size_t expected_iterations = summary.results.size() + estimate.iterations;
    m_progress = [&](size_t calls) {
        bar.Update(
            calls_done + static_cast<double>(calls),
            fmt::format("iteration {}/~{}", summary.results.size() + 1, expected_iterations));
    };

    Stopwatch iteration_clock;
    while(NeedsOptimization(rel_err)) {
        iteration_clock.Reset();
        const size_t ncalls = params.ncalls;
        (*this)(func);
        calls_done += static_cast<double>(ncalls);
        RunTracker::Instance().SetSecondsPerCall(bar.Elapsed() / calls_done);

        StatsData current = summary.results.back();
        rel_err = current.Error() / std::abs(current.Mean());
        const double seconds = iteration_clock.Elapsed();

        if(++params.iteration == params.nrefine) RefineChannels(func);
        estimate = EstimateRemaining(rel_err);
        expected_iterations = summary.results.size() + estimate.iterations;
        bar.SetTotal(calls_done + static_cast<double>(estimate.calls), true);

        const double left = NeedsOptimization(rel_err)
                                ? static_cast<double>(estimate.calls) * bar.Elapsed() / calls_done
                                : 0;
        console::Print(IterationLine(ncalls, seconds, bar.Elapsed(), left));
    }
    bar.Finish();
}

} // namespace achilles

namespace YAML {

template <> struct convert<achilles::MultiChannelSummary> {
    static Node encode(const achilles::MultiChannelSummary &rhs) {
        Node node;
        node["NEntries"] = rhs.results.size();
        for(const auto &entry : rhs.results) node["Entries"].push_back(entry);

        node["NChannels"] = rhs.best_weights.size();
        for(const auto &weight : rhs.best_weights) node["ChannelWeights"].push_back(weight);

        return node;
    }

    static bool decode(const Node &node, achilles::MultiChannelSummary &rhs) {
        // Get the number of entries and ensure that is the number of entries
        auto nentries = node["NEntries"].as<size_t>();
        if(node["Entries"].size() != nentries) return false;

        // Load the entries and keep track of the sum
        for(const auto &entry : node["Entries"]) {
            rhs.results.push_back(entry.as<achilles::StatsData>());
            rhs.sum_results += rhs.results.back();
        }

        // Get the number of channels and ensure that is the number of channels
        auto nchannels = node["NChannels"].as<size_t>();
        if(node["ChannelWeights"].size() != nchannels) return false;

        // Load the best weights
        for(const auto &weight : node["ChannelWeights"])
            rhs.best_weights.push_back(weight.as<double>());

        return true;
    }
};

template <> struct convert<achilles::MultiChannelParams> {
    static Node encode(const achilles::MultiChannelParams &rhs) {
        Node node;

        node["NCalls"] = rhs.ncalls;
        node["NIterations"] = rhs.niterations;
        node["rtol"] = rhs.rtol;
        node["nrefine"] = rhs.nrefine;
        node["beta"] = rhs.beta;
        node["min_alpha"] = rhs.min_alpha;
        node["iteration"] = rhs.iteration;
        node["rescale_factor"] = rhs.rescale_factor;

        return node;
    }

    static bool decode(const Node &node, achilles::MultiChannelParams &rhs) {
        if(node["NCalls"]) rhs.ncalls = node["NCalls"].as<size_t>();
        if(node["NIterations"]) rhs.niterations = node["NIterations"].as<size_t>();
        if(node["rtol"]) rhs.rtol = node["rtol"].as<double>();
        if(node["nrefine"]) rhs.nrefine = node["nrefine"].as<size_t>();
        if(node["beta"]) rhs.beta = node["beta"].as<double>();
        if(node["min_alpha"]) rhs.min_alpha = node["min_alpha"].as<double>();
        if(node["iteration"]) rhs.iteration = node["iteration"].as<size_t>();
        if(node["rescale_factor"]) rhs.rescale_factor = node["rescale_factor"].as<double>();

        return true;
    }
};

template <> struct convert<achilles::MultiChannel> {
    static Node encode(const achilles::MultiChannel &rhs) {
        Node node;
        node["NDims"] = rhs.ndims;
        node["NChannels"] = rhs.best_weights.size();
        node["Parameters"] = rhs.params;
        node["Summary"] = rhs.summary;
        return node;
    }

    static bool decode(const Node &node, achilles::MultiChannel &rhs) {
        if(node.size() != 4) return false;

        rhs.ndims = node["NDims"].as<size_t>();
        rhs.summary = node["Summary"].as<achilles::MultiChannelSummary>();
        rhs.params = node["Parameters"].as<achilles::MultiChannelParams>();

        auto nchannels = node["NChannels"].as<size_t>();
        if(rhs.summary.best_weights.size() != nchannels) return false;
        rhs.channel_weights = rhs.summary.best_weights;
        rhs.best_weights = rhs.summary.best_weights;
        return true;
    }
};

} // namespace YAML

#endif
