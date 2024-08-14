#ifndef MULTICHANNEL_HH
#define MULTICHANNEL_HH

#include "Achilles/Integrand.hh"
#include "Achilles/Vegas.hh"

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
    size_t iteration{};

    static constexpr size_t ncalls_default{1000}, nint_default{10};
    static constexpr double rtol_default{1e-2};
    static constexpr size_t nrefine_default{1};
    static constexpr double beta_default{0.25}, min_alpha_default{1e-5};
    static constexpr size_t nparams = 7;
};

bool operator==(const MultiChannelParams &lhs, const MultiChannelParams &rhs);
bool operator==(const MultiChannelSummary &lhs, const MultiChannelSummary &rhs);

void SaveState(std::ostream &os, const MultiChannelParams &params);
void LoadState(std::istream &is, MultiChannelParams &params);
void SaveState(std::ostream &os, const MultiChannelSummary &summary);
void LoadState(std::istream &is, MultiChannelSummary &summary);

class MultiChannel {
  public:
    MultiChannel() = default;
    MultiChannel(size_t, size_t, MultiChannelParams);

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
    bool NeedsOptimization(double) const;

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
        params.ncalls = static_cast<size_t>(static_cast<double>(params.ncalls) * 1.2);
        for(auto &channel : func.Channels()) {
            if(channel.integrator.Grid().Bins() < 200) channel.integrator.Refine();
        }
    }
    void PrintIteration() const;
    void MaxDifference(const std::vector<double> &);

    bool update_summary{true};
    size_t ndims{};
    MultiChannelParams params{};
    std::vector<double> channel_weights, best_weights;
    double min_diff{lim::infinity()};
    MultiChannelSummary summary;
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
    double rel_err = lim::max();
    while(NeedsOptimization(rel_err)) {
        (*this)(func);
        StatsData current = summary.results.back();
        rel_err = current.Error() / std::abs(current.Mean());

        PrintIteration();
        if(++params.iteration == params.nrefine) RefineChannels(func);
    }
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

        return node;
    }

    static bool decode(const Node &node, achilles::MultiChannelParams &rhs) {
        if(node.size() != rhs.nparams) return false;

        rhs.ncalls = node["NCalls"].as<size_t>();
        rhs.niterations = node["NIterations"].as<size_t>();
        rhs.rtol = node["rtol"].as<double>();
        rhs.nrefine = node["nrefine"].as<size_t>();
        rhs.beta = node["beta"].as<double>();
        rhs.min_alpha = node["min_alpha"].as<double>();
        rhs.iteration = node["iteration"].as<size_t>();

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
