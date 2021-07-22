#ifndef MULTICHANNEL_HH
#define MULTICHANNEL_HH

#include "nuchic/Vegas.hh"
#include "nuchic/Integrand.hh"

namespace nuchic {

struct MultiChannelSummary {
    std::vector<StatsData> results;
    std::vector<double> best_weights;
    StatsData sum_results;

    StatsData Result() const { return sum_results; }
};

struct MultiChannelParams {
    size_t ncalls{ncalls_default}, niterations{nint_default};
    double atol{atol_default}, rtol{rtol_default};
    size_t nrefine{nrefine_default};
    double beta{beta_default}, min_alpha{min_alpha_default};
    size_t iteration{};

    static constexpr size_t ncalls_default{10000}, nint_default{10};
    static constexpr double atol_default{1e-4}, rtol_default{1e-4};
    static constexpr size_t nrefine_default{5};
    static constexpr double beta_default{0.25}, min_alpha_default{1e-5};
    static constexpr size_t nparams = 8;
};

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
        template<typename T>
        void operator()(Integrand<T>&);
        template<typename T>
        void Optimize(Integrand<T>&);

        // Getting results
        MultiChannelSummary Summary();

        // YAML interface
        friend YAML::convert<nuchic::MultiChannel>;

    private:
        void Adapt(const std::vector<double>&);
        void TrainChannels();
        template<typename T>
        void RefineChannels(Integrand<T> &func) {
            params.iteration = 0;
            params.nrefine *= 2;
            for(auto &channel : func.Channels())
                channel.integrator.Refine();
        }
        void PrintIteration() const;
        void MaxDifference(const std::vector<double>&);

        size_t ndims{};
        MultiChannelParams params{};
        std::vector<double> channel_weights, best_weights;
        double min_diff{lim::infinity()};
        MultiChannelSummary summary;
};

template<typename T>
void nuchic::MultiChannel::operator()(Integrand<T> &func) {
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
        double val = func(point, wgt);
        double val2 = val * val;
        func.AddTrainData(ichannel, val2);
        results += val;

        for(size_t j = 0; j < nchannels; ++j) {
            train_data[j] += densities[j] * val2 * wgt;
        }
    }

    Adapt(train_data);
    func.Train();
    MaxDifference(train_data);
    summary.results.push_back(results);
    summary.sum_results += results;
}

template<typename T>
void nuchic::MultiChannel::Optimize(Integrand<T> &func) {
    double abs_err = lim::max(), rel_err = lim::max();
    while((abs_err > params.atol && rel_err > params.rtol) || summary.results.size() < params.niterations) {
        (*this)(func);
        StatsData current = summary.Result();
        abs_err = current.Error();
        rel_err = abs_err / std::abs(current.Mean());

        PrintIteration();
        if(params.iteration++ == params.nrefine) RefineChannels(func);
    }
}

}

namespace YAML {

template<>
struct convert<nuchic::MultiChannelSummary> {
    static Node encode(const nuchic::MultiChannelSummary &rhs) {
        Node node;
        node["NEntries"] = rhs.results.size();
        for(const auto &entry : rhs.results)
            node["Entries"].push_back(entry);

        node["NChannels"] = rhs.best_weights.size();
        for(const auto &weight : rhs.best_weights)
            node["ChannelWeights"].push_back(weight);

        return node;
    }

    static bool decode(const Node &node, nuchic::MultiChannelSummary &rhs) {
        // Get the number of entries and ensure that is the number of entries
        auto nentries = node["NEntries"].as<size_t>();
        if(node["Entries"].size() != nentries) return false;

        // Load the entries and keep track of the sum
        for(const auto &entry : node["Entries"]) {
            rhs.results.push_back(entry.as<nuchic::StatsData>());
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

template<>
struct convert<nuchic::MultiChannelParams> {
    static Node encode(const nuchic::MultiChannelParams &rhs) {
        Node node;

        node["NCalls"] = rhs.ncalls;
        node["NIterations"] = rhs.niterations;
        node["atol"] = rhs.atol;
        node["rtol"] = rhs.rtol;
        node["nrefine"] = rhs.nrefine;
        node["beta"] = rhs.beta;
        node["min_alpha"] = rhs.min_alpha;
        node["iteration"] = rhs.iteration;

        return node;
    }

    static bool decode(const Node &node, nuchic::MultiChannelParams &rhs) {
        if(node.size() != rhs.nparams) return false;

        rhs.ncalls = node["NCalls"].as<size_t>();
        rhs.niterations = node["NIterations"].as<size_t>();
        rhs.atol = node["atol"].as<double>();
        rhs.rtol = node["rtol"].as<double>();
        rhs.nrefine = node["nrefine"].as<size_t>();
        rhs.beta = node["beta"].as<double>();
        rhs.min_alpha = node["min_alpha"].as<double>();
        rhs.iteration = node["iteration"].as<size_t>();

        return true;
    }
};

template<>
struct convert<nuchic::MultiChannel> {
    static Node encode(const nuchic::MultiChannel &rhs) {
        Node node;
        node["NDims"] = rhs.ndims;
        node["NChannels"] = rhs.best_weights.size();
        node["Parameters"] = rhs.params;
        node["Summary"] = rhs.summary;
        return node;
    }

    static bool decode(const Node &node, nuchic::MultiChannel &rhs) {
        if(node.size() != 4) return false; 

        rhs.ndims = node["NDims"].as<size_t>();
        rhs.summary = node["Summary"].as<nuchic::MultiChannelSummary>();
        rhs.params = node["Parameters"].as<nuchic::MultiChannelParams>();

        auto nchannels = node["NChannels"].as<size_t>();
        if(rhs.summary.best_weights.size() != nchannels) return false; 
        rhs.channel_weights = rhs.summary.best_weights;
        rhs.best_weights = rhs.summary.best_weights;
        return true;
    }
};

}

#endif
