#include "Achilles/MultiChannel.hh"

achilles::MultiChannel::MultiChannel(size_t dims, size_t nchannels, MultiChannelParams params_)
    : ndims{std::move(dims)}, params{std::move(params_)} {
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
    std::cout << "Final integral = "
              << fmt::format("{:^8.5e} +/- {:^8.5e} ({:^8.5e} %)", summary.Result().Mean(),
                             summary.Result().Error(),
                             summary.Result().Error() / summary.Result().Mean() * 100)
              << std::endl;
    std::cout << "Channel weights:\n";
    for(size_t i = 0; i < best_weights.size(); ++i) {
        std::cout << "  alpha(" << i << ") = " << best_weights[i] << "\n";
    }
    return summary;
}

void achilles::MultiChannel::PrintIteration() const {
    std::cout << fmt::format("{:3d}   {:^8.5e} +/- {:^8.5e}    {:^8.5e} +/- {:^8.5e}",
                             summary.results.size(), summary.results.back().Mean(),
                             summary.results.back().Error(), summary.Result().Mean(),
                             summary.Result().Error())
              << std::endl;
}
