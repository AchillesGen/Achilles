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
    spdlog::info("{:3d}   {:^8.5e} +/- {:^8.5e}    {:^8.5e} +/- {:^8.5e}", summary.results.size(),
                 summary.results.back().Mean(), summary.results.back().Error(),
                 summary.Result().Mean(), summary.Result().Error());
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
    os << min_diff;
    os << std::setprecision(static_cast<int>(default_precision));
}

void achilles::MultiChannel::LoadState(std::istream &is) {
    is >> ndims;
    achilles::LoadState(is, params);
    achilles::LoadState(is, summary);
    size_t nweights;
    is >> nweights;
    for(size_t i = 0; i < nweights; ++i) {
        double wgt;
        is >> wgt;
        channel_weights.push_back(wgt);
    }
    for(size_t i = 0; i < nweights; ++i) {
        double wgt;
        is >> wgt;
        best_weights.push_back(wgt);
    }
    is >> min_diff;
}

bool achilles::MultiChannel::operator==(const MultiChannel &rhs) const {
    return ndims == rhs.ndims && params == rhs.params && channel_weights == rhs.channel_weights &&
           best_weights == rhs.best_weights && min_diff == rhs.min_diff && summary == rhs.summary;
}
