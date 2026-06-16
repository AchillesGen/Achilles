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

void achilles::MultiChannel::PrintIteration() const {
    std::array<double, 4> values = {summary.results.back().Mean(), summary.results.back().Error(),
                                    summary.Result().Mean(), summary.Result().Error()};
    for(double item : values) {
        if(std::isnan(item)) {
            spdlog::error("Unexpected nan encountered in MultiChannel. Aborting.");
            throw std::runtime_error("Encountered nan in MultiChannel.");
        }
    }

    double scale = UnitScale(m_units);
    std::string unit = ToString(m_units);

    spdlog::info(
        "{:3d}   {:^8.5e} +/- {:^8.5e} {} ( {:^4.2g} %)    {:^8.5e} +/- {:^8.5e} {} ( {:^4.2g} %)",
        summary.results.size(), values[0] * scale, values[1] * scale, unit,
        values[1] / values[0] * 100, values[2] * scale, values[3] * scale, unit,
        values[3] / values[2] * 100);
}

void achilles::MultiChannel::SaveState(std::ostream &os) const {
    const auto default_precision{os.precision()};
    os << std::setprecision(std::numeric_limits<double>::max_digits10 + 1);
    os << static_cast<int>(IntegratorType::MultiChannel) << " ";
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

bool achilles::MultiChannel::LoadState(std::istream &is) {
    // A missing cache, or one written by a different integrator backend, is not an
    // error: leave this integrator empty so it reports NeedsOptimization() and gets
    // optimized from scratch. Do not consume the rest of the stream in that case.
    int tag;
    if(!(is >> tag)) return false;
    if(tag != static_cast<int>(IntegratorType::MultiChannel)) return false;

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
    return true;
}

bool achilles::MultiChannel::operator==(const MultiChannel &rhs) const {
    return ndims == rhs.ndims && params == rhs.params && channel_weights == rhs.channel_weights &&
           best_weights == rhs.best_weights && min_diff == rhs.min_diff && summary == rhs.summary;
}

// Integrator interface: forward to the templated cores instantiated on FourVector.
void achilles::MultiChannel::operator()(Integrand<FourVector> &func) {
    this->operator()<FourVector>(func);
}

void achilles::MultiChannel::Optimize(Integrand<FourVector> &func) {
    this->Optimize<FourVector>(func);
}

std::vector<std::vector<achilles::FourVector>>
achilles::MultiChannel::GeneratePoints(Integrand<FourVector> &func, size_t n) {
    std::vector<std::vector<FourVector>> points;
    points.reserve(n);
    for(size_t i = 0; i < n; ++i) points.push_back(GeneratePoint<FourVector>(func));
    return points;
}

std::vector<double>
achilles::MultiChannel::GenerateWeights(Integrand<FourVector> &func,
                                        const std::vector<std::vector<FourVector>> &points) {
    std::vector<double> weights;
    weights.reserve(points.size());
    for(const auto &point : points) weights.push_back(GenerateWeight<FourVector>(func, point));
    return weights;
}
