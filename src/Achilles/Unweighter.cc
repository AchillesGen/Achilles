#include "Achilles/Unweighter.hh"
#include "Achilles/Random.hh"

using achilles::PercentileUnweighter;

PercentileUnweighter::PercentileUnweighter(const YAML::Node &node)
    : m_percentile{node["percentile"].as<double>() / 100} {}

void PercentileUnweighter::AddEvent(double weight) {
    m_percentile.Add(std::abs(weight));
}

double PercentileUnweighter::AcceptEvent(double weight) {
    double max_wgt = m_percentile.Get();
    double prob = std::abs(weight) / max_wgt;
    m_total++;

    if(prob < achilles::Random::Instance().Uniform(0.0, 1.0)) { return 0; }

    m_accepted++;
    auto sign = weight < 0 ? -1 : 1;
    return prob > 1.0 ? sign * prob : sign * 1.0;
}

std::unique_ptr<achilles::Unweighter> PercentileUnweighter::Construct(const YAML::Node &node) {
    return std::make_unique<PercentileUnweighter>(node);
}
