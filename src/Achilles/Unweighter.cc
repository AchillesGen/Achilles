#include "Achilles/Unweighter.hh"
#include "Achilles/Random.hh"

using achilles::PercentileUnweighter;

PercentileUnweighter::PercentileUnweighter(const YAML::Node &node)
    : m_percentile{node["percentile"].as<double>() / 100} {}

void PercentileUnweighter::AddEvent(double weight) {
    m_percentile.Add(weight);
}

double PercentileUnweighter::AcceptEvent(double weight) {
    double max_wgt = m_percentile.Get();
    double prob = weight / max_wgt;
    m_total++;

    if(prob < achilles::Random::Instance().Uniform(0.0, 1.0)) { return 0; }

    m_accepted++;
    return weight > max_wgt ? weight / max_wgt : 1.0;
}

std::unique_ptr<achilles::Unweighter> PercentileUnweighter::Construct(const YAML::Node &node) {
    return std::make_unique<PercentileUnweighter>(node);
}
