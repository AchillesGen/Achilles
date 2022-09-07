#include "Achilles/Unweighter.hh"
#include "Achilles/Random.hh"

using achilles::PercentileUnweighter;

PercentileUnweighter::PercentileUnweighter(const YAML::Node &node)
    : m_percentile{node["percentile"].as<double>()/100} {}

void PercentileUnweighter::AddEvent(const double &wgt) {
    m_percentile.Add(wgt);
}

bool PercentileUnweighter::AcceptEvent(double &wgt) {
    double max_wgt = m_percentile.Get(); 
    double prob = wgt / max_wgt;
    m_total++;

    if(prob < achilles::Random::Instance().Uniform(0.0, 1.0)) {
        wgt = 0;
        return false;
    }

    m_accepted++;
    wgt = wgt > max_wgt ? wgt : max_wgt;
    return true;
}

std::unique_ptr<achilles::Unweighter> PercentileUnweighter::Construct(const YAML::Node &node) {
    return std::make_unique<PercentileUnweighter>(node);
}
