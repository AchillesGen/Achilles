#include "Achilles/Unweighter.hh"
#include "Achilles/Event.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Random.hh"

using achilles::PercentileUnweighter;

PercentileUnweighter::PercentileUnweighter(const YAML::Node &node)
    : m_percentile{node["percentile"].as<double>() / 100} {}

void PercentileUnweighter::AddEvent(const achilles::Event &event) {
    m_percentile.Add(event.Weight());
}

bool PercentileUnweighter::AcceptEvent(achilles::Event &event) {
    double max_wgt = m_percentile.Get();
    double prob = event.Weight() / max_wgt;
    m_total++;

    if(prob < achilles::Random::Instance().Uniform(0.0, 1.0)) {
        event.Weight() = 0;
        return false;
    }

    m_accepted++;
    event.Weight() = event.Weight() > max_wgt ? event.Weight() : max_wgt;
    return true;
}

std::unique_ptr<achilles::Unweighter> PercentileUnweighter::Construct(const YAML::Node &node) {
    return std::make_unique<PercentileUnweighter>(node);
}
