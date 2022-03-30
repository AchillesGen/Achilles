#include "nuchic/Unweighter.hh"
#include "nuchic/Event.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Random.hh"

using nuchic::PercentileUnweighter;

PercentileUnweighter::PercentileUnweighter(const YAML::Node &node)
    : m_percentile{node["percentile"].as<double>()/100} {}

void PercentileUnweighter::AddEvent(const nuchic::Event &event) {
    m_percentile.Add(event.Weight());
}

bool PercentileUnweighter::AcceptEvent(nuchic::Event &event) {
    double max_wgt = m_percentile.Get(); 
    double prob = event.Weight() / max_wgt;
    m_total++;

    if(prob < nuchic::Random::Instance().Uniform(0.0, 1.0)) {
        event.Weight() = 0;
        return false;
    }

    m_accepted++;
    event.Weight() = event.Weight() > max_wgt ? event.Weight() : max_wgt;
    return true;
}

std::unique_ptr<nuchic::Unweighter> PercentileUnweighter::Construct(const YAML::Node &node) {
    return std::make_unique<PercentileUnweighter>(node);
}
