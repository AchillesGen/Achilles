#include "Achilles/Beams.hh"
#include "Achilles/FluxLoaders.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Utilities.hh"

#include <algorithm>
#include <limits>

using achilles::FlatFlux;
using achilles::PDFBeam;
using achilles::Spectrum;

Spectrum::Spectrum(const FluxData &data) {
    Interp1D interp(data.bin_centers, data.heights, InterpolationType::Polynomial);
    interp.SetPolyOrder(1);

    m_flux = [interp](double x) { return interp(x); };
    m_flux_integral = data.integral;
    m_min_energy = data.min_energy;
    m_max_energy = data.max_energy;
    m_delta_energy = m_max_energy - m_min_energy;
    m_energy_units = data.energy_units;
    m_format_str = data.format;
}

std::string Spectrum::Format() const {
    return m_format_str;
}

achilles::FourVector Spectrum::Flux(const std::vector<double> &ran, double min_energy) const {
    min_energy = std::max(min_energy * m_energy_units, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    double energy = (ran[0] * delta_energy + min_energy) / m_energy_units;
    return {energy, 0, 0, energy};
}

double Spectrum::GenerateWeight(const FourVector &beam, std::vector<double> &ran,
                                double min_energy) const {
    min_energy = std::max(min_energy * m_energy_units, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    ran[0] = (beam.E() * m_energy_units - min_energy) / delta_energy;
    return (delta_energy * m_flux(beam.E() * m_energy_units)) / m_flux_integral;
}

double Spectrum::EvaluateFlux(const FourVector &beam) const {
    return m_flux(beam.E() * m_energy_units);
}

PDFBeam::PDFBeam(const YAML::Node &) {
    throw std::runtime_error("PDFBeam: Not implemented yet!");
}

achilles::FourVector PDFBeam::Flux(const std::vector<double> &, double) const {
    throw std::runtime_error("PDFBeam: Not implemented yet!");
    return {};
}

double PDFBeam::GenerateWeight(const FourVector &, std::vector<double> &, double) const {
    throw std::runtime_error("PDFBeam: Not implemented yet!");
    return {};
}

double PDFBeam::EvaluateFlux(const FourVector &) const {
    throw std::runtime_error("PDFBeam: Not implemented yet!");
    return {};
}

FlatFlux::FlatFlux(const YAML::Node &node) {
    m_min_energy = node["MinEnergy"].as<double>();
    m_max_energy = node["MaxEnergy"].as<double>();
    m_flux_integral = m_max_energy - m_min_energy;
}

achilles::FourVector FlatFlux::Flux(const std::vector<double> &ran, double min_energy) const {
    min_energy = std::max(min_energy, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    double energy = ran[0] * delta_energy + min_energy;
    return {energy, 0, 0, energy};
}

double FlatFlux::GenerateWeight(const FourVector &beam, std::vector<double> &ran,
                                double min_energy) const {
    min_energy = std::max(min_energy, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    ran[0] = (beam.E() - min_energy) / delta_energy;
    return delta_energy / m_flux_integral;
}

achilles::GeometryBeam::GeometryBeam(const YAML::Node &node) {
    m_min_energy = node["EnergyRange"][0].as<double>();
    m_max_energy = node["EnergyRange"][1].as<double>();
    m_flux_integral = m_max_energy - m_min_energy;
}

void achilles::GeometryBeam::SetInjected(const FourVector &lab_ray) {
    static const FourVector zaxis{0, 0, 0, 1};
    // Rotation taking the lab ray's direction onto +z (RotateBack undoes it).
    m_rotation = Poincare(lab_ray, zaxis);
    m_injected = m_rotation * lab_ray;
}

achilles::FourVector achilles::GeometryBeam::Flux(const std::vector<double> &ran,
                                                  double min_energy) const {
    // Generate mode: serve the injected ray (already rotated onto +z).
    if(m_mode == Mode::Generate) return m_injected;

    // Warm-up mode: flat flux over the configured range (like FlatFlux).
    min_energy = std::max(min_energy, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    double energy = ran[0] * delta_energy + min_energy;
    return {energy, 0, 0, energy};
}

double achilles::GeometryBeam::GenerateWeight(const FourVector &beam, std::vector<double> &ran,
                                              double min_energy) const {
    min_energy = std::max(min_energy, m_min_energy);
    double delta_energy = m_max_energy - min_energy;
    ran[0] = (beam.E() - min_energy) / delta_energy;

    // Generate mode: the energy is the injected ray's, not one mapped from
    // ran[0], so it can sit outside [min_energy, m_max_energy] (flux file
    // wider than EnergyRange, or below a process threshold).  The remapped
    // ran[0] then leaves [0, 1] and downstream VEGAS grid lookups read out
    // of bounds, producing garbage (even negative) jacobians.  Clamp the
    // random number so the grid lookup stays in range and return an infinite
    // phase-space volume: the mapper inverts this to a zero density, the
    // event weight vanishes, and the trial is cleanly rejected.
    if(ran[0] < 0.0 || ran[0] > 1.0) {
        spdlog::debug("GeometryBeam: energy {} MeV outside supported range [{}, {}]; "
                      "rejecting trial",
                      beam.E(), min_energy, m_max_energy);
        ran[0] = std::clamp(ran[0], 0.0, 1.0);
        return std::numeric_limits<double>::infinity();
    }
    spdlog::trace("GeometryBeam: E = {} MeV -> ran[0] = {}, phase-space volume = {}", beam.E(),
                  ran[0], delta_energy);
    return delta_energy / m_flux_integral;
}

achilles::Beam::Beam(BeamMap beams) : m_beams{std::move(beams)} {
    n_vars = 0;
    for(const auto &beam : m_beams) {
        if(beam.second->NVariables() > n_vars) n_vars = beam.second->NVariables();
        spdlog::debug("Beam with PID: {} created.", beam.first);
        m_pids.insert(beam.first);
        // Every species must share the same energy support (single EnergyDomain).
        if(m_beams.begin()->second->MinEnergy() != beam.second->MinEnergy() ||
           m_beams.begin()->second->MaxEnergy() != beam.second->MaxEnergy())
            throw std::logic_error(fmt::format(
                "Beam: all species must share the same energy support. Got [{}, {}] and [{}, {}]",
                m_beams.begin()->second->MinEnergy(), m_beams.begin()->second->MaxEnergy(),
                beam.second->MinEnergy(), beam.second->MaxEnergy()));
    }
    if(!m_beams.empty())
        m_domain = EnergyDomain(m_beams.begin()->second->MinEnergy(),
                                m_beams.begin()->second->MaxEnergy());
}

double achilles::Beam::MaxEnergy() const {
    double max = 0;
    for(const auto &beam : m_beams) {
        if(beam.second->MaxEnergy() > max) max = beam.second->MaxEnergy();
    }
    return max;
}

double achilles::Beam::MinEnergy() const {
    double min = std::numeric_limits<double>::max();
    for(const auto &beam : m_beams) {
        if(beam.second->MinEnergy() < min) min = beam.second->MinEnergy();
    }
    return min;
}
