#include "Achilles/Precomputed.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/Process.hh"
#include "Achilles/Random.hh"
#include "Achilles/Settings.hh"
#include "Achilles/Utilities.hh"
#include "Plugins/NuHepMC/NuHepMCWriter.hh"
#include <functional>
#include <sstream>

#ifdef ACHILLES_ENABLE_HEPMC3
#include "Plugins/HepMC3/HepMC3EventWriter.hh"
#endif

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <chrono>

using namespace achilles::Precomputed;

bool iequals(const std::string &a, const std::string &b) {
    return a.size() == b.size() &&
           std::equal(a.begin(), a.end(), b.begin(), [](unsigned char c1, unsigned char c2) {
               return std::tolower(c1) == std::tolower(c2);
           });
}

RawEventReader::RawEventReader(const std::string &filename, std::shared_ptr<Nucleus> nuc)
    : event_stream(filename), m_nucleus(nuc) {
    if(!event_stream.is_open()) {
        auto err_msg = fmt::format("Precomputed: File {} does not exist!", filename);
        throw std::runtime_error(err_msg);
    }

    // Read the preamble
    // NOTE: Currently just the units exist in preamble
    std::string line;
    std::getline(event_stream, line);
    if(iequals("gev", line)) to_mev = 1000;
}

std::optional<achilles::Event> RawEventReader::Next() {
    // Read the line for the event
    std::string line;
    if(!std::getline(event_stream, line)) return std::nullopt;

    return ParseLine(line);
}

achilles::Event RawEventReader::ParseLine(const std::string &line) {
    std::istringstream ss(line);

    // Read in event weight
    double wgt;
    ss >> wgt;

    // Read in all the particles
    std::vector<Particle> particles;
    while(ss) {
        long int pdg;
        double E, px, py, pz;
        ss >> pdg >> E >> px >> py >> pz;
        FourVector mom{E * to_mev, px * to_mev, py * to_mev, pz * to_mev};
        particles.emplace_back(pdg, mom);
    }

    // Setup event
    Event event(m_nucleus, {}, 1);
    event.Weight() = wgt;
    ConvertEvent(event, particles);

    return event;
}
// TODO: Add other file formats?

void achilles::Precomputed::ConvertEvent(Event &event, Particles &particles) {
    // Set the status of the particles and add to the event
    bool first_lepton = true;
    bool first_hadron = true;
    ThreeVector initial_pos;
    for(auto &particle : particles) {
        if(particle.Info().IsLepton()) {
            // Set the first lepton to the incoming one and final state otherwise
            if(first_lepton) {
                particle.Status() = ParticleStatus::initial_state;
                first_lepton = false;
            } else {
                particle.Status() = ParticleStatus::final_state;
            }
            event.Leptons().push_back(particle);
        }

        if(particle.Info().IsHadron()) {
            // Set the first hadron to the one within the nucleus and the rest to propagating
            if(first_hadron) {
                first_hadron = false;

                if(particle.ID() == PID::proton()) {
                    auto protons = event.Protons(ParticleStatus::background);
                    auto sampled_protons = Random::Instance().Sample(1, protons);
                    auto sampled_proton = sampled_protons[0];
                    sampled_proton.get().Momentum() = particle.Momentum();
                    sampled_proton.get().Status() = ParticleStatus::initial_state;
                    initial_pos = sampled_proton.get().Position();
                }

                if(particle.ID() == PID::neutron()) {
                    auto neutrons = event.Neutrons(ParticleStatus::background);
                    auto sampled_neutrons = Random::Instance().Sample(1, neutrons);
                    auto sampled_neutron = sampled_neutrons[0];
                    sampled_neutron.get().Momentum() = particle.Momentum();
                    sampled_neutron.get().Status() = ParticleStatus::initial_state;
                    initial_pos = sampled_neutron.get().Position();
                }
            } else {
                particle.Status() = ParticleStatus::propagating;
                particle.Position() = initial_pos;
                event.Hadrons().push_back(particle);
            }
        }
    }

    // Setup the event history
    auto init_nuc = event.CurrentNucleus()->InitParticle();
    std::vector<Particle> init_parts;
    for(const auto &nucleon : event.Hadrons()) {
        if(nucleon.Status() == ParticleStatus::initial_state) { init_parts.push_back(nucleon); }
    }
    event.History().AddVertex(init_parts[0].Position(), {init_nuc}, init_parts,
                              EventHistory::StatusCode::target);
    // Setup beam in history
    auto init_lep = event.Leptons()[0];
    auto init_beam = init_lep;
    init_beam.Status() = ParticleStatus::beam;
    event.History().AddVertex({}, {init_beam}, {init_lep}, EventHistory::StatusCode::beam);

    // TODO: Figure out how to best handle tracking this with the cascade and decays
    std::vector<Particle> primary_out, propagating;
    for(const auto &part : event.Particles()) {
        if(part.IsFinal()) primary_out.push_back(part);
        if(part.IsPropagating()) {
            primary_out.push_back(part);
            propagating.push_back(part);
        }
    }
    init_parts.push_back(init_lep);
    event.History().AddVertex(init_parts[0].Position(), init_parts, primary_out,
                              EventHistory::StatusCode::primary);
}

RunCascade::RunCascade(const std::string &config_file) {
    auto config = Settings(config_file, "data/Rules_precomputed.yml");

    // Setup random number generator
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config.Exists("Options/Initialize/Seed"))
        if(config.GetAs<int>("Options/Initialize/Seed") > 0)
            seed = config.GetAs<unsigned int>("Options/Initialize/Seed");
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Setup nuclei (currently only supports a single nucleus)
    auto nuclei = config.GetAs<std::vector<std::shared_ptr<Nucleus>>>("Nuclei");
    if(nuclei.size() == 0)
        throw std::runtime_error("Precomputed: Must have at least one nucleus defined");
    nucleus = nuclei[0];

    // Initialize Cascade parameters
    cascade = config.GetAs<Cascade>("Cascade");

    // Setup outputs
    bool zipped =
        config.Exists("Main/Output/Zipped") ? config.GetAs<bool>("Main/Output/Zipped") : true;
    auto format = config.GetAs<std::string>("Main/Output/Format");
    auto name = config.GetAs<std::string>("Main/Output/Name");
    spdlog::trace("Outputing as {} format", format);
    if(format == "Achilles") {
        writer = std::make_unique<AchillesWriter>(name, zipped);
#ifdef ACHILLES_ENABLE_HEPMC3
    } else if(format == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(name, zipped);
    } else if(format == "NuHepMC") {
        writer = std::make_unique<NuHepMCWriter>(name, zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}", format);
        throw std::runtime_error(msg);
    }
    std::vector<ProcessGroup> dummy;
    writer->WriteHeader(config_file, dummy);
    event_filename = config["Main"]["EventFile"].as<std::string>();
}

void RunCascade::RunAll() {
    auto run = [&](Event &event) { Run(event); };
    EventPipeline pipeline(run);
    RawEventReader reader(event_filename, nucleus);
    RunPipeline(reader, pipeline);
}

void RunCascade::Run(Event &event) {
    cascade.Evolve(event, nucleus.get());
    std::vector<Particle> init_part, final_part;
    for(const auto &part : event.Particles()) {
        if(part.IsFinal()) final_part.push_back(part);
        if(part.IsInitial()) init_part.push_back(part);
    }
    event.History().UpdateStatuses(event.Hadrons());
    writer->Write(event);
}
