#include "Achilles/Precomputed.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/Random.hh"
#include "Achilles/Utilities.hh"

#ifdef ACHILLES_ENABLE_HEPMC3
#include "plugins/HepMC3/HepMC3EventWriter.hh"
#endif

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <chrono>

// TODO: Add other file formats
achilles::Event achilles::Precomputed::ParseEvent(const std::string &line,
                                                  std::shared_ptr<Nucleus> nuc) {
    std::vector<std::string> tokens;
    tokens.clear();
    tokenize(line, tokens, ",");

    // Setup the initial state of the event
    Event event;
    event.CurrentNucleus() = nuc;
    event.CurrentNucleus()->GenerateConfig();
    event.Weight() = std::stod(tokens[0]);
    FourVector had_mom{std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3]),
                       std::stod(tokens[4])};
    FourVector lep_mom{std::stod(tokens[5]), std::stod(tokens[6]), std::stod(tokens[7]),
                       std::stod(tokens[8])};
    event.Leptons().push_back(Particle(PID::muon(), lep_mom, {}, ParticleStatus::final_state));

    auto neutrons = event.CurrentNucleus()->NeutronsIDs();
    std::vector<size_t> selected;
    Random::Instance().Sample(1, neutrons, selected);
    event.CurrentNucleus()->Nucleons()[selected[0]].Status() = ParticleStatus::initial_state;

    Particle proton(PID::proton());
    proton.Momentum() = had_mom;
    proton.Status() = ParticleStatus::propagating;
    proton.Position() = event.CurrentNucleus()->Nucleons()[selected[0]].Position();
    event.CurrentNucleus()->Nucleons().push_back(proton);

    return event;
}

achilles::Precomputed::RunCascade::RunCascade(const std::string &config_file) {
    auto config = YAML::LoadFile(config_file);

    // Setup random number generator
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["Seed"])
        if(config["Initialize"]["Seed"].as<int>() > 0)
            seed = config["Initialize"]["Seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Setup nucleus
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());
    auto potential_name = config["Nucleus"]["Potential"]["Name"].as<std::string>();
    auto potential = achilles::PotentialFactory::Initialize(potential_name, nucleus,
                                                            config["Nucleus"]["Potential"]);
    nucleus->SetPotential(std::move(potential));

    // Initialize Cascade parameters
    cascade = config["Cascade"].as<Cascade>();

    // Setup outputs
    auto output = config["Main"]["Output"];
    bool zipped = true;
    if(output["Zipped"]) zipped = output["Zipped"].as<bool>();
    spdlog::trace("Outputing as {} format", output["Format"].as<std::string>());
    if(output["Format"].as<std::string>() == "Achilles") {
        writer = std::make_unique<AchillesWriter>(output["Name"].as<std::string>(), zipped);
#ifdef ACHILLES_ENABLE_HEPMC3
    } else if(output["Format"].as<std::string>() == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(output["Name"].as<std::string>(), zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}",
                                      output["Format"].as<std::string>());
        throw std::runtime_error(msg);
    }
    writer->WriteHeader(config_file);
    event_filename = config["Main"]["EventFile"].as<std::string>();
}

void achilles::Precomputed::RunCascade::RunAll() {
    auto run = [&](Event &event) { Run(event); };
    ProcessEventsFromFile(run, event_filename.c_str(), nucleus);
}

void achilles::Precomputed::RunCascade::Run(Event &event) {
    cascade.Evolve(&event);
    std::vector<Particle> init_part, final_part;
    for(const auto &part : event.Particles()) {
        if(part.IsFinal()) final_part.push_back(part);
        if(part.IsInitial()) init_part.push_back(part);
    }
    event.History().AddVertex(init_part[0].Position(), init_part, final_part,
                              EventHistory::StatusCode::primary);
    writer->Write(event);
}
