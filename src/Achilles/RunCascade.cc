#include "Achilles/RunCascade.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Histogram.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Random.hh"
#include "Achilles/Settings.hh"
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdouble-promotion"
#endif
#include "Plugins/HepMC3/HepMC3EventWriter.hh"
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#include "Plugins/NuHepMC/NuHepMCWriter.hh"

#include "spdlog/spdlog.h"

#include <chrono>
#include <stdexcept>

using achilles::CascadeTest::CascadeRunner;

void achilles::CascadeTest::InitCrossSection(Event &event, PID pid, double mom, double radius,
                                             std::shared_ptr<Nucleus> nuc) {
    // Generate a point in the beam of a given radius
    std::array<double, 2> beam_spot;
    auto beam_r = radius * std::sqrt(Random::Instance().Uniform(0.0, 1.0));
    auto beam_theta = Random::Instance().Uniform(0.0, 2 * M_PI);
    beam_spot[0] = beam_r * cos(beam_theta);
    beam_spot[1] = beam_r * sin(beam_theta);

    // Add test particle from the beam coming in on the z-axis
    // The test particle starts outside the nucleus by 5%
    ThreeVector position{beam_spot[0], beam_spot[1], -1.05 * nuc->Radius()};
    auto mass = achilles::ParticleInfo(pid).Mass();
    FourVector momentum{sqrt(mom * mom + mass * mass), 0, 0, mom}; // Check this
    Particle testPart{pid, momentum, position, ParticleStatus::external_test};

    // Make an event instead, initialize with a nucleus, which generates a configuration
    std::vector<FourVector> momdummy = {};
    event = Event(nuc, momdummy, 1.);
    event.Weight() = M_PI * radius * radius * 10 * 1e6; // radius in fm, result in nb

    // Add the test particle
    event.Hadrons().push_back(testPart);
}

void achilles::CascadeTest::InitTransparency(Event &event, PID pid, double mom,
                                             std::shared_ptr<Nucleus> nuc, bool external) {
    double costheta = Random::Instance().Uniform(-1.0, 1.0);
    double sintheta = sqrt(1 - costheta * costheta);
    double phi = Random::Instance().Uniform(0.0, 2 * M_PI);

    std::vector<FourVector> momdummy = {};
    event = Event(nuc, momdummy, 1.);
    event.Weight() = 1.0 / 1000.0; // Only care about percentage, so remove factor of nb_to_pb

    if(!external) {
        // Want to select a random nucleon with PID == pid
        auto protons = event.Protons(ParticleStatus::background);
        auto neutrons = event.Neutrons(ParticleStatus::background);
        Particle *kicked_particle = nullptr;

        if(pid == PID::proton()) {
            kicked_particle = &(Random::Instance().Sample(1, protons))[0].get();
        }

        if(pid == PID::neutron()) {
            kicked_particle = &(Random::Instance().Sample(1, neutrons))[0].get();
        }

        if(kicked_particle == nullptr) {
            throw std::runtime_error("Achilles: For transparency with internal particle, must "
                                     "select particle that exists in ground state of the nucleus");
        }

        auto mass = kicked_particle->Info().Mass();
        FourVector kick{sqrt(mom * mom + mass * mass), mom * sintheta * cos(phi),
                        mom * sintheta * sin(phi), mom * costheta};
        kicked_particle->SetFormationZone(kicked_particle->Momentum(), kick);
        kicked_particle->Status() = ParticleStatus::internal_test;
        kicked_particle->SetMomentum(kick);
    } else {
        // Randomly select a particle from the nucleus
        size_t idx = Random::Instance().Uniform(0ul, event.Hadrons().size() - 1);
        ThreeVector position = event.Hadrons()[idx].Position();

        // Rotate to some other spot on the sphere so that it is not always on top of another
        // nucleon
        double rth = M_PI * Random::Instance().Uniform(0., 1.);
        double rph = 2. * M_PI * Random::Instance().Uniform(0., 1.);
        position = {position.Px() * sin(rth) * cos(rph), position.Py() * sin(rth) * sin(rph),
                    position.Pz() * cos(rth)};

        auto mass = achilles::ParticleInfo(pid).Mass();
        FourVector kick{sqrt(mom * mom + mass * mass), mom * sintheta * cos(phi),
                        mom * sintheta * sin(phi), mom * costheta};

        Particle testPart{pid, kick, position, ParticleStatus::internal_test};
        //    testPart.SetFormationZone(kicked_particle->Momentum(), kick); -> this is implemented
        //    specifically for nucleons!
        event.Hadrons().push_back(testPart);
    }
}

CascadeRunner::CascadeRunner(const std::string &runcard) : setting{runcard, "data/hA_Rules.yml"} {
    // Read momentum range to run over
    m_mom_range = setting.GetAs<std::pair<double, double>>("KickMomentum");
    if(m_mom_range.first < 0.0) {
        throw std::runtime_error("Achilles: Momentum range must be positive");
    }
    if(m_mom_range.first >= m_mom_range.second) {
        throw std::runtime_error("Achilles: Momentum range must be increasing and non-zero");
    }

    // Initialize cascade
    m_mode = setting.GetAs<CascadeMode>("Cascade/Mode");
    spdlog::debug("Cascade mode: {}", ToString(m_mode));

    m_cascade = std::make_unique<Cascade>(setting.GetAs<Cascade>("Cascade"));

    m_params = setting.GetAs<std::map<std::string, double>>("Cascade/Params");
    m_pid = setting.GetAs<PID>("PID");
    requested_events = setting.GetAs<size_t>("NEvents");

    // Initialize nucleus
    m_nuc = std::make_shared<Nucleus>(setting.GetAs<Nucleus>("Nucleus"));
    auto potential_name = setting.GetAs<std::string>("Nucleus/Potential/Name");
    auto potential = achilles::PotentialFactory::Initialize(potential_name, m_nuc,
                                                            setting["Nucleus/Potential"]);
    m_nuc->SetPotential(std::move(potential));

    // Setting radius for hydrogen here
    if(m_nuc->NProtons() == 1 && m_nuc->NNucleons() == 1) { m_nuc->SetRadius(0.84); }
    if(m_nuc->NProtons() == 0 && m_nuc->NNucleons() == 1) { m_nuc->SetRadius(0.84); }

    // Setup output
    auto output = setting["Output"];
    bool zipped = output["Zipped"] ? output["Zipped"].as<bool>() : true;

    spdlog::trace("Outputing as {} format", output["Format"].as<std::string>());
    m_output_name = output["Name"].as<std::string>();
    m_output_name.erase(m_output_name.length() - 6);
    const std::string format = output["Format"].as<std::string>();

    if(format == "Achilles") {
        m_writer = std::make_unique<AchillesWriter>(output["Name"].as<std::string>(), zipped);
    } else if(format == "HepMC3") {
        m_writer = std::make_unique<HepMC3Writer>(output["Name"].as<std::string>(), zipped);
    } else if(format == "NuHepMC") {
        m_writer = std::make_unique<NuHepMCWriter>(output["Name"].as<std::string>(), zipped);
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}",
                                      output["Format"].as<std::string>());
        throw std::runtime_error(msg);
    }

    // Setup dummy process group for writing the header
    ProcessGroup group;
    std::unique_ptr<Unweighter> unweighter = nullptr;
    Process process(ProcessInfo(), std::move(unweighter));
    group.AddProcess(std::move(process));
    std::vector<ProcessGroup> process_groups;
    process_groups.push_back(std::move(group));

    // Write header
    m_writer->WriteHeader(runcard, process_groups);
}

void CascadeRunner::GenerateEvent(double mom, Histogram &hits, Histogram &no_hits) {
    // Generate a point based on the run mode
    Event event;
    switch(m_mode) {
    case CascadeMode::CrossSection:
        InitCrossSection(event, m_pid, mom, m_params["radius"], m_nuc);
        break;
    case CascadeMode::Transparency:
        InitTransparency(event, m_pid, mom, m_nuc, static_cast<bool>(m_params["external"]));
        break;
    }
    // Set kicked indices
    for(size_t i = 0; i < event.Hadrons().size(); ++i) {
        if(event.Hadrons()[i].Status() == ParticleStatus::external_test ||
           event.Hadrons()[i].Status() == ParticleStatus::internal_test) {
            m_cascade->SetKicked(i);
        }
    }
    // Scale event weight by momentum range
    event.Weight() *= (m_mom_range.second - m_mom_range.first);

    // Cascade
    m_cascade->Evolve(event, m_nuc.get());
    // Write the event to file if an interaction happened
    if(event.History().size() > 0) {
        // Set status of the first interaction as the primary interaction
        event.History().Node(0)->Status() = EventHistoryNode::StatusCode::primary;
        // TODO: This is ugly, and the history should store references to the particles
        // However, the first attempt to do this had issues when event.Hadrons() needed to
        // be resized
        event.History().UpdateStatuses(event.Hadrons());
        hits.Fill(mom, event.Weight());
        generated_events++;
    } else {
        no_hits.Fill(mom, event.Weight());
        event.Weight() = 0.0;
    }
    m_writer->Write(event);
}

void achilles::CascadeTest::CascadeRunner::run() {
    fmt::print("Cascade running in {} mode\n", ToString(m_mode));
    fmt::print("  Generating {} events in range [{}, {}] MeV\n", requested_events,
               m_mom_range.first, m_mom_range.second);

    Histogram h_hits(50, m_mom_range.first, m_mom_range.second, "hits");
    Histogram h_nohits(50, m_mom_range.first, m_mom_range.second, "no hits");

    // Generate events
    while(NeedsEvents()) {
        static constexpr size_t statusUpdate = 100;
        if(generated_events % statusUpdate == 0) {
            fmt::print("Generated {} / {} events\r", generated_events, requested_events);
        }
        auto current_mom = Random::Instance().Uniform(m_mom_range.first, m_mom_range.second);
        GenerateEvent(current_mom, h_hits, h_nohits);
    }

    if(m_mode == CascadeMode::Transparency) {
        h_hits.Save(m_output_name + "_hits");
        h_nohits.Save(m_output_name + "_nohits");
    }
}

void achilles::CascadeTest::RunCascade(const std::string &runcard) {
    // TODO: Use settings class to allow included files
    auto config = YAML::LoadFile(runcard);
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["seed"]) seed = config["Initialize"]["seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Initialize CascadeRunner
    CascadeTest::CascadeRunner generator(runcard);
    generator.run();
}
