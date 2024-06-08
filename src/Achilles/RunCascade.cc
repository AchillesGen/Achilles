#include "Achilles/RunCascade.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Random.hh"

#include "spdlog/spdlog.h"

#include <chrono>
#include <fstream>
#include <stdexcept>

using achilles::CascadeTest::CalcCrossSection;
using achilles::CascadeTest::CalcCrossSectionMFP;
using achilles::CascadeTest::CalcMeanFreePath;
using achilles::CascadeTest::CalcTransparency;
using achilles::CascadeTest::CalcTransparencyMFP;
using achilles::CascadeTest::CalcTransparencyExternal;

void CalcCrossSection::GenerateEvent(double mom) {
    // Generate a point in the beam of a given radius
    std::array<double, 2> beam_spot;
    while(true) {
        beam_spot[0] = m_radius * Random::Instance().Uniform(-1.0, 1.0);
        beam_spot[1] = m_radius * Random::Instance().Uniform(-1.0, 1.0);
        if(beam_spot[0] * beam_spot[0] + beam_spot[1] * beam_spot[1] < m_radius * m_radius) break;
    }

    // Add test particle from the beam coming in on the z-axis
    // The test particle starts outside the nucleus by 5%
    ThreeVector position{beam_spot[0], beam_spot[1], -1.05 * m_nuc->Radius()};
    auto mass = achilles::ParticleInfo(m_pid).Mass();
    FourVector momentum{ sqrt(mom * mom + mass * mass), 0, 0, mom}; //Check this
    Particle testPart{m_pid, momentum, position, ParticleStatus::external_test};

    
    //Make an event instead, initialize with a nucleus, which generates a configuration
    std::vector<FourVector> momdummy = {};
    Event event(m_nuc.get(), momdummy, 1.);

    //Add the test particle
    event.Hadrons().push_back(testPart);
    
    // Cascade
    m_cascade.SetKicked(event.Hadrons().size() - 1);

    //Use the mean-free path algortihm ? Once it hits it hits
//    m_cascade.Evolve(event, m_nuc.get());
    m_cascade.MeanFreePath(event, m_nuc.get());

    // Analyze output
    spdlog::debug("Final Nucleons:");
    for(const auto &part : event.Hadrons()) { spdlog::debug("  - {}", part); }

    nevents++;
    for(const auto &part : event.Hadrons()) {
        if(part.Status() != ParticleStatus::background &&
	part.Status() != ParticleStatus::external_test) {
            nhits++;
            break;
        }
    }
}

void CalcCrossSection::PrintResults(std::ofstream &out) const {
    double xsec = nhits * M_PI * m_radius * m_radius / nevents * 10;
    double error = sqrt(nhits) * M_PI * m_radius * m_radius / nevents * 10;
    fmt::print("  Calculated xsec: {} +/- {} mb\n", xsec, error);
    out << fmt::format("{},{}\n", xsec, error);
}

void CalcCrossSectionMFP::GenerateEvent(double mom) {

    // Generate a point in the beam of a given radius
    std::array<double, 2> beam_spot;
    while(true) {
        beam_spot[0] = m_radius * Random::Instance().Uniform(-1.0, 1.0);
        beam_spot[1] = m_radius * Random::Instance().Uniform(-1.0, 1.0);
        if(beam_spot[0] * beam_spot[0] + beam_spot[1] * beam_spot[1] < m_radius * m_radius) break;
    }

    // Add test particle from the beam coming in on the z-axis
    // The test particle starts outside the nucleus by 5%
    ThreeVector position{beam_spot[0], beam_spot[1], -1.05 * m_nuc->Radius()};
    auto mass = achilles::ParticleInfo(m_pid).Mass();
    FourVector momentum{0, 0, mom, sqrt(mom * mom + mass * mass)}; //Check this, different than before!
    Particle testPart{m_pid, momentum, position, ParticleStatus::external_test};

    //Todo: fix all the dummy stuff
    std::vector<FourVector> momdummy;
    Event event(m_nuc.get(), momdummy , 1.);

    //Add the test particle
    event.Hadrons().push_back(testPart);
    
    // Cascade
    // Todo : fix the call to NuWro or delete altogheter
    m_cascade.SetKicked(event.Hadrons().size() - 1);
//    m_cascade.NuWro(event, m_nuc.get()); //No implementation for this ? Declared anyway
    m_cascade.Evolve(event, m_nuc.get()); //No implementation for this ? Declared anyway


    // Analyze output
    spdlog::debug("Final Nucleons:");
    for(const auto &part : event.Hadrons()) { spdlog::debug("  - {}", part); }

    nevents++;

    for(const auto &part : event.Hadrons()) {
        if(part.Status() == ParticleStatus::final_state || 
	part.Status() == ParticleStatus::captured) {
            nhits++;
            break;
        }
    }
}

void CalcCrossSectionMFP::PrintResults(std::ofstream &out) const {
    double xsec = nhits * M_PI * m_radius * m_radius / nevents * 10;
    double error = sqrt(nhits) * M_PI * m_radius * m_radius / nevents * 10;
    fmt::print("  Calculated xsec: {} +/- {} mb\n", xsec, error);
    out << fmt::format("{},{}\n", xsec, error);
}

CalcMeanFreePath::CalcMeanFreePath(int pid, std::shared_ptr<Nucleus> nuc, Cascade cascade)
    : RunMode(nuc, std::move(cascade)), m_pid{pid} {
    m_hist = Histogram(100, 0.0, 2 * m_nuc->Radius(), "mfp");
}

void CalcMeanFreePath::GenerateEvent(double kick_mom) {
    double costheta = Random::Instance().Uniform(-1.0, 1.0);
    double sintheta = sqrt(1 - costheta * costheta);
    double phi = Random::Instance().Uniform(0.0, 2 * M_PI);

    ThreeVector position{}; //Always in the center 
    FourVector kick{kick_mom * sintheta * cos(phi), kick_mom * sintheta * sin(phi),
                    kick_mom * costheta, sqrt(kick_mom * kick_mom + Constant::mN * Constant::mN)};
    Particle testPart{m_pid, kick, position, ParticleStatus::internal_test};

    std::vector<FourVector> momdummy;
    Event event(m_nuc.get(),momdummy, 1.);
    //Add the test particle
    event.Hadrons().push_back(testPart);
    // Cascade
    m_cascade.SetKicked(event.Hadrons().size() - 1);

    m_cascade.MeanFreePath(event, m_nuc.get());

    for(const auto &part : event.Hadrons()) {
        if(part.Status() == ParticleStatus::internal_test) {
            m_hist.Fill(part.GetDistanceTraveled());
            break;
        }
    }
}

void CalcMeanFreePath::PrintResults(std::ofstream &out) const {
    m_hist.Save(&out);
    fmt::print("  Histogram saved\n");
}

void CalcTransparency::GenerateEvent(double kick_mom) {
    double costheta = Random::Instance().Uniform(-1.0, 1.0);
    double sintheta = sqrt(1 - costheta * costheta);
    double phi = Random::Instance().Uniform(0.0, 2 * M_PI);
    
    std::vector<FourVector> momdummy; 
    Event event(m_nuc.get(), momdummy, 1.);

    size_t idx = Random::Instance().Uniform(0ul, event.Hadrons().size() - 1);
    m_cascade.SetKicked(idx);

    auto kicked_particle = &(event.Hadrons()[idx]);
    auto mass = kicked_particle->Info().Mass();
    FourVector kick{sqrt(kick_mom * kick_mom + mass * mass), kick_mom * sintheta * cos(phi),
                    kick_mom * sintheta * sin(phi), kick_mom * costheta};
    kicked_particle->SetFormationZone(kicked_particle->Momentum(), kick);
    kicked_particle->Status() = ParticleStatus::internal_test;
    kicked_particle->SetMomentum(kick);

    spdlog::debug("Initial Nucleons:");
    for(const auto &part : event.Hadrons()) { spdlog::debug("  - {}", part); }

    m_cascade.MeanFreePath(event, m_nuc.get());
    nevents++;

    spdlog::debug("Final Nucleons:");
    for(const auto &part : event.Hadrons()) { spdlog::debug("  - {}", part); }

    for(const auto &part : event.Hadrons()) {
        if(part.Status() == ParticleStatus::internal_test) {
            ninteract++;
            distance += part.GetDistanceTraveled();
        } else if(part.Status() == ParticleStatus::captured) {
            ncaptured++;
            ninteract++;
        }
    }
}



void CalcTransparency::PrintResults(std::ofstream &out) const {
    double transparency = 1 - ninteract / nevents;
    double error = sqrt(ninteract / nevents / nevents);
    fmt::print("  Calculated transparency: {} +/- {}\n", transparency, error);
    fmt::print("  Average distance to interact: {}\n", distance / ninteract);
    out << fmt::format("{},{}\n", transparency, error);
}


void CalcTransparencyExternal::GenerateEvent(double kick_mom) {
    double costheta = Random::Instance().Uniform(-1.0, 1.0);
    double sintheta = sqrt(1 - costheta * costheta);
    double phi = Random::Instance().Uniform(0.0, 2 * M_PI);
    
    std::vector<FourVector> momdummy; 
    Event event(m_nuc.get(), momdummy, 1.);

    //We put the external particle at the position of a nucleon -> distribute according to density
    //Note: all this skips the configurations ? 
    size_t idx = Random::Instance().Uniform(0ul, event.Hadrons().size() - 1);
    ThreeVector position = event.Hadrons()[idx].Position();

    //Optional: rotate to some other spot on the sphere so that it is not always on top of another nucleon:
    double rth = M_PI* Random::Instance().Uniform(0.,1.);
    double rph = 2. * M_PI * Random::Instance().Uniform(0.,1.);
    position = { position.Px()*sin(rth)*cos(rph), position.Py()*sin(rth)*sin(rph), position.Pz()*cos(rth) };
    
    auto mass = achilles::ParticleInfo(m_pid).Mass();
    FourVector kick{sqrt(kick_mom * kick_mom + mass * mass), kick_mom * sintheta * cos(phi),
                    kick_mom * sintheta * sin(phi), kick_mom * costheta};

    
    Particle testPart{m_pid, kick, position, ParticleStatus::internal_test};
//    testPart.SetFormationZone(kicked_particle->Momentum(), kick); -> this is implemented specifically for nucleons!
    event.Hadrons().push_back(testPart);
    m_cascade.SetKicked(event.Hadrons().size() - 1);


    spdlog::debug("Initial Hadrons:");
    for(const auto &part : event.Hadrons()) { spdlog::debug("  - {}", part); }

    m_cascade.MeanFreePath(event, m_nuc.get());
    nevents++;

    spdlog::debug("Final Nucleons:");
    for(const auto &part : event.Hadrons()) { spdlog::debug("  - {}", part); }

    for(const auto &part : event.Hadrons()) {
        if(part.Status() == ParticleStatus::interacted) {
            ninteract++;
            distance += part.GetDistanceTraveled();
	    break; //Break to only count 1
        } else if(part.Status() == ParticleStatus::captured) {
            ncaptured++;
            ninteract++;
        }
    }
}

void CalcTransparencyExternal::PrintResults(std::ofstream &out) const {
    double transparency = 1 - ninteract / nevents;
    fmt::print( " interact {}  with totevents {} ", ninteract, nevents);
    double error = sqrt(ninteract / nevents / nevents);
    fmt::print("  Calculated transparency: {} +/- {}\n", transparency, error);
    fmt::print("  Average distance to interact: {}\n", distance / ninteract);
    out << fmt::format("{},{}\n", transparency, error);
}

void CalcTransparencyMFP::GenerateEvent(double kick_mom) {
    double costheta = Random::Instance().Uniform(-1.0, 1.0);
    double sintheta = sqrt(1 - costheta * costheta);
    double phi = Random::Instance().Uniform(0.0, 2 * M_PI);

    
    std::vector<FourVector> momdummy; 
    Event event(m_nuc.get(), momdummy, 1.);

    size_t idx = Random::Instance().Uniform(0ul, event.Hadrons().size() - 1);
    m_cascade.SetKicked(idx);

    auto kicked_particle = &(event.Hadrons()[idx]);
    auto mass = kicked_particle->Info().Mass();
    FourVector kick{sqrt(kick_mom * kick_mom + mass * mass), kick_mom * sintheta * cos(phi),
                    kick_mom * sintheta * sin(phi), kick_mom * costheta};
    kicked_particle->SetFormationZone(kicked_particle->Momentum(), kick);
    kicked_particle->Status() = ParticleStatus::internal_test;
    kicked_particle->SetMomentum(kick);

    spdlog::debug("Initial Nucleons:");
    for(const auto &part : event.Hadrons()) { spdlog::debug("  - {}", part); }

    m_cascade.MeanFreePath_NuWro(event, m_nuc.get(), 10000);
    nevents++;

    spdlog::debug("Final Nucleons:");
    for(const auto &part : event.Hadrons()) { spdlog::debug("  - {}", part); }

    for(const auto &part : event.Hadrons()) {
        if(part.Status() == ParticleStatus::internal_test) {
            ninteract++;
            distance += part.GetDistanceTraveled();
        } else if(part.Status() == ParticleStatus::captured) {
            ncaptured++;
            ninteract++;
        }
    }
}

void CalcTransparencyMFP::PrintResults(std::ofstream &out) const {
    double transparency = 1 - ninteract / nevents;
    double error = sqrt(ninteract / nevents / nevents);
    fmt::print("  Calculated transparency: {} +/- {}\n", transparency, error);
    fmt::print("  Average distance to interact: {}\n", distance / ninteract);
    fmt::print("  Number of captured nucleons: {}\n", ncaptured);
    out << fmt::format("{},{}\n", transparency, error);
}

void achilles::CascadeTest::RunCascade(const std::string &runcard) {
    auto config = YAML::LoadFile(runcard);
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["seed"]) seed = config["Initialize"]["seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load setup
    auto nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());
    auto potential_name = config["Nucleus"]["Potential"]["Name"].as<std::string>();
    auto potential = achilles::PotentialFactory::Initialize(potential_name, nucleus,
                                                            config["Nucleus"]["Potential"]);
    nucleus->SetPotential(std::move(potential));

    auto kick_mom = config["KickMomentum"].as<std::vector<double>>();
    auto nevents = config["NEvents"].as<size_t>();

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Mode"].as<std::string>());
    auto mode = config["Cascade"]["Mode"].as<CascadeMode>();
    auto cascade = config["Cascade"].as<Cascade>();
    std::unique_ptr<RunMode> generator = nullptr;
    switch(mode) {
    case CascadeMode::CrossSection:
        generator = std::make_unique<CalcCrossSection>(config["PID"].as<int>(), nucleus,
                                                       std::move(cascade));
        break;
    case CascadeMode::CrossSectionMFP:
        generator = std::make_unique<CalcCrossSectionMFP>(config["PID"].as<int>(), nucleus,
                                                          std::move(cascade));
        break;
    case CascadeMode::MeanFreePath:
        generator = std::make_unique<CalcMeanFreePath>(config["PID"].as<int>(), nucleus,
                                                       std::move(cascade));
        break;
    case CascadeMode::Transparency:
        generator = std::make_unique<CalcTransparency>(nucleus, std::move(cascade));
        break;
    case CascadeMode::TransparencyMFP:
        generator = std::make_unique<CalcTransparencyMFP>(nucleus, std::move(cascade));
        break;
    case CascadeMode::TransparencyExternal:
        generator = std::make_unique<CalcTransparencyExternal>(config["PID"].as<int>(), nucleus, std::move(cascade));
        break;
    }

    // Open results file
    std::string filename = fmt::format("{}.dat", config["SaveAs"].as<std::string>());
    std::ofstream results(filename);

    //Setting raius for hydrogen here
    if (nucleus->NProtons() == 1 && nucleus->NNucleons() == 1) { nucleus->SetRadius(0.84);}

    // Generate events
    fmt::print("Cascade running in {} mode\n", config["Cascade"]["Mode"].as<std::string>());
    fmt::print("  Generating {} events per momentum point\n", nevents);
    double current_mom = kick_mom[0];
    while(current_mom <= kick_mom[1]) {
        generator->Reset();
        for(size_t i = 0; i < nevents; ++i) {
            generator->GenerateEvent(current_mom);
        }

        fmt::print("  Kick momentum: {} MeV\n", current_mom);
        results << fmt::format("{},", current_mom);
        generator->PrintResults(results);
        current_mom += kick_mom[2];
    }

    results.close();
}
