#include "nuchic/RunCascade.hh"
#include "nuchic/Histogram.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Cascade.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Random.hh"

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <chrono>
#include <fstream>

namespace nuchic {

enum class CascadeMode {
    CrossSection,
    MeanFreePath,
    Transparency,
    CrossSectionMFP,
    TransparencyMFP,
};

}

namespace YAML {

template<>
struct convert<nuchic::CascadeMode> {
    static bool decode(const YAML::Node &node, nuchic::CascadeMode &mode) {
        std::string name = node.as<std::string>();
        if(name == "CrossSection") mode = nuchic::CascadeMode::CrossSection;
        else if(name == "MeanFreePath") mode = nuchic::CascadeMode::MeanFreePath;
        else if(name == "Transparency") mode = nuchic::CascadeMode::Transparency;
        else if(name == "CrossSectionMFP") mode = nuchic::CascadeMode::CrossSectionMFP;
        else if(name == "TransparencyMFP") mode = nuchic::CascadeMode::TransparencyMFP;
        else return false;

        return true;
    }
};

}

namespace nuchic {

class RunMode {
    public:
        RunMode(std::shared_ptr<Nucleus> nuc, Cascade cascade) 
            : m_nuc{nuc}, m_cascade{std::move(cascade)} {}
        virtual ~RunMode() = default;
        virtual void GenerateEvent(double) = 0;
        virtual void PrintResults(std::ofstream&) const = 0;
    protected:
        std::shared_ptr<Nucleus> m_nuc;
        Cascade m_cascade;
};

class CalcCrossSection : public RunMode {
    public:
        CalcCrossSection(int pid, std::shared_ptr<Nucleus> nuc, Cascade cascade, double radius=10)
            : RunMode(nuc, std::move(cascade)), m_radius{std::move(radius)}, m_pid{pid} {}
        void GenerateEvent(double mom) override {
            auto particles = m_nuc -> Nucleons();
           
            // Generate a point in the beam of a given radius
            std::array<double, 2> beam_spot;
            while(true) {
                beam_spot[0] = m_radius*Random::Instance().Uniform(-1.0, 1.0);
                beam_spot[1] = m_radius*Random::Instance().Uniform(-1.0, 1.0);
                if(beam_spot[0]*beam_spot[0] + beam_spot[1]*beam_spot[1] < m_radius*m_radius)
                    break;
            }

            // Add test particle from the beam coming in on the z-axis
            // The test particle starts outside the nucleus by 5%
            ThreeVector position{beam_spot[0], beam_spot[1], -1.05*m_nuc->Radius()};
            auto mass = nuchic::ParticleInfo(m_pid).Mass();
            FourVector momentum{0, 0, mom, sqrt(mom*mom + mass*mass)}; 
            Particle testPart{m_pid, momentum, position, ParticleStatus::external_test};

            // Cascade
            m_cascade.SetKicked(particles.size());
            particles.push_back(testPart);
            m_nuc->SetNucleons(particles);
            m_cascade.Evolve(m_nuc);

            // Analyze output
            spdlog::debug("Final Nucleons:");
            for(const auto &part : m_nuc -> Nucleons()) {
                spdlog::debug("  - {}", part);
            }

            nevents++;
            for(const auto &part : m_nuc->Nucleons()) {
                if(part.Status() == ParticleStatus::escaped) {
                    nhits++; 
                    break;
                }
            }
        }
        void PrintResults(std::ofstream &out) const override {
            double xsec = nhits*M_PI*m_radius*m_radius/nevents*10;
            double error = sqrt(nhits)*M_PI*m_radius*m_radius/nevents*10;
            fmt::print("  Calculated xsec: {} +/- {} mb\n", xsec, error);
            out << fmt::format("{},{}\n",xsec,error);
        }
    private:
        double m_radius;
        int m_pid;
        double nhits{}, nevents{};
};

class CalcCrossSectionMFP : public RunMode {
    public:
        CalcCrossSectionMFP(int pid, std::shared_ptr<Nucleus> nuc, Cascade cascade, double radius=10)
            : RunMode(nuc, std::move(cascade)), m_radius{std::move(radius)}, m_pid{pid} {}
        void GenerateEvent(double mom) override {
            auto particles = m_nuc -> Nucleons();
           
            // Generate a point in the beam of a given radius
            std::array<double, 2> beam_spot;
            while(true) {
                beam_spot[0] = m_radius*Random::Instance().Uniform(-1.0, 1.0);
                beam_spot[1] = m_radius*Random::Instance().Uniform(-1.0, 1.0);
                if(beam_spot[0]*beam_spot[0] + beam_spot[1]*beam_spot[1] < m_radius*m_radius)
                    break;
            }

            // Add test particle from the beam coming in on the z-axis
            // The test particle starts outside the nucleus by 5%
            ThreeVector position{beam_spot[0], beam_spot[1], -1.05*m_nuc->Radius()};
            auto mass = nuchic::ParticleInfo(m_pid).Mass();
            FourVector momentum{0, 0, mom, sqrt(mom*mom + mass*mass)}; 
            Particle testPart{m_pid, momentum, position, ParticleStatus::external_test};

            // Cascade
            m_cascade.SetKicked(particles.size());
            particles.push_back(testPart);
            m_nuc->SetNucleons(particles);
            m_cascade.NuWro(m_nuc);

            // Analyze output
            spdlog::debug("Final Nucleons:");
            for(const auto &part : m_nuc -> Nucleons()) {
                spdlog::debug("  - {}", part);
            }

            nevents++;
            for(const auto &part : m_nuc->Nucleons()) {
                if(part.Status() == ParticleStatus::escaped) {
                    nhits++; 
                    break;
                }
            }
        }
        void PrintResults(std::ofstream &out) const override {
            double xsec = nhits*M_PI*m_radius*m_radius/nevents*10;
            double error = sqrt(nhits)*M_PI*m_radius*m_radius/nevents*10;
            fmt::print("  Calculated xsec: {} +/- {} mb\n", xsec, error);
            out << fmt::format("{},{}\n",xsec,error);
        }
    private:
        double m_radius;
        int m_pid;
        double nhits{}, nevents{};
};

class CalcMeanFreePath : public RunMode {
    public:
        CalcMeanFreePath(int pid, std::shared_ptr<Nucleus> nuc, Cascade cascade) 
            : RunMode(nuc, std::move(cascade)), m_pid{pid} {
            m_hist = Histogram(100, 0.0, 2*m_nuc->Radius(), "mfp");
        }
        void GenerateEvent(double kick_mom) override {
            double costheta = Random::Instance().Uniform(-1.0, 1.0);
            double sintheta = sqrt(1-costheta*costheta);
            double phi = Random::Instance().Uniform(0.0, 2*M_PI);
            auto particles = m_nuc->Nucleons();
            ThreeVector position{};
            FourVector kick{kick_mom*sintheta*cos(phi),
                            kick_mom*sintheta*sin(phi),
                            kick_mom*costheta,
                            sqrt(kick_mom*kick_mom + Constant::mN*Constant::mN)};
            Particle testPart{m_pid, kick, position, ParticleStatus::internal_test};
            m_cascade.SetKicked(particles.size());
            particles.push_back(testPart);
            m_nuc->SetNucleons(particles);
            m_cascade.MeanFreePath(m_nuc);
            for(const auto &part : m_nuc -> Nucleons()) {
                if(part.Status() == ParticleStatus::internal_test) {
                    m_hist.Fill(part.GetDistanceTraveled());
                    break;
                }
            }
        }
        void PrintResults(std::ofstream &out) const override {
            m_hist.Save(&out);
            fmt::print("  Histogram saved\n");
        }

    private:
        int m_pid;
        Histogram m_hist;
};

class CalcTransparency : public RunMode {
    public:
        CalcTransparency(std::shared_ptr<Nucleus> nuc, Cascade cascade) 
            : RunMode(nuc, std::move(cascade)) {}
        void GenerateEvent(double kick_mom) override {
            double costheta = Random::Instance().Uniform(-1.0, 1.0); 
            double sintheta = sqrt(1-costheta*costheta);
            double phi = Random::Instance().Uniform(0.0, 2*M_PI);
            auto particles = m_nuc->Nucleons();
            size_t idx = Random::Instance().Uniform(0ul, particles.size()-1);
            m_cascade.SetKicked(idx);
            auto kicked_particle = &particles[idx];
            auto mass = kicked_particle -> Info().Mass(); 
            FourVector kick{kick_mom*sintheta*cos(phi),
                            kick_mom*sintheta*sin(phi),
                            kick_mom*costheta,
                            sqrt(kick_mom*kick_mom + mass*mass)};
            kicked_particle->SetFormationZone(kicked_particle->Momentum(), kick);
            kicked_particle->Status() = ParticleStatus::internal_test;
            kicked_particle->SetMomentum(kick);
            m_nuc -> SetNucleons(particles);

            spdlog::debug("Initial Nucleons:");
            for(const auto &part : m_nuc -> Nucleons()) {
                spdlog::debug("  - {}", part);
            }

            m_cascade.MeanFreePath(m_nuc);
            nevents++;

            spdlog::debug("Final Nucleons:");
            for(const auto &part : m_nuc -> Nucleons()) {
                spdlog::debug("  - {}", part);
            }

            for(const auto &part : m_nuc -> Nucleons()) {
                if(part.Status() == ParticleStatus::internal_test) {
                    ninteract++;
                    break;
                }
            }
        }
        void PrintResults(std::ofstream &out) const override {
            double transparency = 1-ninteract/nevents;
            double error = sqrt(ninteract/nevents/nevents);
            fmt::print("  Calculated transparency: {} +/- {}\n", transparency, error);
            out << fmt::format("{},{}\n", transparency, error);
        }

    private:
        double nevents{};
        double ninteract{};
};

class CalcTransparencyMFP : public RunMode {
    public:
        CalcTransparencyMFP(std::shared_ptr<Nucleus> nuc, Cascade cascade) 
            : RunMode(nuc, std::move(cascade)) {}
        void GenerateEvent(double kick_mom) override {
            double costheta = Random::Instance().Uniform(-1.0, 1.0); 
            double sintheta = sqrt(1-costheta*costheta);
            double phi = Random::Instance().Uniform(0.0, 2*M_PI);
            auto particles = m_nuc->Nucleons();
            size_t idx = Random::Instance().Uniform(0ul, particles.size()-1);
            m_cascade.SetKicked(idx);
            auto kicked_particle = &particles[idx];
            auto mass = kicked_particle -> Info().Mass(); 
            FourVector kick{kick_mom*sintheta*cos(phi),
                            kick_mom*sintheta*sin(phi),
                            kick_mom*costheta,
                            sqrt(kick_mom*kick_mom + mass*mass)};
            kicked_particle->SetFormationZone(kicked_particle->Momentum(), kick);
            kicked_particle->Status() = ParticleStatus::internal_test;
            kicked_particle->SetMomentum(kick);
            m_nuc -> SetNucleons(particles);

            spdlog::debug("Initial Nucleons:");
            for(const auto &part : m_nuc -> Nucleons()) {
                spdlog::debug("  - {}", part);
            }

            m_cascade.MeanFreePath_NuWro(m_nuc);
            nevents++;

            spdlog::debug("Final Nucleons:");
            for(const auto &part : m_nuc -> Nucleons()) {
                spdlog::debug("  - {}", part);
            }

            for(const auto &part : m_nuc -> Nucleons()) {
                if(part.Status() == ParticleStatus::internal_test) {
                    ninteract++;
                    break;
                }
            }
        }
        void PrintResults(std::ofstream &out) const override {
            double transparency = 1-ninteract/nevents;
            double error = sqrt(ninteract/nevents/nevents);
            fmt::print("  Calculated transparency: {} +/- {}\n", transparency, error);
            out << fmt::format("{},{}\n", transparency, error);
        }

    private:
        double nevents{};
        double ninteract{};
};

}

void nuchic::RunCascade(const std::string &runcard) {
    auto config = YAML::LoadFile(runcard);
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["seed"])
        seed = config["Initialize"]["seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load setup
    auto nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());
    auto kick_mom = config["KickMomentum"].as<std::vector<double>>();
    auto nevents = config["NEvents"].as<size_t>();

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Mode"].as<std::string>());
    auto mode = config["Cascade"]["Mode"].as<CascadeMode>();
    auto cascade = config["Cascade"].as<Cascade>();
    std::unique_ptr<RunMode> generator = nullptr; 
    switch(mode) {
        case CascadeMode::CrossSection:
            generator = std::make_unique<CalcCrossSection>(config["PID"].as<int>(),
                                                           nucleus, std::move(cascade)); 
            break;
        case CascadeMode::CrossSectionMFP:
            generator = std::make_unique<CalcCrossSectionMFP>(config["PID"].as<int>(),
                                                           nucleus, std::move(cascade)); 
            break;
        case CascadeMode::MeanFreePath:
            generator = std::make_unique<CalcMeanFreePath>(config["PID"].as<int>(),
                                                           nucleus, std::move(cascade)); 
            break;
        case CascadeMode::Transparency:
            generator = std::make_unique<CalcTransparency>(nucleus, std::move(cascade)); 
            break;
        case CascadeMode::TransparencyMFP:
            generator = std::make_unique<CalcTransparencyMFP>(nucleus, std::move(cascade)); 
            break;
    }

    // Open results file
    std::string filename = fmt::format("{}.dat", config["SaveAs"].as<std::string>());
    std::ofstream results(filename);

    // Generate events
    fmt::print("Cascade running in {} mode\n", config["Cascade"]["Mode"].as<std::string>());
    fmt::print("  Generating {} events per momentum point\n", nevents);
    double current_mom = kick_mom[0];
    while(current_mom <= kick_mom[1]) {
        for(size_t i = 0; i < nevents; ++i) {
            nucleus -> GenerateConfig();
            generator -> GenerateEvent(current_mom);
        }

        fmt::print("  Kick momentum: {} MeV\n", current_mom);
        results << fmt::format("{},", current_mom);
        generator -> PrintResults(results);
        current_mom += kick_mom[2];
    }

    results.close();
}
