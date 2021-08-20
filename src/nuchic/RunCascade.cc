#include "nuchic/RunCascade.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Cascade.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Random.hh"

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <chrono>

namespace nuchic {

enum class CascadeMode {
    CrossSection,
    MeanFreePath,
    Transparency,
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
        virtual void PrintResults() = 0;
    protected:
        std::shared_ptr<Nucleus> m_nuc;
        Cascade m_cascade;
};

class CalcCrossSection : public RunMode {
    public:
        CalcCrossSection(std::shared_ptr<Nucleus> nuc, Cascade cascade, double radius=10)
            : RunMode(nuc, std::move(cascade)), m_radius{std::move(radius)} {}
        void GenerateEvent(double) override {
        
        }
        void PrintResults() override {}
    private:
        double m_radius;
};

class CalcMeanFreePath : public RunMode {
    public:
        CalcMeanFreePath(std::shared_ptr<Nucleus> nuc, Cascade cascade) 
            : RunMode(nuc, std::move(cascade)) {}
        void GenerateEvent(double) override {

        }
        void PrintResults() override {}
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
            FourVector kick{kick_mom*sintheta*cos(phi),
                            kick_mom*sintheta*sin(phi),
                            kick_mom*costheta,
                            sqrt(kick_mom*kick_mom + Constant::mN*Constant::mN)};
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
        void PrintResults() override {
            double transparency = 1-ninteract/nevents;
            double error = sqrt(ninteract/nevents/nevents);
            fmt::print("  Calculated transparency: {} +/- {}\n", transparency, error);
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
    auto kick_mom = config["KickMomentum"].as<double>();
    auto nevents = config["NEvents"].as<size_t>();

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Mode"].as<std::string>());
    auto mode = config["Cascade"]["Mode"].as<CascadeMode>();
    auto cascade = config["Cascade"].as<Cascade>();
    std::unique_ptr<RunMode> generator = nullptr; 
    switch(mode) {
        case CascadeMode::CrossSection:
            generator = std::make_unique<CalcCrossSection>(nucleus, std::move(cascade)); 
            break;
        case CascadeMode::MeanFreePath:
            generator = std::make_unique<CalcMeanFreePath>(nucleus, std::move(cascade)); 
            break;
        case CascadeMode::Transparency:
            generator = std::make_unique<CalcTransparency>(nucleus, std::move(cascade)); 
            break;
    }

    // Generate events
    for(size_t i = 0; i < nevents; ++i) {
        nucleus -> GenerateConfig();
        generator -> GenerateEvent(kick_mom);
    }

    fmt::print("Cascade ran in {} mode\n", config["Cascade"]["Mode"].as<std::string>());
    fmt::print("  Kick momentum: {} MeV\n", kick_mom);
    fmt::print("  Generated {} events\n", nevents);
    generator -> PrintResults();
}
