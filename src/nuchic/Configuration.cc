#include <utility>

#include "nuchic/Configuration.hh"
#include "nuchic/Particle.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/Random.hh"
#include "nuchic/Utilities.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "gzstream/gzstream.h"
#pragma GCC diagnostic pop

nuchic::DensityConfiguration::DensityConfiguration(const std::string &filename) {
    // Load configuration
    igzstream configs(filename.c_str());
    std::string line;
    std::getline(configs, line);
    std::vector<std::string> tokens;
    tokenize(line, tokens);
    m_nnucleons = std::stoul(tokens[0]);
    m_nconfigs = std::stoul(tokens[1]);
    m_maxWgt = std::stod(tokens[2]);
    m_minWgt = std::stod(tokens[3]);

    for(size_t iconfig = 0; iconfig < m_nconfigs; ++iconfig) {
        Configuration config;
        for(size_t inucleon = 0; inucleon < m_nnucleons; ++inucleon) {
            tokens.clear();
            std::getline(configs, line);
            tokenize(line, tokens);
            auto pid = tokens[0] == "1" ? PID::proton() : PID::neutron();
            auto pos = ThreeVector(std::stod(tokens[1]), std::stod(tokens[2]),
                                   std::stod(tokens[3]));
            config.nucleons.emplace_back(pid, FourVector(), pos);
        }
        std::getline(configs, line);
        config.wgt = std::stod(line);
        std::getline(configs, line);
        m_configurations.push_back(config);
    }
}

std::vector<nuchic::Particle> nuchic::DensityConfiguration::GetConfiguration() {
    Configuration config;
    while(true) {
        config = Random::Instance().Pick(m_configurations);

        if(config.wgt/m_maxWgt > Random::Instance().Uniform(0.0, 1.0))
            break;
    }

    std::array<double, 3> angles{};
    Random::Instance().Generate(angles, 0.0, 2*M_PI);
    angles[1] /= 2;

    for(auto& part : config.nucleons) {
        part.SetPosition(part.Position().Rotate(angles));
    }

    return config.nucleons;
}
