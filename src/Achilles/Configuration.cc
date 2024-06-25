#include <utility>

#include "Achilles/Configuration.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Random.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

#ifdef GZIP
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "gzstream/gzstream.h"
#pragma GCC diagnostic pop
#else
#include <fstream>
#endif

achilles::DensityConfiguration::DensityConfiguration(const std::string &filename) {
    spdlog::debug("Loading density configurations from {}", filename);
    // Load configuration
#ifdef GZIP
    igzstream configs(filename.c_str());
#else
    std::ifstream configs(filename.c_str());
#endif
    if(!configs.good()) throw std::runtime_error("Invalid file " + filename);
    std::string line;
    std::getline(configs, line);
    std::vector<std::string> tokens;
    tokenize(line, tokens);
    m_nnucleons = std::stoul(tokens[0]);
    m_nconfigs = std::stoul(tokens[1]);
    m_maxWgt = std::stod(tokens[2]);
    m_minWgt = std::stod(tokens[3]);

    spdlog::trace("Loading {} configurations with {} nucleons", m_nconfigs, m_nnucleons);

    for(size_t iconfig = 0; iconfig < m_nconfigs; ++iconfig) {
        Configuration config;
        for(size_t inucleon = 0; inucleon < m_nnucleons; ++inucleon) {
            tokens.clear();
            std::getline(configs, line);
            tokenize(line, tokens);
#ifdef ACHILLES_LOW_MEMORY
            auto is_proton = tokens[0] == "1" ? true : false;
            std::array<double, 3> pos{std::stod(tokens[1]), std::stod(tokens[2]),
                                      std::stod(tokens[3])};
            config.nucleons.emplace_back(is_proton, pos);
#else
            const auto pid = tokens[0] == "1" ? PID::proton() : PID::neutron();
            const ThreeVector position{std::stod(tokens[1]), std::stod(tokens[2]),
                                       std::stod(tokens[3])};
            config.nucleons.emplace_back(pid, FourVector{}, position);
#endif
        }
        std::getline(configs, line);
        config.wgt = std::stod(line);
        std::getline(configs, line);
        m_configurations.push_back(config);
    }

    configs.close();
}

std::vector<achilles::Particle> achilles::DensityConfiguration::GetConfiguration() {
#ifdef ACHILLES_LOW_MEMORY
    std::vector<achilles::Particle> particles;
#endif
    while(true) {
        auto index = Random::Instance().Uniform<std::size_t>(0, m_configurations.size());

        if(m_configurations[index].wgt / m_maxWgt > Random::Instance().Uniform(0.0, 1.0)) {
            Configuration &config = m_configurations[index];
            std::array<double, 3> angles{};
            Random::Instance().Generate(angles, 0.0, 2 * M_PI);
            angles[1] /= 2;

#ifdef ACHILLES_LOW_MEMORY
            for(auto &part : config.nucleons) {
                const auto pid = part.is_proton ? PID::proton() : PID::neutron();
                const auto position = ThreeVector(part.position).Rotate(angles);
                particles.emplace_back(pid, FourVector{}, position);
            }

            return particles;
#else
            for(auto &part : config.nucleons) { part.SetPosition(part.Position().Rotate(angles)); }

            return config.nucleons;
#endif
        }
    }
}
