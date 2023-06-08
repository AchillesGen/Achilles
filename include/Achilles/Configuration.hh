#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <array>
#include <string>
#include <vector>

namespace achilles {

class Particle;

#ifdef ACHILLES_LOW_MEMORY
struct ConfigurationParticle {
    ConfigurationParticle(bool proton, std::array<double, 3> _position)
        : is_proton(proton), position{std::move(_position)} {}
    bool is_proton;
    std::array<double, 3> position;
};
#endif

struct Configuration {
    std::vector<ConfigurationParticle> nucleons;
    double wgt;
};
#else
struct Configuration {
    std::vector<Particle> nucleons;
    double wgt;
};
#endif

class Density {
  public:
    Density() = default;
    Density(const Density &) = default;
    Density(Density &&) = default;
    Density &operator=(const Density &) = default;
    Density &operator=(Density &&) = default;
    virtual ~Density() = default;
    virtual std::vector<Particle> GetConfiguration() = 0;
};

class DensityConfiguration : public Density {
  public:
    DensityConfiguration(std::string);
    std::vector<Particle> GetConfiguration() override;

  private:
    size_t m_nconfigs, m_nnucleons;
    double m_maxWgt, m_minWgt;
    std::vector<Configuration> m_configurations;
};

} // namespace achilles

#endif
