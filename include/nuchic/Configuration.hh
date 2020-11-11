#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include "nuchic/Random.hh"
#include <string>
#include <vector>

namespace nuchic {

class Particle;

struct Configuration {
    std::vector<Particle> nucleons; 
    double wgt;
};

class Density {
    public:
        Density() = default;
        virtual ~Density() = default;
        virtual std::vector<Particle> GetConfiguration() = 0;
};

class DensityConfiguration : public Density {
    public:
        DensityConfiguration(const std::string&, std::shared_ptr<randutils::mt19937_rng>);
        std::vector<Particle> GetConfiguration() override;

    private:
        size_t m_nconfigs, m_nnucleons;
        double m_maxWgt, m_minWgt;
        std::vector<Configuration> m_configurations;
        std::shared_ptr<randutils::mt19937_rng> m_rng;
};

}

#endif
