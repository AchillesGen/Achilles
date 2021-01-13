#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

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
        Density(const Density&) = default;
        Density(Density&&) = default;
        Density& operator=(const Density&) = default;
        Density& operator=(Density&&) = default;
        virtual ~Density() = default;
        virtual std::vector<Particle> GetConfiguration() = 0;
};

class DensityConfiguration : public Density {
    public:
        DensityConfiguration(const std::string&);
        std::vector<Particle> GetConfiguration() override;

    private:
        size_t m_nconfigs, m_nnucleons;
        double m_maxWgt, m_minWgt;
        std::vector<Configuration> m_configurations;
};

}

#endif
