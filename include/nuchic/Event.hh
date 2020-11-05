#ifndef EVENT_HH
#define EVENT_HH

#include <vector>

namespace nuchic {

class Particle;

using vParticles = std::vector<Particle>;

class Event {
    public:
        Event() = default;
        Event(vParticles&& particles) : m_particles{std::move(particles)} {}
        Event(vParticles&& particles, double wgt) : m_particles{std::move(particles)},
            m_wgt{std::move(wgt)} {}

        void SetWgt(double wgt) { m_wgt = std::move(wgt); }
        void SetParticles(vParticles &&particles) { m_particles = std::move(particles); }
        void AddParticle(const Particle&);

        const vParticles& Particles() const { return m_particles; }
        const double& Weight() const { return m_wgt; }

    private:
        vParticles m_particles{};
        double m_wgt{};
};

}

#endif
