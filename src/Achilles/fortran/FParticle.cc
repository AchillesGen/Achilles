#include "Achilles/Particle.hh"

extern "C" {

    using namespace achilles; 

    Particle* CreateParticle(const long int pid, const FourVector *momentum,
                             const ThreeVector *position,
                             const long int istatus) {
        auto status = static_cast<ParticleStatus>(istatus);
        return new Particle(pid, *momentum, *position, status);
    }

    Particle* CopyParticle(Particle *self) {
        return new Particle(*self);
    }

    void DeleteParticle(Particle *self) {
        delete self;
    }

    long int GetParticleStatus(const Particle *self) {
        return static_cast<long int>(self -> Status());
    }

    ParticleInfo* GetParticleInfo(const Particle *self) {
        return new ParticleInfo(self -> Info());
    } 

    FourVector* GetParticleMomentum(const Particle *self) {
        return new FourVector(self -> Momentum());
    }

    ThreeVector* GetParticlePosition(const Particle *self) {
        return new ThreeVector(self -> Position());
    }

    Particle* SetParticleMomentum(Particle *self, const FourVector *momentum) {
        self -> SetMomentum(*momentum);
        return self;
    }

    void SetParticlePosition(Particle *self, const ThreeVector *position) {
        self -> SetPosition(*position);
    }

}
