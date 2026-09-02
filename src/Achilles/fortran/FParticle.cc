// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Particle.hh"

extern "C" {

using namespace achilles;

Particle *CreateParticle(const long int pid, const FourVector *momentum,
                         const ThreeVector *position, const long int istatus) {
    // Positions cross the Fortran boundary in fm.
    auto status = static_cast<ParticleStatus>(istatus);
    const auto p = position->Native();
    return new Particle(pid, *momentum,
                        ThreePosition{p[0].native() * units::fm, p[1].native() * units::fm,
                                      p[2].native() * units::fm},
                        status);
}

Particle *CopyParticle(Particle *self) {
    return new Particle(*self);
}

void DeleteParticle(Particle *self) {
    delete self;
}

long int GetParticleStatus(const Particle *self) {
    return static_cast<long int>(self->Status());
}

ParticleInfo *GetParticleInfo(const Particle *self) {
    return new ParticleInfo(self->Info());
}

FourVector *GetParticleMomentum(const Particle *self) {
    return new FourVector(self->Momentum());
}

ThreeVector *GetParticlePosition(const Particle *self) {
    return new ThreeVector(ThreeVector::FromNative(self->Position().ToArray(units::fm)));
}

Particle *SetParticleMomentum(Particle *self, const FourVector *momentum) {
    self->SetMomentum(*momentum);
    return self;
}

void SetParticlePosition(Particle *self, const ThreeVector *position) {
    const auto p = position->Native();
    self->SetPosition(
        {p[0].native() * units::fm, p[1].native() * units::fm, p[2].native() * units::fm});
}
}
