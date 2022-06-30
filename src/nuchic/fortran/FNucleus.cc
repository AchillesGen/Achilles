#include "nuchic/Nucleus.hh"

extern "C" {
    void GetNucleon(nuchic::Nucleus *nuc, size_t i, nuchic::Particle *part) {
        part = &(nuc -> Nucleons()[i]);
    }

    void AddParticle(nuchic::Nucleus *nuc, nuchic::Particle *part) {
        nuchic::Particle tmp(*part);
        nuc -> Nucleons().push_back(tmp);
    }
}
