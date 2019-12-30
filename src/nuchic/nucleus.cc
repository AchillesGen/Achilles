#include "../../include/nuchic/nucleus.hh"

using namespace nuchic;

Nucleus::Nucleus(const int& Z, const int& A, const double& binding_, const double& fermiEnergy_, const ConfigMode& mode) :
    nProtons(Z), nNucleons(A), binding(binding_), fermiEnergy(fermiEnergy_) {

    density = make_density(Z, A, mode) 
}
