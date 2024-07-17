#include "Achilles/SpectralFunction.hh"

extern "C" {
achilles::SpectralFunction *LoadSpectralFunction(char *filename) {
    std::string filename_str(filename);
    auto spectral = new achilles::SpectralFunction(filename_str);
    return spectral;
}

void DeleteSpectralFunction(achilles::SpectralFunction *spectral) {
    delete spectral;
}

double SpectralNormalization(achilles::SpectralFunction *spectral) {
    return spectral->Normalization();
}

double SpectralFunction(achilles::SpectralFunction *spectral, double p, double E) {
    return spectral->operator()(p, E);
}

double MomentumDistribution(achilles::SpectralFunction *spectral, double p) {
    return spectral->operator()(p);
}
}
