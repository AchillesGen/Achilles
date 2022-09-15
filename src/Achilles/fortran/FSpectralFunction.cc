#include "Achilles/SpectralFunction.hh"

extern "C" {
    achilles::SpectralFunction *LoadSpectralFunction(const char *filename) {
        std::string filename_str(filename); 
        achilles::SpectralFunction *spectral = new achilles::SpectralFunction(filename_str);
        return spectral;
    }

    void DeleteSpectralFunction(achilles::SpectralFunction *spectral) {
        delete spectral;
    }

    double SpectralNormalization(achilles::SpectralFunction *spectral) {
        return spectral -> Normalization();
    }

    double SpectralFunction(achilles::SpectralFunction *spectral, double p, double E) {
        return spectral -> operator()(p, E);
    }
}
