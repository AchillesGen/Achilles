#include "Achilles/SpectralFunction.hh"

extern "C" {
    achilles::SpectralFunction* LoadSpectralFunction(char *filename, size_t len) {
        std::string filename_str(filename, filename+len); 
        return new achilles::SpectralFunction(filename_str);
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
