#include "catch2/catch.hpp"
#include <fstream>
#include <iostream>

TEST_CASE("Overestimate of Spectral Function", "[spectral]") {
    const std::string filename = "data/pke12_tot.data";
    std::ifstream data(filename);
    size_t nep{}, np{};
    data >> nep >> np;
    std::vector<double> p(np), pke(nep * np), dp_p(np), xe_p(nep);
    for(size_t i = 0; i < np; ++i) {
        data >> p[i];
        for(size_t j = 0; j < nep; ++j) { data >> xe_p[j] >> pke[i * nep + j]; }
    }
    data.close();

    double norm{};
    double hp = p[1] - p[0];
    double he_p = xe_p[1] - xe_p[0];
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < nep; ++j) {
            dp_p[i] += pke[i * nep + j] * he_p;
            pke[i * nep + j] *= 1e8;
        }
    }
    for(size_t i = 0; i < np; ++i) { norm += p[i] * p[i] * dp_p[i] * 4 * M_PI * hp; }

    CHECK(norm == Approx(6));

    // Find maximums for each p
    std::vector<double> maxS(np);
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < nep; ++j) {
            if(maxS[i] < pke[i * nep + j]) maxS[i] = pke[i * nep + j];
        }
    }

    // Write out data to plot
    std::ofstream result("overestimate.txt");
    result << "p,E,S,Max(S)\n";
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < nep; ++j) {
            result << p[i] << "," << xe_p[j] << "," << pke[i * nep + j] << "," << maxS[i] << "\n";
        }
    }
    result.close();
}
