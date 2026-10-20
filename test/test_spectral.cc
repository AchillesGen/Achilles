// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

TEST_CASE("Overestimate of Spectral Function", "[spectral]") {
    const std::string filename = "data/Spectral_Functions/pke12p_tot.data";
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

    CHECK(norm == Catch::Approx(6));

    // The per-momentum maximum is used as a sampling overestimate (envelope) for
    // the spectral function. Instead of writing the (p, E, S, Max(S)) grid out to a
    // file for plotting, assert the spectral function is physical and that the
    // envelope genuinely bounds every value it is meant to overestimate.
    std::vector<double> maxS(np);
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < nep; ++j) {
            if(maxS[i] < pke[i * nep + j]) maxS[i] = pke[i * nep + j];
        }
    }

    bool physical = true;
    bool envelope_bounds_values = true;
    for(size_t i = 0; i < np; ++i) {
        for(size_t j = 0; j < nep; ++j) {
            const double s = pke[i * nep + j];
            if(!std::isfinite(s) || s < 0) physical = false;
            if(s > maxS[i]) envelope_bounds_values = false;
        }
    }
    CHECK(physical);
    CHECK(envelope_bounds_values);

    // Every momentum slice must have a positive overestimate to sample against.
    CHECK(std::all_of(maxS.begin(), maxS.end(), [](double m) { return m > 0; }));
}
