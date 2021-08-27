#include "nuchic/SpectralFunction.hh"
#include "nuchic/Constants.hh"
#include "nuchic/Utilities.hh"
#include "spdlog/spdlog.h"

#include <fstream>

nuchic::SpectralFunction::SpectralFunction(const std::string &filename) {
    std::ifstream data(filename); 
    std::vector<double> mom, rho_mom;
    std::string line;

    // Remove header
    std::getline(data, line);
    std::getline(data, line);

    // Read each line in
    while(std::getline(data, line)) {
        std::vector<std::string> tokens;
        tokenize(line, tokens);
        mom.emplace_back(std::stod(tokens[0]));
        rho_mom.emplace_back(std::stod(tokens[3]));
    }

    m_spectral = Interp1D(mom, rho_mom); 
    m_spectral.CubicSpline();

    data.close();
}

double nuchic::SpectralFunction::operator()(double mom) const {
    try {
        return m_spectral(mom/(Constant::HBARC/1.0_GeV))/pow(Constant::HBARC/1.0_GeV, 3);
    } catch(std::domain_error &) {
        return 0;
    }
}
