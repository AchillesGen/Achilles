#include "nuchic/PhaseSpaceMapper.hh"
#include "nuchic/FourVector.hh"

void nuchic::PSMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) const {
    // Get the dimensions for each component of the PSMapper
    const int initialVars = static_cast<int>(initial -> NDims());
    const int leptonVars = static_cast<int>(leptonic -> NDims());

    // Separate the random numbers into the appropriate components
    const std::vector<double> initialRans(rans.begin(), rans.begin() + initialVars);
    const std::vector<double> leptonRans(rans.begin() + initialVars, rans.begin() + initialVars + leptonVars);
    const std::vector<double> hadronRans(rans.begin() + initialVars + leptonVars, rans.end());

    // Generate the phase space
    // The momentum are given in the following order:
    // 1. Momentum of the probe
    // 2. Momentum of all the leptons (with incoming lepton first)
    // 3. Momentum of all the hadrons (with incoming hadrons first)
    std::vector<FourVector> momentum(nleptons+nhadrons+1);
    initial -> GeneratePoint(momentum, initialRans);
    leptonic -> GeneratePoint(momentum, leptonRans);
    hadronic -> GeneratePoint(momentum, hadronRans);
    swap(point, momentum);
}

double nuchic::PSMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) const {
    // Get the dimensions for each component of the PSMapper and create temporary rans vector for each
    const size_t initialVars = initial -> NDims();
    const size_t leptonVars = leptonic -> NDims();
    const size_t hadronVars = hadronic -> NDims();
    std::vector<double> initialRans(initialVars), leptonRans(leptonVars), hadronRans(hadronVars);

    // Separate the momentum into the appropriate components
    const std::vector<FourVector> initialMom(point.begin(), point.begin() + 2);
    const std::vector<FourVector> leptonMom(point.begin(), point.begin() + 1 + static_cast<int>(nleptons));
    std::vector<FourVector> hadronMom(point.end() - static_cast<int>(nhadrons), point.end());
    hadronMom.push_back(point.front());
    
    // Calculate the weights
    double wgt = 1.0;
    wgt = initial -> GenerateWeight(initialMom, initialRans);
    wgt = leptonic -> GenerateWeight(leptonMom, leptonRans);
    wgt = hadronic -> GenerateWeight(hadronMom, hadronRans);

    // Merge the random numbers
    initialRans.insert(
            initialRans.end(),
            std::make_move_iterator(leptonRans.begin()),
            std::make_move_iterator(leptonRans.end()));
    initialRans.insert(
            initialRans.end(),
            std::make_move_iterator(hadronRans.begin()),
            std::make_move_iterator(hadronRans.end()));
    swap(rans, initialRans);

    return wgt;
}
