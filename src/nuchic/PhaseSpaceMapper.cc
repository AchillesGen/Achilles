#include "nuchic/PhaseSpaceMapper.hh"
#include "nuchic/FourVector.hh"

void nuchic::PSMapper::GeneratePoint(std::vector<FourVector> &momentum, const std::vector<double> &rans) const {
    // Get the dimensions for each component of the PSMapper
    const int hbeamVars = static_cast<int>(hbeam -> NDims());
    const int lbeamVars = static_cast<int>(lbeam -> NDims());
    const int mainVars = static_cast<int>(main -> NDims());

    // Separate the random numbers into the appropriate components
    const std::vector<double> hbeamRans(rans.begin(), rans.begin() + hbeamVars);
    const std::vector<double> lbeamRans(rans.begin() + hbeamVars, rans.begin() + hbeamVars + lbeamVars);
    const std::vector<double> mainRans(rans.end() - mainVars, rans.end());

    // Generate the phase space
    // The momentum are given in the following order:
    // 1. Momentum of the initial hadron
    // 2. Momentum of the initial lepton
    // 3. Momentum of all outgoing parts of the leptonic tensor
    // 4. Momentum of all outgoing hadrons
    momentum.resize(nleptons + nhadrons);
    hbeam -> GeneratePoint(momentum, hbeamRans);
    lbeam -> GeneratePoint(momentum, lbeamRans);
    main -> GeneratePoint(momentum, mainRans);
}

double nuchic::PSMapper::GenerateWeight(const std::vector<FourVector> &momentum,
                                        std::vector<double> &rans) const {
    // Get the dimensions for each component of the PSMapper and create temporary rans vector for each
    const auto hbeamVars = hbeam -> NDims();
    const auto lbeamVars = lbeam -> NDims();
    const auto mainVars = main -> NDims();
    std::vector<double> hbeamRans(hbeamVars), lbeamRans(lbeamVars), mainRans(mainVars);
    
    // Calculate the weights
    double wgt = 1.0;
    wgt /= hbeam -> GenerateWeight(momentum, hbeamRans);
    wgt /= lbeam -> GenerateWeight(momentum, lbeamRans);
    wgt /= main -> GenerateWeight(momentum, mainRans);

    // Merge the random numbers
    hbeamRans.insert(
            hbeamRans.end(),
            std::make_move_iterator(lbeamRans.begin()),
            std::make_move_iterator(lbeamRans.end()));
    hbeamRans.insert(
            hbeamRans.end(),
            std::make_move_iterator(mainRans.begin()),
            std::make_move_iterator(mainRans.end()));
    swap(rans, hbeamRans);

    return 1.0/wgt;
}
