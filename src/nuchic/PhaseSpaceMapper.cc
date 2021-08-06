#include "nuchic/PhaseSpaceMapper.hh"
#include "nuchic/FourVector.hh"
#include "spdlog/spdlog.h"

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
    lbeam -> GeneratePoint(momentum, lbeamRans);
    hbeam -> GeneratePoint(momentum, hbeamRans);
    main -> GeneratePoint(momentum, mainRans);

    // Debugging
    spdlog::trace("PSMapper::GeneratePoint");
    size_t idx = 0;
    for(const auto &mom : momentum) {
        spdlog::trace(" - {}: {}", idx++, mom);
    }
}

double nuchic::PSMapper::GenerateWeight(const std::vector<FourVector> &momentum,
                                        std::vector<double> &rans) const {
    // Get the dimensions for each component of the PSMapper and create temporary rans vector for each
    const auto hbeamVars = hbeam -> NDims();
    const auto lbeamVars = lbeam -> NDims();
    const auto mainVars = main -> NDims();
    std::vector<double> hbeamRans(hbeamVars), lbeamRans(lbeamVars), mainRans(mainVars);
    
    // Calculate the weights
    double lwgt = lbeam -> GenerateWeight(momentum, lbeamRans);
    double hwgt = hbeam -> GenerateWeight(momentum, hbeamRans);
    double mwgt = main -> GenerateWeight(momentum, mainRans);
    double wgt = 1.0/lwgt/hwgt/mwgt;

    spdlog::trace("PSMapper::GenerateWeight:");
    spdlog::trace("  lbeam wgt = {}", lwgt);
    spdlog::trace("  hbeam wgt = {}", hwgt);
    spdlog::trace("  main wgt = {}", mwgt);

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
