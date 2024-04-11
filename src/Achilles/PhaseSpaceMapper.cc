#include "Achilles/PhaseSpaceMapper.hh"
#include "Achilles/FourVector.hh"
#include "spdlog/spdlog.h"

void achilles::PSMapper::GeneratePoint(std::vector<FourVector> &momentum,
                                       const std::vector<double> &rans) {
    // Get the dimensions for each component of the PSMapper
    const int hbeamVars = static_cast<int>(hbeam->NDims());
    const int lbeamVars = static_cast<int>(lbeam->NDims());
    const int mainVars = static_cast<int>(main->NDims());

    // Separate the random numbers into the appropriate components
    const std::vector<double> hbeamRans(rans.begin(), rans.begin() + hbeamVars);
    const std::vector<double> lbeamRans(rans.begin() + hbeamVars,
                                        rans.begin() + hbeamVars + lbeamVars);
    const std::vector<double> mainRans(rans.end() - mainVars, rans.end());

    // Generate the phase space
    // The momentum are given in the following order:
    // 1. Momentum of the initial lepton
    // 2. Momentum of the initial hadron
    // 3. Momentum of all outgoing parts of the leptonic tensor
    // 4. Momentum of all outgoing hadrons
    // 5. Momentum of all spectators
    momentum.resize(nleptons + nhadrons + nspectators);
    lbeam->GeneratePoint(momentum, lbeamRans);
    hbeam->GeneratePoint(momentum, hbeamRans);
    main->GeneratePoint(momentum, mainRans);

    // Debugging
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, momentum, rans);
}

double achilles::PSMapper::GenerateWeight(const std::vector<FourVector> &momentum,
                                          std::vector<double> &rans) {
    // Get the dimensions for each component of the PSMapper and create
    // temporary rans vector for each
    const auto hbeamVars = hbeam->NDims();
    const auto lbeamVars = lbeam->NDims();
    const auto mainVars = main->NDims();
    std::vector<double> hbeamRans(hbeamVars), lbeamRans(lbeamVars), mainRans(mainVars);

    // Calculate the weights
    double lwgt = lbeam->GenerateWeight(momentum, lbeamRans);
    double hwgt = hbeam->GenerateWeight(momentum, hbeamRans);
    double mwgt = main->GenerateWeight(momentum, mainRans);
    double wgt = lwgt * hwgt * mwgt;

    // Merge the random numbers
    hbeamRans.insert(hbeamRans.end(), std::make_move_iterator(lbeamRans.begin()),
                     std::make_move_iterator(lbeamRans.end()));
    hbeamRans.insert(hbeamRans.end(), std::make_move_iterator(mainRans.begin()),
                     std::make_move_iterator(mainRans.end()));
    swap(rans, hbeamRans);

    // Debugging
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, momentum, rans);
    spdlog::trace("  Lepton Beam Weight = {}", lwgt);
    spdlog::trace("  Hadron Beam Weight = {}", hwgt);
    spdlog::trace("  Phase Space Weight = {}", mwgt);
    spdlog::trace("  Weight = {}", wgt);

    return wgt;
}
