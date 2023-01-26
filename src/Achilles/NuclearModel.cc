#include "Achilles/NuclearModel.hh"
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Spinor.hh"
#include "Achilles/Particle.hh"

using achilles::NuclearModel;
using achilles::Coherent;
using achilles::QESpectral;

NuclearModel::NuclearModel(const YAML::Node& config, const std::shared_ptr<Nucleus> &nuc,
                           FormFactorBuilder &ffbuilder = FormFactorBuilder::Instance()) : m_nucleus{nuc} {
    spdlog::debug("Setting up form factors");
    const auto vectorFF = config["vector"].as<std::string>();
    const auto axialFF = config["axial"].as<std::string>();
    const auto coherentFF = config["coherent"].as<std::string>();
    m_form_factor = ffbuilder.Vector(vectorFF, config[vectorFF])
                             .AxialVector(axialFF, config[axialFF])
                             .Coherent(coherentFF, config[coherentFF])
                             .build();
}

NuclearModel::FormFactorArray NuclearModel::CouplingsFF(const FormFactor::Values &formFactors,
                                                        const std::vector<FormFactorInfo> &ffInfo) const {
    FormFactorArray results{};

    for(const auto & ff : ffInfo) {
        spdlog::trace("Form Factor: {}, Coupling: {}", ff.form_factor, ff.coupling);
        switch(ff.form_factor) {
            case FormFactorInfo::Type::F1p:
                results[0] += formFactors.F1p*ff.coupling;
                break;
            case FormFactorInfo::Type::F1n:
                results[0] += formFactors.F1n*ff.coupling;
                break;
            case FormFactorInfo::Type::F2p:
                results[1] += formFactors.F2p*ff.coupling;
                break;
            case FormFactorInfo::Type::F2n:
                results[1] += formFactors.F2n*ff.coupling;
                break;
            case FormFactorInfo::Type::FA:
                results[2] += formFactors.FA*ff.coupling;
                break;
            case FormFactorInfo::Type::FCoh:
                results[3] += formFactors.Fcoh*ff.coupling;
                break;
        }
    }

    return results;
}

YAML::Node NuclearModel::LoadFormFactor(const YAML::Node &config) {
    return YAML::LoadFile(config["NuclearModel"]["FormFactorFile"].as<std::string>());
}

achilles::Process_Group NuclearModel::AllowedStates(Process_Info info) {
    // TODO: Fill info and split for each nuclear state
    // TODO: Figure out how to combine with multiple incoming beams
    // Check for charge conservation
    int charge = -ParticleInfo(info.ids[0]).IntCharge();
    for(size_t i = 1; i < info.ids.size(); ++i) {
        charge += ParticleInfo(info.ids[i]).IntCharge();
    }
    charge /= 3;
    spdlog::debug("Charge = {}", charge);

    switch(Mode()) {
        case NuclearMode::None:
            throw std::runtime_error("NuclearModel: Invalid mode. Define custom AllowedStates");
            break;
        case NuclearMode::Coherent:
            if(charge != 0)
                throw std::runtime_error(fmt::format("Coherent: Requires charge 0, but found charge {}", charge));
          
            info.state = {{m_nucleus->ID()}, {m_nucleus->ID()}};
            m_group.AddProcess(info);
            break;
        case NuclearMode::Quasielastic:
            if(std::abs(charge) > 1)
                throw std::runtime_error(fmt::format("Quasielastic: Requires |charge| < 2, but found |charge| {}", std::abs(charge)));
            
            switch(charge) {
                case -1: // Final state has less charge than initial
                    info.state = {{PID::neutron()}, {PID::proton()}}; 
                    m_group.AddProcess(info);
                    break; 
                case 0: // Same charge in inital and final
                    info.state = {{PID::proton()}, {PID::proton()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::neutron()}, {PID::neutron()}}; 
                    m_group.AddProcess(info);
                    break;
                case 1: // Final state has more charge than initial
                    info.state = {{PID::proton()}, {PID::neutron()}}; 
                    m_group.AddProcess(info);
                    break;
            }
            break;
        case NuclearMode::MesonExchangeCurrent:
            if(std::abs(charge) > 1)
                throw std::runtime_error(fmt::format("{}: Requires |charge| < 2, but found |charge| {}",
                                         ToString(Mode()), std::abs(charge)));
            
            switch(charge) {
                case -1: // Final state has less charge than initial
                    info.state = {{PID::neutron(), PID::proton()},
                                  {PID::proton(), PID::proton()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::proton(), PID::neutron()},
                                  {PID::proton(), PID::proton()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::neutron(), PID::neutron()},
                                  {PID::proton(), PID::neutron()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::neutron(), PID::neutron()},
                                  {PID::neutron(), PID::proton()}}; 
                    m_group.AddProcess(info);
                    break; 
                case 0: // Same charge in inital and final
                    info.state = {{PID::neutron(), PID::neutron()},
                                  {PID::neutron(), PID::neutron()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::proton(), PID::neutron()},
                                  {PID::proton(), PID::neutron()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::neutron(), PID::proton()},
                                  {PID::neutron(), PID::proton()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::proton(), PID::proton()},
                                  {PID::proton(), PID::proton()}}; 
                    m_group.AddProcess(info);
                    break;
                case 1: // Final state has more charge than initial
                    info.state = {{PID::proton(), PID::neutron()},
                                  {PID::neutron(), PID::neutron()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::neutron(), PID::proton()},
                                  {PID::neutron(), PID::neutron()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::proton(), PID::proton()},
                                  {PID::neutron(), PID::proton()}}; 
                    m_group.AddProcess(info);
                    info.state = {{PID::proton(), PID::proton()},
                                  {PID::proton(), PID::neutron()}}; 
                    m_group.AddProcess(info);
                    break;
            }
            break;
        // TODO: Implement remaining cases
        case NuclearMode::Interference_QE_MEC:
        case NuclearMode::Resonance:
        case NuclearMode::ShallowInelastic:
        case NuclearMode::DeepInelastic:
            throw std::runtime_error(fmt::format("NuclearModel: Allowed states for {} not implemented yet",
                                                 ToString(Mode())));
    }
    return m_group;
}

size_t NuclearModel::NSpins() const {
    size_t nspins = 1;
    // TODO: Make less dependent on process 0
    for(const auto &init_state : m_group.Processes()[0].state.first) {
        nspins *= ParticleInfo(init_state).NSpins();
    }
    for(const auto &final_state : m_group.Processes()[0].state.second) {
        nspins *= ParticleInfo(final_state).NSpins();
    }
    return nspins;
}

// TODO: Clean this up such that the nucleus isn't loaded twice
Coherent::Coherent(const YAML::Node&, const YAML::Node &form_factor,
                   const std::shared_ptr<Nucleus> &nuc,
                   FormFactorBuilder &builder = FormFactorBuilder::Instance()) 
        : NuclearModel(form_factor, nuc, builder) {}

std::vector<achilles::NuclearModel::Currents> Coherent::CalcCurrents(const Event &event,
                                                                   const std::vector<FFInfoMap> &ff) const {
    auto pIn = event.Momentum().front();
    auto pOut = event.Momentum()[2];
    auto qVec = event.Momentum()[1];
    for(size_t  i = 3; i < event.Momentum().size(); ++i) {
        qVec -= event.Momentum()[i];
    }

    auto ffVals = EvalFormFactor(-qVec.M2());

    // Calculate coherent contributions
    std::vector<Currents> results(1);
    for(const auto &formFactor : ff[2]) {
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        spdlog::trace("fcoh = {}", ffVal[3]);
        Current current;
        std::vector<std::complex<double>> subcur(4);
        for(size_t i = 0; i < subcur.size(); ++i) {
            subcur[i] = (pIn[i] + pOut[i])*ffVal[3];
        }
        current.push_back(subcur);
        results[0][formFactor.first] = current;
        spdlog::trace("HadronicCurrent[{}] = [{}, {}, {}, {}]", formFactor.first,
                      subcur[0], subcur[1], subcur[2], subcur[3]);
    }

    return results;
}

// void Coherent::AllowedStates(Process_Info &info) {
//     // Check for charge conservation
//     int charge = -ParticleInfo(info.ids[0]).IntCharge();
//     for(size_t i = 1; i < info.ids.size(); ++i) {
//         charge += ParticleInfo(info.ids[i]).IntCharge();
//     }
//     charge /= 3;
//     spdlog::debug("Charge = {}", charge);
//     if(charge != 0)
//         throw std::runtime_error(fmt::format("Coherent: Requires charge 0, but found charge {}", charge));
// 
//     info.states[{nucleus_pid}] = {nucleus_pid}; 
// }

// bool Coherent::FillNucleus(Event &event, const std::vector<double> &xsecs) const {
//     // Calculate total cross section
//     if(xsecs[0] == 0) return false;
//     event.SetMEWeight(xsecs[0]);
// 
//     // Remove all nucleons
//     event.CurrentNucleus() -> Nucleons().clear();
// 
//     // Setup initial and final state nucleus
//     Particle initial = Particle(nucleus_pid, event.Momentum().front());
//     initial.Status() = ParticleStatus::initial_state;
//     event.CurrentNucleus() -> Nucleons().push_back(initial);
//     Particle final(nucleus_pid, event.Momentum()[2]);
//     final.Status() = ParticleStatus::final_state;
//     event.CurrentNucleus() -> Nucleons().push_back(final);
// 
//     return true;
// }

std::unique_ptr<NuclearModel> Coherent::Construct(const YAML::Node &config, const std::shared_ptr<Nucleus> &nuc) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<Coherent>(config, form_factor, nuc);
}

// TODO: Clean this interface up
QESpectral::QESpectral(const YAML::Node &config, const YAML::Node &form_factor,
                       const std::shared_ptr<Nucleus> &nuc,
                       FormFactorBuilder &builder = FormFactorBuilder::Instance()) 
        : NuclearModel(form_factor, nuc, builder),
          spectral_proton{config["NuclearModel"]["SpectralP"].as<std::string>()},
          spectral_neutron{config["NuclearModel"]["SpectralN"].as<std::string>()} {
    b_ward = config["NuclearModel"]["Ward"].as<bool>();
}

std::vector<NuclearModel::Currents> QESpectral::CalcCurrents(const Event &event,
                                                             const std::vector<FFInfoMap> &ff) const {
    auto pIn = event.Momentum().front();
    auto pOut = event.Momentum()[2];
    auto qVec = event.Momentum()[1];
    for(size_t i = 3; i < event.Momentum().size(); ++i) {
        qVec -= event.Momentum()[i];
    }
    auto removal_energy = Constant::mN - pIn.E();
    auto free_energy = sqrt(pIn.P2() + Constant::mN2);
    auto ffVals = EvalFormFactor(-qVec.M2()/1.0_GeV/1.0_GeV);
    auto omega = qVec.E();
    qVec.E() = qVec.E() + pIn.E() - free_energy;

    std::vector<Currents> results(2);
    // TODO: Move to the phase space definition
    std::vector<double> spectral(2);
    spectral[0] = spectral_proton(pIn.P(), removal_energy);
    spectral[1] = spectral_neutron(pIn.P(), removal_energy);
    spdlog::debug("Spectral function: S_p({}, {}) = {}, S_n({}, {}) = {}",
                  removal_energy, pIn.P(), spectral[0],
                  removal_energy, pIn.P(), spectral[1]);
    // Setup spinors
    pIn.E() = free_energy;
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pOut);
    ubar[1] = UBarSpinor(1, pOut);
    u[0] = USpinor(-1, -pIn);
    u[1] = USpinor(1, -pIn);

    // Loop over proton and neutron
    for(size_t i = 0; i < results.size(); ++i) {
        // Calculate nucleon contributions
        for(const auto &formFactor : ff[i]) {
            auto ffVal = CouplingsFF(ffVals, formFactor.second);
            spdlog::debug("{}: f1 = {}, f2 = {}, fa = {}", i, ffVal[0], ffVal[1], ffVal[2]);
            auto current = HadronicCurrent(ubar, u, qVec, ffVal);
            for(auto &subcur : current) {
                for(auto &val : subcur) {
                    // TODO: Move this to phase space 
                    val *= sqrt(spectral[i]);
                }
                // Correct the Ward identity
                if(b_ward) subcur[3] = omega/qVec.P()*subcur[0];
            }
            results[i][formFactor.first] = current;
        }
    }
    return results;
}

// void QESpectral::AllowedStates(Process_Info &info) {
//     // Check for charge conservation
//     int charge = -ParticleInfo(info.ids[0]).IntCharge();
//     for(size_t i = 1; i < info.ids.size(); ++i) {
//         charge += ParticleInfo(info.ids[i]).IntCharge();
//     }
//     charge /= 3;
//     if(std::abs(charge) > 1)
//         throw std::runtime_error(fmt::format("Quasielastic: Requires |charge| < 2, but found |charge| {}", std::abs(charge)));
// 
//     switch(charge) {
//         case -1: // Final state has less charge than initial
//             info.states[{PID::neutron()}] = {PID::proton()}; 
//             break; 
//         case 0: // Same charge in inital and final
//             info.states[{PID::proton()}] = {PID::proton()}; 
//             info.states[{PID::neutron()}] = {PID::neutron()}; 
//             break;
//         case 1: // Final state has more charge than initial
//             info.states[{PID::proton()}] = {PID::neutron()}; 
//             break;
//     }
// }

// bool QESpectral::FillNucleus(Event &event, const std::vector<double> &xsecs) const {
//     // Calculate total cross section
//     for(size_t i = 0; i < event.CurrentNucleus() -> Nucleons().size(); ++i) {
//         auto current_nucleon = event.CurrentNucleus() -> Nucleons()[i];
//         if(current_nucleon.ID() == PID::proton()) {
//             event.MatrixElementWgt(i) = xsecs[0];
//         } else {
//             event.MatrixElementWgt(i) = xsecs[1];
//         }
//     }
//     if(!event.TotalCrossSection())
//         return false;
// 
//     return true;
// }

std::unique_ptr<NuclearModel> QESpectral::Construct(const YAML::Node &config,
                                                    const std::shared_ptr<Nucleus> &nuc) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<QESpectral>(config, form_factor, nuc);
}

NuclearModel::Current QESpectral::HadronicCurrent(const std::array<Spinor, 2> &ubar,
                                                  const std::array<Spinor, 2> &u,
                                                  const FourVector &qVec,
                                                  const FormFactorArray &ffVal) const {
    Current result;
    std::array<SpinMatrix, 4> gamma{};
    for(size_t mu = 0; mu < 4; ++mu) {
        gamma[mu] = ffVal[0]*SpinMatrix::GammaMu(mu) + ffVal[2]*SpinMatrix::GammaMu(mu)*SpinMatrix::Gamma_5();
        double sign = 1;
        for(size_t nu = 0; nu < 4; ++nu) {
            gamma[mu] += std::complex<double>(0, 1)*(ffVal[1]*SpinMatrix::SigmaMuNu(mu, nu)*sign*qVec[nu]/(2*Constant::mN));
            sign = -1;
        }
    }

    for(size_t i = 0; i < 2; ++i) {
        for(size_t j = 0; j < 2; ++j) {
            std::vector<std::complex<double>> subcur(4);
            for(size_t mu = 0; mu < 4; ++mu) {
                subcur[mu] = ubar[i]*gamma[mu]*u[j];
            }
            result.push_back(subcur);
        }
    }
    return result;
}
