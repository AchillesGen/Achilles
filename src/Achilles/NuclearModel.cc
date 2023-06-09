#include "Achilles/NuclearModel.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/Poincare.hh"
#include "Achilles/Spinor.hh"

using achilles::Coherent;
using achilles::NuclearModel;
using achilles::QESpectral;
using Type = achilles::FormFactorInfo::Type;

NuclearModel::NuclearModel(const YAML::Node &config,
                           FormFactorBuilder &ffbuilder = FormFactorBuilder::Instance()) {
    spdlog::debug("Setting up form factors");
    const auto vectorFF = config["vector"].as<std::string>();
    const auto axialFF = config["axial"].as<std::string>();
    const auto coherentFF = config["coherent"].as<std::string>();
    m_form_factor = ffbuilder.Vector(vectorFF, config[vectorFF])
                        .AxialVector(axialFF, config[axialFF])
                        .Coherent(coherentFF, config[coherentFF])
                        .build();
    ffbuilder.Reset();
}

NuclearModel::FormFactorMap
NuclearModel::CouplingsFF(const FormFactor::Values &formFactors,
                          const std::vector<FormFactorInfo> &ffInfo) const {
    FormFactorMap results{};

    for(const auto &ff : ffInfo) {
        spdlog::trace("Form Factor: {}, Coupling: {}", ff.form_factor, ff.coupling);
        switch(ff.form_factor) {
        case Type::F1p:
            results[Type::F1] += formFactors.F1p * ff.coupling;
            break;
        case Type::F1n:
            results[Type::F1] += formFactors.F1n * ff.coupling;
            break;
        case Type::F2p:
            results[Type::F2] += formFactors.F2p * ff.coupling;
            break;
        case Type::F2n:
            results[Type::F2] += formFactors.F2n * ff.coupling;
            break;
        case Type::FA:
            results[Type::FA] += formFactors.FA * ff.coupling;
            break;
        case Type::FCoh:
            results[Type::FCoh] += formFactors.Fcoh * ff.coupling;
            break;
        case Type::F1:
        case Type::F2:
            throw std::runtime_error("Types F1 and F2 are reserved for nuclear models");
        }
    }

    return results;
}

YAML::Node NuclearModel::LoadFormFactor(const YAML::Node &config) {
    return YAML::LoadFile(config["NuclearModel"]["FormFactorFile"].as<std::string>());
}

NuclearModel::ModelMap achilles::LoadModels(const YAML::Node &node) {
    NuclearModel::ModelMap models;

    std::set<std::string> model_names;
    for(const auto &model_config : node["NuclearModels"]) {
        const auto name = model_config["NuclearModel"]["Model"].as<std::string>();
        if(model_names.find(name) != model_names.end()) {
            auto msg = fmt::format(
                "NuclearModels: Multiple definitions for model {}, skipping second model", name);
            spdlog::warn(msg);
            continue;
        }
        model_names.insert(name);
        auto model = NuclearModelFactory::Initialize(name, model_config);
        if(models.find(model->Mode()) != models.end()) {
            auto msg = fmt::format("NuclearModels: Multiple nuclear models for mode {} defined!",
                                   ToString(model->Mode()));
            throw std::runtime_error(msg);
        }
        models[model->Mode()] = std::move(model);
    }

    return models;
}

// TODO: Clean this up such that the nucleus isn't loaded twice
Coherent::Coherent(const YAML::Node &config, const YAML::Node &form_factor,
                   FormFactorBuilder &builder = FormFactorBuilder::Instance())
    : NuclearModel(form_factor, builder) {
    nucleus_pid = config["NuclearModel"]["Nucleus"].as<PID>();
}

achilles::NuclearModel::Currents Coherent::CalcCurrents(const std::vector<FourVector> &had_in,
                                                        const std::vector<FourVector> &had_out,
                                                        const FourVector &qVec,
                                                        const FFInfoMap &ff) const {
    auto pIn = had_in[0];
    auto pOut = had_out[0];

    auto ffVals = EvalFormFactor(-qVec.M2());

    // Calculate coherent contributions
    Currents results;
    for(const auto &formFactor : ff) {
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        spdlog::trace("fcoh = {}", ffVal[Type::FCoh]);

        Current current;
        VCurrent subcur;
        for(size_t i = 0; i < subcur.size(); ++i) {
            subcur[i] = (pIn[i] + pOut[i]) * ffVal[Type::FCoh];
        }
        current.push_back(subcur);
        results[formFactor.first] = current;
        spdlog::trace("HadronicCurrent[{}] = [{}, {}, {}, {}]", formFactor.first, subcur[0],
                      subcur[1], subcur[2], subcur[3]);
    }

    return results;
}

// TODO: Should return a process group
std::vector<achilles::ProcessInfo> Coherent::AllowedStates(const ProcessInfo &info) const {
    // Check for charge conservation
    const auto charge = info.LeptonicCharge();
    spdlog::debug("Charge = {}", charge);
    if(charge != 0)
        throw std::runtime_error(
            fmt::format("Coherent: Requires charge 0, but found charge {}", charge));

    auto result = info;
    result.m_hadronic = {{nucleus_pid}, {nucleus_pid}};
    return {result};
}

std::unique_ptr<NuclearModel> Coherent::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<Coherent>(config, form_factor);
}

// TODO: Clean this interface up
QESpectral::QESpectral(const YAML::Node &config, const YAML::Node &form_factor,
                       FormFactorBuilder &builder = FormFactorBuilder::Instance())
    : NuclearModel(form_factor, builder),
      m_ward{ToEnum(config["NuclearModel"]["Ward"].as<std::string>())},
      spectral_proton{config["NuclearModel"]["SpectralP"].as<std::string>()},
      spectral_neutron{config["NuclearModel"]["SpectralN"].as<std::string>()} {}

NuclearModel::Currents QESpectral::CalcCurrents(const std::vector<FourVector> &had_in,
                                                const std::vector<FourVector> &had_out,
                                                const FourVector &q, const FFInfoMap &ff) const {
    auto pIn = had_in[0];
    auto pOut = had_out[0];
    auto qVec = q;
    auto removal_energy = Constant::mN - pIn.E();
    auto free_energy = sqrt(pIn.P2() + Constant::mN2);
    auto ffVals = EvalFormFactor(-qVec.M2() / 1.0_GeV / 1.0_GeV);
    auto omega = qVec.E();
    qVec.E() = qVec.E() + pIn.E() - free_energy;

    Currents results;
    // TODO: Move to the phase space definition
    std::vector<double> spectral(2);
    spectral[0] = spectral_proton(pIn.P(), removal_energy);
    spectral[1] = spectral_neutron(pIn.P(), removal_energy);
    spdlog::debug("Spectral function: S_p({}, {}) = {}, S_n({}, {}) = {}", removal_energy, pIn.P(),
                  spectral[0], removal_energy, pIn.P(), spectral[1]);
    // Setup spinors
    pIn.E() = free_energy;
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pOut);
    ubar[1] = UBarSpinor(1, pOut);
    u[0] = USpinor(-1, -pIn);
    u[1] = USpinor(1, -pIn);

    // Calculate nucleon contributions
    for(const auto &formFactor : ff) {
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        spdlog::debug("f1 = {}, f2 = {}, fa = {}", ffVal[Type::F1], ffVal[Type::F2],
                      ffVal[Type::FA]);
        auto current = HadronicCurrent(ubar, u, qVec, ffVal);
        for(auto &subcur : current) {
            // Correct the Ward identity
            switch(m_ward) {
            case WardGauge::None:
                continue;
                break;
            case WardGauge::Coulomb:
                CoulombGauge(subcur, q, omega);
                break;
            case WardGauge::Weyl:
                WeylGauge(subcur, q, omega);
                break;
            case WardGauge::Landau:
                LandauGauge(subcur, q);
                break;
            }
        }
        results[formFactor.first] = current;
    }
    return results;
}

void QESpectral::CoulombGauge(VCurrent &cur, const FourVector &q, double omega) const {
    FourVector cur4_real{cur[0].real(), cur[1].real(), cur[2].real(), cur[3].real()};
    FourVector cur4_imag{cur[0].imag(), cur[1].imag(), cur[2].imag(), cur[3].imag()};
    FourVector ref{0, 0, 0, 1};
    Poincare poincare(q, ref, 0);

    poincare.Rotate(cur4_real);
    poincare.Rotate(cur4_imag);

    cur4_real[3] = omega / q.P() * cur[0].real();
    cur4_imag[3] = omega / q.P() * cur[0].imag();
    poincare.RotateBack(cur4_real);
    poincare.RotateBack(cur4_imag);
    // spdlog::info("Before: {}, {}, {}, {}",
    //              cur[0], cur[1], cur[2], cur[3]);
    for(size_t i = 0; i < 4; ++i) cur[i] = {cur4_real[i], cur4_imag[i]};
    // spdlog::info("After: {}, {}, {}, {}",
    //              cur[0], cur[1], cur[2], cur[3]);
}

void QESpectral::WeylGauge(VCurrent &cur, const FourVector &q, double omega) const {
    FourVector cur4_real{cur[0].real(), cur[1].real(), cur[2].real(), cur[3].real()};
    FourVector cur4_imag{cur[0].imag(), cur[1].imag(), cur[2].imag(), cur[3].imag()};
    FourVector ref{0, 0, 0, 1};
    Poincare poincare(q, ref, 0);

    poincare.Rotate(cur4_real);
    poincare.Rotate(cur4_imag);

    cur4_real[0] = q.P() / omega * cur4_real[3];
    cur4_imag[0] = q.P() / omega * cur4_imag[3];
    poincare.RotateBack(cur4_real);
    poincare.RotateBack(cur4_imag);
    // spdlog::info("Before: {}, {}, {}, {}",
    //              cur[0], cur[1], cur[2], cur[3]);
    for(size_t i = 0; i < 4; ++i) cur[i] = {cur4_real[i], cur4_imag[i]};
    // spdlog::info("After: {}, {}, {}, {}",
    //              cur[0], cur[1], cur[2], cur[3]);
}

void QESpectral::LandauGauge(VCurrent &cur, const FourVector &q) const {
    // spdlog::info("Before: {}, {}, {}, {}",
    //              cur[0], cur[1], cur[2], cur[3]);
    auto jdotq = cur[0] * q[0] - cur[1] * q[1] - cur[2] * q[2] - cur[3] * q[3];
    auto Q2 = -q.M2();
    for(size_t i = 0; i < 4; ++i) cur[i] += jdotq / Q2 * q[i];
    // spdlog::info("After: {}, {}, {}, {}",
    //              cur[0], cur[1], cur[2], cur[3]);
}

// TODO: Should return a process group
std::vector<achilles::ProcessInfo> QESpectral::AllowedStates(const ProcessInfo &info) const {
    // Check for charge conservation
    const auto charge = info.LeptonicCharge();
    if(std::abs(charge) > 1)
        throw std::runtime_error(fmt::format(
            "Quasielastic: Requires |charge| < 2, but found |charge| {}", std::abs(charge)));

    std::vector<ProcessInfo> results;
    auto local = info;
    switch(charge) {
    case -1: // Final state has less charge than initial
        local.m_hadronic = {{PID::neutron()}, {PID::proton()}};
        results.push_back(local);
        break;
    case 0: // Same charge in inital and final
        local.m_hadronic = {{PID::neutron()}, {PID::neutron()}};
        results.push_back(local);
        local.m_hadronic = {{PID::proton()}, {PID::proton()}};
        results.push_back(local);
        break;
    case 1: // Final state has more charge than initial
        local.m_hadronic = {{PID::proton()}, {PID::neutron()}};
        results.push_back(local);
        break;
    }

    return results;
}

std::unique_ptr<NuclearModel> QESpectral::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<QESpectral>(config, form_factor);
}

double QESpectral::InitialStateWeight(const std::vector<PID> &nucleons,
                                      const std::vector<FourVector> &mom) const {
    const double removal_energy = Constant::mN - mom[0].E();
    return nucleons[0] == PID::proton() ? spectral_proton(mom[0].P(), removal_energy)
                                        : spectral_neutron(mom[0].P(), removal_energy);
}

NuclearModel::Current QESpectral::HadronicCurrent(const std::array<Spinor, 2> &ubar,
                                                  const std::array<Spinor, 2> &u,
                                                  const FourVector &qVec,
                                                  const FormFactorMap &ffVal) const {
    Current result;
    std::array<SpinMatrix, 4> gamma{};
    auto mpi2 = pow(ParticleInfo(211).Mass(), 2);
    auto ffAP = 2.0 * Constant::mN2 / (-qVec.M2() + mpi2) * ffVal.at(Type::FA);
    for(size_t mu = 0; mu < 4; ++mu) {
        gamma[mu] = ffVal.at(Type::F1) * SpinMatrix::GammaMu(mu);
        gamma[mu] += ffVal.at(Type::FA) * SpinMatrix::GammaMu(mu) * SpinMatrix::Gamma_5();
        gamma[mu] += ffAP * SpinMatrix::Gamma_5() * qVec[mu] / Constant::mN;
        double sign = 1;
        for(size_t nu = 0; nu < 4; ++nu) {
            gamma[mu] +=
                std::complex<double>(0, 1) * (ffVal.at(Type::F2) * SpinMatrix::SigmaMuNu(mu, nu) *
                                              sign * qVec[nu] / (2 * Constant::mN));
            sign = -1;
        }
    }

    for(size_t i = 0; i < 2; ++i) {
        for(size_t j = 0; j < 2; ++j) {
            VCurrent subcur;
            for(size_t mu = 0; mu < 4; ++mu) { subcur[mu] = ubar[i] * gamma[mu] * u[j]; }
            result.push_back(subcur);
        }
    }

    return result;
}
