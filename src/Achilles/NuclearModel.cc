#include "Achilles/Settings.hh"
#include "Achilles/System.hh"
#include <stdexcept>
#ifdef ACHILLES_EVENT_DETAILS
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#endif

#include "Achilles/Exception.hh"
#include "Achilles/Exceptions.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Poincare.hh"
#include "Achilles/Process.hh"
#include "Achilles/Spinor.hh"

#include "fmt/std.h"

using achilles::Coherent;
using achilles::HyperonSpectral;
using achilles::NuclearModel;
using achilles::QESpectral;
using Type = achilles::FormFactorInfo::Type;

NuclearModel::NuclearModel(const YAML::Node &config,
                           FormFactorBuilder &ffbuilder = FormFactorBuilder::Instance()) {
    spdlog::debug("Setting up form factors");
    const auto vectorFF = config["vector"].as<std::string>();
    const auto axialFF = config["axial"].as<std::string>();
    const auto coherentFF = config["coherent"].as<std::string>();
    const auto resvectorFF = config["resonancevector"].as<std::string>();
    const auto resaxialFF = config["resonanceaxial"].as<std::string>();
    const auto mecvectorFF = config["mecvector"].as<std::string>();
    const auto mecaxialFF = config["mecaxial"].as<std::string>();
    const auto hyperonFF = config["hyperon"].as<std::string>();
    m_form_factor = ffbuilder.Vector(vectorFF, config[vectorFF])
                        .AxialVector(axialFF, config[axialFF])
                        .Coherent(coherentFF, config[coherentFF])
                        .ResonanceVector(resvectorFF, config[resvectorFF])
                        .ResonanceAxial(resaxialFF, config[resaxialFF])
                        .MesonExchangeVector(mecvectorFF, config[mecvectorFF])
                        .MesonExchangeAxial(mecaxialFF, config[mecaxialFF])
                        .Hyperon(hyperonFF, config[hyperonFF])
                        .build();
    ffbuilder.Reset();
}

void NuclearModel::SetTransform() {
    spdlog::debug("Setting up frame transformation: {}", static_cast<int>(Frame()));
    switch(Frame()) {
    case NuclearFrame::Lab:
        transform = std::bind(&NuclearModel::TransformLab, this, std::placeholders::_1,
                              std::placeholders::_2, std::placeholders::_3);
        break;
    case NuclearFrame::QZ:
        transform = std::bind(&NuclearModel::TransformQZ, this, std::placeholders::_1,
                              std::placeholders::_2, std::placeholders::_3);
        break;
    case NuclearFrame::Custom:
        transform = std::bind(&NuclearModel::TransformCustom, this, std::placeholders::_1,
                              std::placeholders::_2, std::placeholders::_3);
        break;
    default:
        throw std::runtime_error("NuclearModel: Invalid frame");
    }
}

NuclearModel::FormFactorMap
NuclearModel::CouplingsFF(const FormFactor::Values &formFactors,
                          const std::vector<FormFactorInfo> &ffInfo) const {
    FormFactorMap results{};
    results[Type::F1] = 0;
    results[Type::F2] = 0;
    results[Type::FA] = 0;
    results[Type::FAP] = 0;
    results[Type::FCoh] = 0;

    for(const auto &ff : ffInfo) {
        SPDLOG_TRACE("Form Factor: {}, Coupling: {}", ff.form_factor, ff.coupling);
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
        case Type::FAP:
            results[Type::FA] += formFactors.FA * ff.coupling;
            results[Type::FAP] += formFactors.FAP * ff.coupling;
            break;
        case Type::FCoh:
            results[Type::FCoh] += formFactors.Fcoh * ff.coupling;
            break;
        case Type::FResV:
            results[Type::FResV] += formFactors.FresV * ff.coupling;
            break;
        case Type::FResA:
            results[Type::FResA] += formFactors.FresA * ff.coupling;
            break;
        case Type::FPiEM:
            results[Type::FPiEM] += formFactors.Fpiem * ff.coupling;
            break;
        case Type::FMecV3:
            results[Type::FMecV3] += formFactors.FmecV3 * ff.coupling;
            break;
        case Type::FMecV4:
            results[Type::FMecV4] += formFactors.FmecV4 * ff.coupling;
            break;
        case Type::FMecV5:
            results[Type::FMecV5] += formFactors.FmecV5 * ff.coupling;
            break;
        case Type::FMecA5:
            results[Type::FMecA5] += formFactors.FmecA5 * ff.coupling;
            break;
        case Type::F1Lam:
            results[Type::F1Lam] += formFactors.F1lam * ff.coupling;
            break;
        case Type::F2Lam:
            results[Type::F2Lam] += formFactors.F2lam * ff.coupling;
            break;
        case Type::FALam:
            results[Type::FALam] += formFactors.FAlam * ff.coupling;
            break;
        case Type::F1Sigm:
            results[Type::F1Sigm] += formFactors.F1sigm * ff.coupling;
            break;
        case Type::F2Sigm:
            results[Type::F2Sigm] += formFactors.F2sigm * ff.coupling;
            break;
        case Type::FASigm:
            results[Type::FASigm] += formFactors.FAsigm * ff.coupling;
            break;
        case Type::F1Sig0:
            results[Type::F1Sig0] += formFactors.F1sig0 * ff.coupling;
            break;
        case Type::F2Sig0:
            results[Type::F2Sig0] += formFactors.F2sig0 * ff.coupling;
            break;
        case Type::FASig0:
            results[Type::FASig0] += formFactors.FAsig0 * ff.coupling;
            break;
        case Type::F1:
        case Type::F2:
            throw std::runtime_error("Types F1 and F2 are reserved for nuclear models");
        }
    }

    return results;
}

void NuclearModel::TransformFrame(Event &event, const Process &process, bool forward) const {
    transform(event, process, forward);
}

void NuclearModel::TransformQZ(Event &event, const Process &process, bool forward) {
    if(forward) {
        auto q = process.ExtractQ(event);
        rotation = q.AlignZ();
        for(auto &mom : event.Momentum()) { mom = mom.Rotate(rotation); }
    } else {
        for(auto &mom : event.Momentum()) { mom = mom.RotateBack(rotation); }
    }
}

YAML::Node NuclearModel::LoadFormFactor(const YAML::Node &config) {
    std::string filename = config["NuclearModel"]["FormFactorFile"].as<std::string>();
    try {
        return YAML::LoadFile(Filesystem::FindFile(filename, "NuclearModel"));
    } catch(const AchillesLoadError &e) {
        spdlog::warn("NuclearModel: Copying and using default Form Factors file from {} as "
                     "FormFactorsDefault.yml",
                     PathVariables::installDefaults / "FormFactors.yml");
        fs::copy(PathVariables::installDefaults / "FormFactors.yml", "FormFactorsDefault.yml");
        return YAML::LoadFile("FormFactorsDefault.yml");
    }
}

YAML::Node NuclearModel::LoadModelParams(const YAML::Node &config) {
    std::string filename = Filesystem::FindFile(
        config["NuclearModel"]["ModelParamsFile"].as<std::string>(), "NuclearModel");
    return YAML::LoadFile(filename);
}

NuclearModel::ModelMap achilles::LoadModels(const Settings &settings) {
    NuclearModel::ModelMap models;

    std::set<std::string> model_names;
    for(const auto &model_config : settings["NuclearModels"]) {
        const auto name = model_config["NuclearModel"]["Model"].as<std::string>();
        // TODO: Move this logic to be better contained in the Settings or Factory classes
        // TODO: Figure out how to handle this for fortran models???
        try {
            auto model = NuclearModelFactory::Initialize(name, model_config);
            if(model_names.find(model->GetName()) != model_names.end()) {
                auto msg = fmt::format(
                    "NuclearModels: Multiple definitions for model {}, skipping second model",
                    model->GetName());
                spdlog::warn(msg);
                continue;
            }
            model_names.insert(model->GetName());
            if(models.find(model->Mode()) != models.end()) {
                auto msg =
                    fmt::format("NuclearModels: Multiple nuclear models for mode {} defined!",
                                ToString(model->Mode()));
                throw std::runtime_error(msg);
            }
            model->SetTransform();
            models[model->Mode()] = std::move(model);
        } catch(std::out_of_range &e) {
            spdlog::error("NuclearModel: Requested model \"{}\", did you mean \"{}\"", name,
                          achilles::GetSuggestion(achilles::NuclearModelFactory::List(), name));
            spdlog::error("NuclearModel: Run `achilles --display-nuc-models` to see all options");
            exit(-1);
        }
    }

    return models;
}

// TODO: Rewrite to match process grouping
std::vector<achilles::ProcessInfo> NuclearModel::AllowedStates(const ProcessInfo &info) const {
    // Check for charge conservation
    const auto charge = info.LeptonicCharge();
    spdlog::debug("Charge = {}", charge);
    std::vector<ProcessInfo> results;
    auto local = info;

    switch(Mode()) {
    case NuclearMode::None:
        throw std::runtime_error("NuclearModel: Invalid mode. Define custom AllowedStates.");

    case NuclearMode::Coherent:
        throw std::runtime_error("NuclearModel: Coherent needs custom AllowedStates.");

    case NuclearMode::Quasielastic:
        if(std::abs(charge) > 1)
            throw std::runtime_error(fmt::format(
                "Quasielastic: Requires |charge| < 2, but found |charge| {}", std::abs(charge)));

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

    case NuclearMode::Resonance:
        if(std::abs(charge) > 1)
            throw std::runtime_error(
                fmt::format("{}}: Requires |charge| < 2, but found |charge| {}", ToString(Mode()),
                            std::abs(charge)));

        switch(charge) {
        case -1: // Final state has less charge than initial
            local.m_hadronic = {{PID::neutron()}, {PID::proton(), PID::pion0()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron()}, {PID::neutron(), PID::pionp()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton()}, {PID::proton(), PID::pionp()}};
            results.push_back(local);
            break;
        case 0: // Same charge in inital and final
            local.m_hadronic = {{PID::proton()}, {PID::proton(), PID::pion0()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton()}, {PID::neutron(), PID::pionp()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron()}, {PID::neutron(), PID::pion0()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron()}, {PID::proton(), -PID::pionp()}};
            results.push_back(local);
            break;
        case 1: // Final state has more charge than initial
            local.m_hadronic = {{PID::proton()}, {PID::neutron(), PID::pion0()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton()}, {PID::proton(), -PID::pionp()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron()}, {PID::neutron(), -PID::pionp()}};
            results.push_back(local);
            break;
        }
        return results;

    case NuclearMode::MesonExchangeCurrent:
        if(std::abs(charge) > 1)
            throw std::runtime_error(fmt::format("{}: Requires |charge| < 2, but found |charge| {}",
                                                 ToString(Mode()), std::abs(charge)));

        switch(charge) {
        case -1: // Final state has less charge than initial
            local.m_hadronic = {{PID::neutron(), PID::proton()}, {PID::proton(), PID::proton()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton(), PID::neutron()}, {PID::proton(), PID::proton()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron(), PID::neutron()}, {PID::proton(), PID::neutron()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron(), PID::neutron()}, {PID::neutron(), PID::proton()}};
            results.push_back(local);
            break;
        case 0: // Same charge in inital and final
            local.m_hadronic = {{PID::neutron(), PID::neutron()}, {PID::neutron(), PID::neutron()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton(), PID::neutron()}, {PID::proton(), PID::neutron()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron(), PID::proton()}, {PID::neutron(), PID::proton()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton(), PID::proton()}, {PID::proton(), PID::proton()}};
            results.push_back(local);
            break;
        case 1: // Final state has more charge than initial
            local.m_hadronic = {{PID::proton(), PID::neutron()}, {PID::neutron(), PID::neutron()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron(), PID::proton()}, {PID::neutron(), PID::neutron()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton(), PID::proton()}, {PID::neutron(), PID::proton()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton(), PID::proton()}, {PID::proton(), PID::neutron()}};
            results.push_back(local);
            break;
        }
        return results;
    case NuclearMode::Interference_QE_MEC:
        if(std::abs(charge) > 1)
            throw std::runtime_error(fmt::format(
                "Interference: Requires |charge| < 2, but found |charge| {}", std::abs(charge)));

        switch(charge) {
        case -1: // Final state has less charge than initial
            local.m_hadronic = {{PID::neutron()}, {PID::proton()}};
            local.m_spectator = {PID::neutron()};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron()}, {PID::proton()}};
            local.m_spectator = {PID::proton()};
            results.push_back(local);
            break;
        case 0: // Same charge in inital and final
            local.m_hadronic = {{PID::neutron()}, {PID::neutron()}};
            local.m_spectator = {PID::neutron()};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron()}, {PID::neutron()}};
            local.m_spectator = {PID::proton()};
            results.push_back(local);
            local.m_hadronic = {{PID::proton()}, {PID::proton()}};
            local.m_spectator = {PID::neutron()};
            results.push_back(local);
            local.m_hadronic = {{PID::proton()}, {PID::proton()}};
            local.m_spectator = {PID::proton()};
            results.push_back(local);
            break;
        case 1: // Final state has more charge than initial
            local.m_hadronic = {{PID::proton()}, {PID::neutron()}};
            local.m_spectator = {PID::neutron()};
            results.push_back(local);
            local.m_hadronic = {{PID::proton()}, {PID::neutron()}};
            local.m_spectator = {PID::proton()};
            results.push_back(local);
            break;
        }

        return results;

    case NuclearMode::Hyperon:
        if(std::abs(charge) > 1)
            throw std::runtime_error(fmt::format(
                "Hyperon: Requires |charge| < 2, but found |charge| {}", std::abs(charge)));

        switch(charge) {
        case -1: // Final state has less charge than initial
            break;
        case 0: // Same charge in inital and final
            break;
        case 1: // Final state has more charge than initial
            local.m_hadronic = {{PID::proton()}, {PID::lambda0()}};
            results.push_back(local);
            local.m_hadronic = {{PID::proton()}, {PID::sigma0()}};
            results.push_back(local);
            local.m_hadronic = {{PID::neutron()}, {PID::sigmam()}};
            results.push_back(local);
            break;
        }

        return results;

    // TODO: Implement remaining cases
    case NuclearMode::ShallowInelastic:
    case NuclearMode::DeepInelastic:
        throw std::runtime_error(fmt::format(
            "NuclearModel: Allowed states for {} not implemented yet", ToString(Mode())));
    }

    return results;
}

// TODO: Clean this up so it isn't hardcoded
size_t NuclearModel::NSpins() const {
    switch(Mode()) {
    case NuclearMode::None:
        throw std::runtime_error("NuclearModel: Invalid mode. Define custom NSpins.");
    case NuclearMode::Coherent:
        return 1;
    case NuclearMode::Quasielastic:
    case NuclearMode::Interference_QE_MEC:
    case NuclearMode::Resonance:
    case NuclearMode::Hyperon:
        return 4;
    case NuclearMode::MesonExchangeCurrent:
        return 16;
    case NuclearMode::ShallowInelastic:
    case NuclearMode::DeepInelastic:
        throw std::runtime_error(
            fmt::format("NuclearModel: NSpins for {} not implemented yet", ToString(Mode())));
    }

    return 0;
}

void NuclearModel::CoulombGauge(VCurrent &cur, const FourVector &q, double omega) const {
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
    for(size_t i = 0; i < 4; ++i) cur[i] = {cur4_real[i], cur4_imag[i]};
}

void NuclearModel::WeylGauge(VCurrent &cur, const FourVector &q, double omega) const {
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
    for(size_t i = 0; i < 4; ++i) cur[i] = {cur4_real[i], cur4_imag[i]};
}

void NuclearModel::LandauGauge(VCurrent &cur, const FourVector &q) const {
    auto jdotq = cur[0] * q[0] - cur[1] * q[1] - cur[2] * q[2] - cur[3] * q[3];
    auto Q2 = -q.M2();
    for(size_t i = 0; i < 4; ++i) cur[i] += jdotq / Q2 * q[i];
}

// TODO: Clean this up such that the nucleus isn't loaded twice, and that it works with multiple
// nuclei
Coherent::Coherent(const YAML::Node &config, const YAML::Node &form_factor,
                   FormFactorBuilder &builder = FormFactorBuilder::Instance())
    : NuclearModel(form_factor, builder) {
    nucleus_pid = config["NuclearModel"]["Nucleus"].as<PID>();
}

achilles::NuclearModel::Currents Coherent::CalcCurrents(const std::vector<Particle> &had_in,
                                                        const std::vector<Particle> &had_out,
                                                        const std::vector<Particle> &,
                                                        const FourVector &qVec,
                                                        const FFInfoMap &ff) const {
    auto pIn = had_in[0].Momentum();
    auto pOut = had_out[0].Momentum();
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

std::string Coherent::PhaseSpace(PID nuc_id) const {
    if(nuc_id == PID::hydrogen() || nuc_id == PID::free_neutron()) throw achilles::InvalidChannel("Hydrogen or free neutron nucleus");
    if(nuc_id == nucleus_pid) return PSName();
    throw achilles::InvalidChannel(
        fmt::format("Nucleus don't match: Model {}, Phasespace {}", nucleus_pid, nuc_id));
}

// TODO: Clean this interface up
QESpectral::QESpectral(const YAML::Node &config, const YAML::Node &form_factor,
                       FormFactorBuilder &builder = FormFactorBuilder::Instance())
    : NuclearModel(form_factor, builder),
      m_ward{ToEnum(config["NuclearModel"]["Ward"].as<std::string>())},
      spectral_proton{config["NuclearModel"]["SpectralP"].as<std::string>()},
      spectral_neutron{config["NuclearModel"]["SpectralN"].as<std::string>()} {}

NuclearModel::Currents QESpectral::CalcCurrents(const std::vector<Particle> &had_in,
                                                const std::vector<Particle> &had_out,
                                                const std::vector<Particle> &, const FourVector &q,
                                                const FFInfoMap &ff) const {
    if(had_in[0].ID() == PID::neutron() && is_hydrogen) return {};
    if(had_in[0].ID() == PID::proton() && is_free_neutron) return {};

    auto pIn = had_in[0].Momentum();
    auto pOut = had_out[0].Momentum();
    auto qVec = q;
    auto free_energy = sqrt(pIn.P2() + Constant::mN2);
    auto ffVals = EvalFormFactor(-qVec.M2() / 1.0_GeV / 1.0_GeV);
    auto omega = qVec.E();
    qVec.E() = qVec.E() + pIn.E() - free_energy;

    Currents results;

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

std::unique_ptr<NuclearModel> QESpectral::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<QESpectral>(config, form_factor);
}

double QESpectral::InitialStateWeight(const std::vector<Particle> &nucleons,
                                      const std::vector<Particle> &, size_t nprotons,
                                      size_t nneutrons) const {
    if(is_hydrogen) return nucleons[0].ID() == PID::proton() ? 1 : 0;
    if(is_free_neutron) return nucleons[0].ID() == PID::neutron() ? 1 : 0;
    const double removal_energy = Constant::mN - nucleons[0].E();
    return nucleons[0].ID() == PID::proton()
               ? static_cast<double>(nprotons) *
                     spectral_proton(nucleons[0].Momentum().P(), removal_energy)
               : static_cast<double>(nneutrons) *
                     spectral_neutron(nucleons[0].Momentum().P(), removal_energy);
}

NuclearModel::Current QESpectral::HadronicCurrent(const std::array<Spinor, 2> &ubar,
                                                  const std::array<Spinor, 2> &u,
                                                  const FourVector &qVec,
                                                  const FormFactorMap &ffVal) const {
    Current result;
    std::array<SpinMatrix, 4> gamma{};
    for(size_t mu = 0; mu < 4; ++mu) {
        gamma[mu] = ffVal.at(Type::F1) * SpinMatrix::GammaMu(mu);
        gamma[mu] += ffVal.at(Type::FA) * SpinMatrix::GammaMu(mu) * SpinMatrix::Gamma_5();
        gamma[mu] += ffVal.at(Type::FAP) * SpinMatrix::Gamma_5() * qVec[mu] / Constant::mN;
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

std::string QESpectral::PhaseSpace(PID nuc_id) const {
    if(nuc_id != PID::hydrogen() && nuc_id != PID::free_neutron())  return PSName(); 
    if(nuc_id == PID::hydrogen()) is_hydrogen = true;
    else is_free_neutron = true;
    return Coherent::Name();
}

// TODO: Clean this interface up
HyperonSpectral::HyperonSpectral(const YAML::Node &config, const YAML::Node &form_factor,
                                 FormFactorBuilder &builder = FormFactorBuilder::Instance())
    : NuclearModel(form_factor, builder),
      m_ward{ToEnum(config["NuclearModel"]["Ward"].as<std::string>())},
      spectral_proton{config["NuclearModel"]["SpectralP"].as<std::string>()},
      spectral_neutron{config["NuclearModel"]["SpectralN"].as<std::string>()} {}

NuclearModel::Currents HyperonSpectral::CalcCurrents(const std::vector<Particle> &had_in,
                                                     const std::vector<Particle> &had_out,
                                                     const std::vector<Particle> &,
                                                     const FourVector &q,
                                                     const FFInfoMap &ff) const {
    if(had_in[0].ID() == PID::neutron() && is_hydrogen) return {};
    if(had_in[0].ID() == PID::proton() && is_free_neutron) return {};

    auto pIn = had_in[0].Momentum();
    auto pOut = had_out[0].Momentum();
    auto qVec = q;
    auto free_energy = sqrt(pIn.P2() + Constant::mN2);
    auto ffVals = EvalFormFactor(-qVec.M2() / 1.0_GeV / 1.0_GeV);
    auto omega = qVec.E();
    qVec.E() = qVec.E() + pIn.E() - free_energy;

    Currents results;

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
        auto current = HadronicCurrent(ubar, u, qVec, ffVal, had_out[0]);
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

std::unique_ptr<NuclearModel> HyperonSpectral::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<HyperonSpectral>(config, form_factor);
}

double HyperonSpectral::InitialStateWeight(const std::vector<Particle> &nucleons,
                                           const std::vector<Particle> &, size_t nprotons,
                                           size_t nneutrons) const {
    if(is_hydrogen) return nucleons[0].ID() == PID::proton() ? 1 : 0;
    if(is_free_neutron) return nucleons[0].ID() == PID::neutron() ? 1 : 0;
    const double removal_energy = Constant::mN - nucleons[0].E();
    return nucleons[0].ID() == PID::proton()
               ? static_cast<double>(nprotons) *
                     spectral_proton(nucleons[0].Momentum().P(), removal_energy)
               : static_cast<double>(nneutrons) *
                     spectral_neutron(nucleons[0].Momentum().P(), removal_energy);
}

NuclearModel::Current HyperonSpectral::HadronicCurrent(const std::array<Spinor, 2> &ubar,
                                                       const std::array<Spinor, 2> &u,
                                                       const FourVector &qVec,
                                                       const FormFactorMap &ffVal,
                                                       const Particle &hyperon) const {
    std::complex<double> F1hyp, F2hyp, FAhyp;
    if(hyperon.ID() == PID::lambda0()) {
        F1hyp = ffVal.at(Type::F1Lam);
        F2hyp = ffVal.at(Type::F2Lam);
        FAhyp = ffVal.at(Type::FALam);
    } else if(hyperon.ID() == PID::sigmam()) {
        F1hyp = ffVal.at(Type::F1Sigm);
        F2hyp = ffVal.at(Type::F2Sigm);
        FAhyp = ffVal.at(Type::FASigm);
    } else if(hyperon.ID() == PID::sigma0()) {
        F1hyp = ffVal.at(Type::F1Sig0);
        F2hyp = ffVal.at(Type::F2Sig0);
        FAhyp = ffVal.at(Type::FASig0);
    }

    Current result;
    std::array<SpinMatrix, 4> gamma{};
    for(size_t mu = 0; mu < 4; ++mu) {
        gamma[mu] = F1hyp * SpinMatrix::GammaMu(mu);
        gamma[mu] += FAhyp * SpinMatrix::GammaMu(mu) * SpinMatrix::Gamma_5();
        double sign = 1;
        for(size_t nu = 0; nu < 4; ++nu) {
            gamma[mu] += std::complex<double>(0, 1) * (F2hyp * SpinMatrix::SigmaMuNu(mu, nu) *
                                                       sign * qVec[nu] / (hyperon.Mass()));
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

std::string HyperonSpectral::PhaseSpace(PID nuc_id) const {
    if(nuc_id != PID::hydrogen() && nuc_id != PID::free_neutron())  return PSName(); 
    if(nuc_id == PID::hydrogen()) is_hydrogen = true;
    else is_free_neutron = true;
    return Coherent::Name();
}
