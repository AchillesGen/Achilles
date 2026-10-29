// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/QEFermiGas.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ReferenceHandler.hh"
#include "Achilles/Spinor.hh"
#include "Achilles/Units.hh"

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <cmath>

using achilles::QEFermiGas;

QEFermiGas::QEFermiGas(const YAML::Node &config, const YAML::Node &form_factor,
                       FormFactorBuilder &builder)
    : NuclearModel(form_factor, builder),
      m_ward{ToEnum(config["NuclearModel"]["Ward"].as<std::string>())} {
    const auto fg = config["NuclearModel"]["FermiGas"];
    if(!fg) throw std::runtime_error("QEFermiGas: missing 'FermiGas' block");

    const auto dist = fg["Distribution"].as<std::string>();
    if(dist == "RFG")
        m_distribution = QEFermiGas::MomentumDistribution::RFG;
    else if(dist == "CFG")
        m_distribution = QEFermiGas::MomentumDistribution::CFG;
    else
        throw std::runtime_error(
            fmt::format("QEFermiGas: unknown Distribution '{}' (RFG|CFG)", dist));

    m_kf = fg["kF"].as<double>();
    m_eb = fg["Binding"].as<double>();
    if(fg["PauliBlocking"]) m_pauli = fg["PauliBlocking"].as<bool>();

    if(m_distribution == QEFermiGas::MomentumDistribution::CFG) {
        m_c0 = fg["c0"].as<double>();
        m_lambda = fg["lambda"].as<double>();
        if(m_lambda <= 1) throw std::runtime_error("QEFermiGas: CFG requires lambda > 1");
        m_alpha0 = 1 - (1 - 1 / m_lambda) * m_c0 / (M_PI * M_PI);
        m_alpha1 = m_c0 * std::pow(m_kf, 4) / (3 * M_PI * M_PI);
    }

    ReferenceHandler::Handle().AddReference(
        {ReferenceType::inspire, "Smith:1972xh", "Relativistic Fermi gas"});
    if(m_distribution == QEFermiGas::MomentumDistribution::CFG)
        ReferenceHandler::Handle().AddReference(
            {ReferenceType::inspire, "Hen:2014yfa", "Correlated Fermi gas"});

    spdlog::info("QEFermiGas: {} with kF = {} MeV, Eb = {} MeV, Pauli blocking {}", dist, m_kf,
                 m_eb, m_pauli ? "on" : "off");
    if(m_distribution == QEFermiGas::MomentumDistribution::CFG)
        spdlog::info("QEFermiGas: c0 = {}, lambda = {}, alpha0 = {}, tail fraction = {}", m_c0,
                     m_lambda, m_alpha0, 1 - m_alpha0);
}

std::unique_ptr<achilles::NuclearModel> QEFermiGas::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<QEFermiGas>(config, form_factor, FormFactorBuilder::Instance());
}

std::string QEFermiGas::InspireHEP() const {
    return m_distribution == QEFermiGas::MomentumDistribution::CFG ? "Hen:2014yfa" : "Smith:1972xh";
}

YAML::Node QEFermiGas::PhaseSpaceSettings() const {
    YAML::Node settings;
    settings["PCut"] = PCut();
    settings["Binding"] = m_eb;
    return settings;
}

double QEFermiGas::Occupation(double k) const {
    if(k < m_kf) return m_alpha0; // = 1 for RFG
    if(m_distribution == QEFermiGas::MomentumDistribution::CFG && k < m_lambda * m_kf)
        return m_alpha1 / std::pow(k, 4);
    return 0;
}

double QEFermiGas::Density(double p) const {
    return 3.0 / (4 * M_PI * std::pow(m_kf, 3)) * Occupation(p);
}

std::string QEFermiGas::PhaseSpace(PID nuc_id) const {
    if(nuc_id != PID::hydrogen() && nuc_id != PID::free_neutron()) return PSName();
    if(nuc_id == PID::hydrogen())
        is_hydrogen = true;
    else
        is_free_neutron = true;
    return Coherent::Name();
}

// NOTE: the second argument is the list of spectators (see XSecBackend::InitialStateFactor),
// not the outgoing nucleon, so Pauli blocking is applied in CalcCurrents instead.
double QEFermiGas::InitialStateWeight(const std::vector<Particle> &had_in,
                                      const std::vector<Particle> &, size_t nprotons,
                                      size_t nneutrons) const {
    if(is_hydrogen) return had_in[0].ID() == PID::proton() ? 1 : 0;
    if(is_free_neutron) return had_in[0].ID() == PID::neutron() ? 1 : 0;

    const double nucleons = had_in[0].ID() == PID::proton() ? static_cast<double>(nprotons)
                                                            : static_cast<double>(nneutrons);
    return nucleons * Density(had_in[0].Momentum().P());
}

double QEFermiGas::PauliFactor(double pout) const {
    if(!m_pauli || is_hydrogen || is_free_neutron) return 1;
    return 1 - Occupation(pout);
}

// Identical to QESpectral::CalcCurrents: nucleon put on shell, binding moved into q0
achilles::NuclearModel::Currents QEFermiGas::CalcCurrents(const std::vector<Particle> &had_in,
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
    // Pauli blocking: |J|^2 -> [1 - occupation(|p'|)] |J|^2, so scale the amplitude by the sqrt
    const double pauli = PauliFactor(pOut.P());
    if(pauli <= 0) return results;
    const double pauli_amp = std::sqrt(pauli);

    pIn.E() = free_energy;
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pOut);
    ubar[1] = UBarSpinor(1, pOut);
    u[0] = USpinor(-1, -pIn);
    u[1] = USpinor(1, -pIn);

    for(const auto &formFactor : ff) {
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        auto current = HadronicCurrent(ubar, u, qVec, ffVal);
        for(auto &subcur : current) {
            for(auto &comp : subcur) comp *= pauli_amp;
            switch(m_ward) {
            case WardGauge::None:
                continue;
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

achilles::NuclearModel::Current QEFermiGas::HadronicCurrent(const std::array<Spinor, 2> &ubar,
                                                            const std::array<Spinor, 2> &u,
                                                            const FourVector &qVec,
                                                            const FormFactorMap &ffVal) const {
    Current result;
    std::array<SpinMatrix, 4> gamma{};
    for(size_t mu = 0; mu < 4; ++mu) {
        gamma[mu] = ffVal.at(FormFactorInfo::Type::F1) * SpinMatrix::GammaMu(mu);
        gamma[mu] +=
            ffVal.at(FormFactorInfo::Type::FA) * SpinMatrix::GammaMu(mu) * SpinMatrix::Gamma_5();
        gamma[mu] +=
            ffVal.at(FormFactorInfo::Type::FAP) * SpinMatrix::Gamma_5() * qVec[mu] / Constant::mN;
        double sign = 1;
        for(size_t nu = 0; nu < 4; ++nu) {
            gamma[mu] += std::complex<double>(0, 1) *
                         (ffVal.at(FormFactorInfo::Type::F2) * SpinMatrix::SigmaMuNu(mu, nu) *
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
