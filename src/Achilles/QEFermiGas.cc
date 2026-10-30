// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/QEFermiGas.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ReferenceHandler.hh"
#include "Achilles/Spinor.hh"
#include "Achilles/Units.hh"

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <algorithm>
#include <cmath>

using achilles::LocalFermiSea;
using achilles::QEFermiGas;

LocalFermiSea::LocalFermiSea(const std::vector<double> &edges, const std::vector<double> &rho_p,
                             const std::vector<double> &rho_n)
    : m_edges{edges}, m_rho_p{rho_p}, m_rho_n{rho_n} {
    if(edges.size() != rho_p.size() + 1 || rho_p.size() != rho_n.size())
        throw std::runtime_error("LocalFermiSea: need n+1 bin edges and n densities");
    const auto kf = [](double rho) { return Constant::HBARC * std::cbrt(3 * M_PI * M_PI * rho); };
    for(size_t i = 0; i < rho_p.size(); ++i) {
        m_dv.push_back(4 * M_PI / 3 * (std::pow(edges[i + 1], 3) - std::pow(edges[i], 3)));
        m_kf_p.push_back(kf(std::max(rho_p[i], 0.0)));
        m_kf_n.push_back(kf(std::max(rho_n[i], 0.0)));
        m_kf_max = std::max({m_kf_max, m_kf_p.back(), m_kf_n.back()});
    }
}

const std::vector<double> &LocalFermiSea::FermiMomenta(PID nucleon) const {
    return nucleon == PID::proton() ? m_kf_p : m_kf_n;
}

double LocalFermiSea::Volume(double k, PID nucleon) const {
    const auto &kf = FermiMomenta(nucleon);
    double volume = 0;
    for(size_t i = 0; i < m_dv.size(); ++i)
        if(kf[i] > k) volume += m_dv[i];
    return volume;
}

// Fraction of the volume occupied by (in, p) in which (out, pout) is above the local Fermi sea
double LocalFermiSea::PauliFactor(double p, PID in, double pout, PID out) const {
    const auto &kf_in = FermiMomenta(in);
    const auto &kf_out = FermiMomenta(out);
    double occupied = 0, allowed = 0;
    for(size_t i = 0; i < m_dv.size(); ++i) {
        if(kf_in[i] <= p) continue;
        occupied += m_dv[i];
        if(kf_out[i] <= pout) allowed += m_dv[i];
    }
    return occupied > 0 ? allowed / occupied : 0;
}

// Innermost and outermost bin edges where (in, p) is occupied and (out, pout) is free
std::optional<std::pair<double, double>> LocalFermiSea::AllowedRadii(double p, PID in, double pout,
                                                                     PID out) const {
    const auto &kf_in = FermiMomenta(in);
    const auto &kf_out = FermiMomenta(out);
    std::optional<std::pair<double, double>> radii;
    for(size_t i = 0; i < m_dv.size(); ++i) {
        if(kf_in[i] <= p || kf_out[i] > pout) continue;
        if(!radii) radii = std::make_pair(m_edges[i], m_edges[i + 1]);
        radii->second = m_edges[i + 1];
    }
    return radii;
}

double LocalFermiSea::NNucleons(PID nucleon) const {
    const auto &rho = nucleon == PID::proton() ? m_rho_p : m_rho_n;
    double n = 0;
    for(size_t i = 0; i < m_dv.size(); ++i) n += rho[i] * m_dv[i];
    return n;
}

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

    if(fg["kF"].as<std::string>() == "Local") {
        if(m_distribution != QEFermiGas::MomentumDistribution::RFG)
            throw std::runtime_error("QEFermiGas: kF: Local is only implemented for RFG");
        m_local = true;
    } else {
        m_kf = fg["kF"].as<double>();
    }
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

    spdlog::info("QEFermiGas: {} with kF = {}, Eb = {} MeV, Pauli blocking {}", dist,
                 m_local ? "local" : fmt::format("{} MeV", m_kf), m_eb, m_pauli ? "on" : "off");
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

// Local Fermi sea from the nucleus' proton and neutron densities, on 0.01 fm radial bins
void QEFermiGas::SetNucleus(const Nucleus &nucleus) {
    if(!m_local || nucleus.ID() == PID::hydrogen() || nucleus.ID() == PID::free_neutron()) return;
    constexpr double dr = 0.01;
    const auto nbins = static_cast<size_t>(std::ceil(nucleus.Radius() / dr));
    std::vector<double> edges(nbins + 1), rho_p(nbins), rho_n(nbins);
    for(size_t i = 0; i <= nbins; ++i) edges[i] = static_cast<double>(i) * dr;
    for(size_t i = 0; i < nbins; ++i) {
        const double r = 0.5 * (edges[i] + edges[i + 1]);
        rho_p[i] = nucleus.ProtonRho(r);
        rho_n[i] = nucleus.NeutronRho(r);
    }
    m_sea = LocalFermiSea(edges, rho_p, rho_n);
    spdlog::info("QEFermiGas: local Fermi sea of {}: max kF = {:.1f} MeV, Z = {:.3f}, N = {:.3f} "
                 "from the densities",
                 nucleus.ID(), m_sea.MaxFermiMomentum(), m_sea.NNucleons(PID::proton()),
                 m_sea.NNucleons(PID::neutron()));
}

std::optional<std::pair<double, double>> QEFermiGas::StruckNucleonRadii(double p, PID in,
                                                                        double pout,
                                                                        PID out) const {
    if(!m_local || is_hydrogen || is_free_neutron) return std::nullopt;
    return m_sea.AllowedRadii(p, in, pout, out);
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

    if(m_local)
        return 2 * m_sea.Volume(had_in[0].Momentum().P(), had_in[0].ID()) /
               std::pow(2 * M_PI * Constant::HBARC, 3);

    const double nucleons = had_in[0].ID() == PID::proton() ? static_cast<double>(nprotons)
                                                            : static_cast<double>(nneutrons);
    return nucleons * Density(had_in[0].Momentum().P());
}

double QEFermiGas::PauliFactor(double pout) const {
    if(!m_pauli || is_hydrogen || is_free_neutron) return 1;
    return 1 - Occupation(pout);
}

double QEFermiGas::PauliFactor(double p, PID in, double pout, PID out) const {
    if(!m_local) return PauliFactor(pout);
    if(!m_pauli || is_hydrogen || is_free_neutron) return 1;
    return m_sea.PauliFactor(p, in, pout, out);
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
    const double pauli = PauliFactor(pIn.P(), had_in[0].ID(), pOut.P(), had_out[0].ID());
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
        auto current = OneBodyCurrent(ubar, u, qVec, ffVal);
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
