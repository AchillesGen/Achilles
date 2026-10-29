// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef QE_FERMI_GAS_HH
#define QE_FERMI_GAS_HH

#include "Achilles/NuclearModel.hh"

#include <vector>

namespace achilles {

// Local Fermi sea on radial shells: kF_N(r) = hbar c (3 pi^2 rho_N(r))^{1/3}, N = p or n.
// A nucleon of type N and momentum p can sit wherever kF_N(r) > p, a volume V_N(p). The number
// of such nucleons per d^3p is 2 V_N(p) / (2 pi hbar c)^3. Its final state (type N', momentum
// p') is free where kF_N'(r) <= p' at the same r, so the Pauli factor is the fraction of V_N(p)
// with kF_N'(r) <= p'; for N' = N this is 1 - V_N(p') / V_N(p).
class LocalFermiSea {
  public:
    LocalFermiSea() = default;
    // Shell edges [fm] (n+1, increasing) and the proton and neutron densities [fm^-3] at the
    // n shell centres
    LocalFermiSea(const std::vector<double> &edges, const std::vector<double> &rho_p,
                  const std::vector<double> &rho_n);

    double Volume(double k, PID nucleon) const; // fm^3 where kF_nucleon(r) > k [MeV]
    double PauliFactor(double p, PID in, double pout, PID out) const;
    double NNucleons(PID nucleon) const; // integral of rho_N d^3r, for checks
    double MaxFermiMomentum() const { return m_kf_max; }

  private:
    const std::vector<double> &FermiMomenta(PID nucleon) const;

    std::vector<double> m_dv, m_kf_p, m_kf_n, m_rho_p, m_rho_n;
    double m_kf_max{};
};

// Quasielastic scattering on a (correlated) relativistic Fermi gas.
//
//   RFG: n(p) = 3/(4 pi kF^3) theta(kF - p),        occupation(k) = theta(kF - k)
//   CFG: n(p) = 3/(4 pi kF^3) n_CFG(p),             occupation(k) = n_CFG(k)
//        n_CFG = alpha0 (k < kF), alpha1/k^4 (kF < k < lambda kF), 0 otherwise
//        alpha0 = 1 - (1 - 1/lambda) c0/pi^2,  alpha1 = c0 kF^4 / (3 pi^2)
//
// Bound nucleon: E_in(p) = sqrt(p^2 + mN^2) - Eb. With the off-shell prescription of
// QESpectral this gives the Smith-Moniz shift  q0 -> omega - Eb  in the hadronic current.
// Weight: A_N * n(|p|) (InitialStateWeight, A_N = Z or N) times the Pauli factor
// [1 - occupation(|p'|)], applied to |J|^2 in CalcCurrents (the only hook that sees p').
//
// kF: Local (RFG only) replaces the global kF by the local Fermi sea of the target nucleus
// (LocalFermiSea, built in SetNucleus from the nucleus' proton and neutron densities):
// weight 2 V_N(p) / (2 pi hbar c)^3, Pauli factor from the same radial shells.
class QEFermiGas : public NuclearModel, RegistrableNuclearModel<QEFermiGas> {
  public:
    enum class MomentumDistribution { RFG, CFG };

    QEFermiGas(const YAML::Node &, const YAML::Node &, FormFactorBuilder &);

    NuclearMode Mode() const override { return NuclearMode::Quasielastic; }
    std::string PhaseSpace(PID) const override;
    Currents CalcCurrents(const std::vector<Particle> &, const std::vector<Particle> &,
                          const std::vector<Particle> &, const FourVector &,
                          const FFInfoMap &) const override;
    size_t NSpins() const override { return 4; }
    double InitialStateWeight(const std::vector<Particle> &, const std::vector<Particle> &, size_t,
                              size_t) const override;
    std::string GetName() const override { return QEFermiGas::Name(); }
    std::string InspireHEP() const override;
    std::string PSName() const override { return "OneBodyMomentum"; }
    YAML::Node PhaseSpaceSettings() const override;
    void SetNucleus(const Nucleus &) override;

    // Required factory methods
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &);
    static std::string Name() { return "QEFermiGas"; }

    // Exposed for testing
    double Density(double p) const;        // normalized: int d^3p Density = 1
    double Occupation(double k) const;     // 0 <= occupation <= 1
    double PauliFactor(double pout) const; // global kF: 1 - occupation(|p'|), or 1 if disabled
    // Pauli factor for a struck nucleon (in, |p|) going to (out, |p'|); local or global kF
    double PauliFactor(double p, PID in, double pout, PID out) const;
    double PCut() const {
        if(m_local) return m_sea.MaxFermiMomentum();
        return m_distribution == MomentumDistribution::CFG ? m_lambda * m_kf : m_kf;
    }
    bool IsLocal() const { return m_local; }
    const LocalFermiSea &Sea() const { return m_sea; }

  private:
    const WardGauge m_ward;
    MomentumDistribution m_distribution{MomentumDistribution::RFG};
    double m_kf{}, m_eb{}, m_c0{}, m_lambda{1};
    double m_alpha0{1}, m_alpha1{0};
    bool m_pauli{true};
    bool m_local{false};
    LocalFermiSea m_sea;
    mutable bool is_hydrogen{false};
    mutable bool is_free_neutron{false};
};

} // namespace achilles

#endif
