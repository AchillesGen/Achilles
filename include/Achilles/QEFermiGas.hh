// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef QE_FERMI_GAS_HH
#define QE_FERMI_GAS_HH

#include "Achilles/NuclearModel.hh"

namespace achilles {

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

    // Required factory methods
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &);
    static std::string Name() { return "QEFermiGas"; }

    // Exposed for testing
    double Density(double p) const;        // normalized: int d^3p Density = 1
    double Occupation(double k) const;     // 0 <= occupation <= 1
    double PauliFactor(double pout) const; // 1 - occupation(|p'|), or 1 if disabled
    double PCut() const {
        return m_distribution == MomentumDistribution::CFG ? m_lambda * m_kf : m_kf;
    }

  private:
    Current HadronicCurrent(const std::array<Spinor, 2> &, const std::array<Spinor, 2> &,
                            const FourVector &, const FormFactorMap &) const;

    const WardGauge m_ward;
    MomentumDistribution m_distribution{MomentumDistribution::RFG};
    double m_kf{}, m_eb{}, m_c0{}, m_lambda{1};
    double m_alpha0{1}, m_alpha1{0};
    bool m_pauli{true};
    mutable bool is_hydrogen{false};
    mutable bool is_free_neutron{false};
};

} // namespace achilles

#endif
