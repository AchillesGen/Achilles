// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef HADRONIC_MAPPER_HH
#define HADRONIC_MAPPER_HH

#include <cmath>

#include "Achilles/Factory.hh"
#include "Achilles/Mapper.hh"
#include "Achilles/ProcessInfo.hh"

#include "yaml-cpp/node/node.h"

namespace achilles {

class FourVector;

class HadronicBeamMapper : public Mapper<FourVector> {
  public:
    HadronicBeamMapper(const ProcessInfo &info, size_t idx) : m_info{info}, m_idx{idx} {}

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override = 0;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override = 0;
    size_t NDims() const override = 0;
    // Model-specific settings (NuclearModel::PhaseSpaceSettings); ignored by default
    virtual void Configure(const YAML::Node &) {}
    static std::string Name() { return "Hadronic Initial State"; }

  protected:
    size_t HadronIdx() const { return m_idx; }
    ProcessInfo m_info;

  private:
    size_t m_idx;
};

class QESpectralMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, QESpectralMapper, const ProcessInfo &, size_t> {
  public:
    QESpectralMapper(const ProcessInfo &info, size_t idx);
    static std::string Name() { return "OneBodySpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info, size_t idx) {
        return std::make_unique<QESpectralMapper>(info, idx);
    }

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return 4; }

  private:
    // static constexpr double dCos = 2;
    static constexpr double dPhi = 2 * M_PI;
    // static constexpr double dp = 800;
    // static constexpr double dE = 400;
};

// Initial-state mapper for models with an analytic momentum distribution (RFG, CFG): the bound
// nucleon has a fixed energy for a given momentum, E_in(p) = sqrt(p^2 + mN^2) - Eb, so only
// (|p|, cos, phi) are sampled. The weight is the density in d^3p (no dE integral).
class MomentumMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, MomentumMapper, const ProcessInfo &, size_t> {
  public:
    MomentumMapper(const ProcessInfo &info, size_t idx);
    static std::string Name() { return "OneBodyMomentum"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info, size_t idx) {
        return std::make_unique<MomentumMapper>(info, idx);
    }

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return 3; }

    // Reads PCut (largest sampled momentum) and Binding from QEFermiGas::PhaseSpaceSettings
    void Configure(const YAML::Node &) override;
    double PCut() const { return m_pcut; }
    double Binding() const { return m_binding; }
    double InitialEnergy(double p) const;

  private:
    struct Limits {
        double pmin, pmax;
        bool valid;
    };
    Limits MomentumLimits(const FourVector &) const;
    double CosMax(const FourVector &, double) const;
    double SMax(const FourVector &, double) const;
    static constexpr double dPhi = 2 * M_PI;
    double m_pcut{250}, m_binding{0};
};

class CoherentMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, CoherentMapper, const ProcessInfo &, size_t> {
  public:
    CoherentMapper(const ProcessInfo &info, size_t idx);
    static std::string Name() { return "Coherent"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info, size_t idx) {
        return std::make_unique<CoherentMapper>(info, idx);
    }

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return 0; }

  private:
    double m_mass;
};

class IntfSpectralMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, IntfSpectralMapper, const ProcessInfo &, size_t> {
  public:
    IntfSpectralMapper(const ProcessInfo &info, size_t idx);
    static std::string Name() { return "IntfSpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info, size_t idx) {
        return std::make_unique<IntfSpectralMapper>(info, idx);
    }

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return 7; }

  private:
    static constexpr double dCos2 = 2;
    static constexpr double dPhi = 2 * M_PI;
    static constexpr double dp2 = 400; // reduced b/c MF SF has no strength at high P
    // static constexpr double dE = 400;
};

} // namespace achilles

#endif
