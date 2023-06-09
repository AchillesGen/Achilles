#ifndef NUCLEAR_MODEL_HH
#define NUCLEAR_MODEL_HH

#include "Achilles/Current.hh"
#include "Achilles/Event.hh"
#include "Achilles/Factory.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/SpectralFunction.hh"

#include "yaml-cpp/node/node.h"

#include <complex>
#include <stdexcept>

namespace achilles {

class PID;
class PSBuilder;
class Spinor;

enum class NuclearMode {
    None = -1,
    Coherent,
    Quasielastic,
    MesonExchangeCurrent,
    Interference_QE_MEC,
    Resonance,
    ShallowInelastic,
    DeepInelastic
};

inline std::string ToString(NuclearMode mode) {
    switch(mode) {
    case NuclearMode::None:
        return "None";
    case NuclearMode::Coherent:
        return "Coherent";
    case NuclearMode::Quasielastic:
        return "Quasielastic";
    case NuclearMode::MesonExchangeCurrent:
        return "MesonExchangeCurrent";
    case NuclearMode::Interference_QE_MEC:
        return "Interference_QE_MEC";
    case NuclearMode::Resonance:
        return "Resonance";
    case NuclearMode::ShallowInelastic:
        return "ShallowInelastic";
    case NuclearMode::DeepInelastic:
        return "DeepInelastic";
    }

    throw std::runtime_error("Invalid NuclearMode");
}

enum class WardGauge {
    None = -1,
    Coulomb,
    Weyl,
    Landau,
};

inline WardGauge ToEnum(const std::string &ward) {
    if(ward == "None")
        return WardGauge::None;
    else if(ward == "Coulomb")
        return WardGauge::Coulomb;
    else if(ward == "Weyl")
        return WardGauge::Weyl;
    else if(ward == "Landau")
        return WardGauge::Landau;
    auto msg = fmt::format("WardGauge: Invalid option {}", ward);
    throw std::runtime_error(msg);
}

class NuclearModel {
  public:
    using Current = std::vector<VCurrent>;
    using Currents = std::map<int, Current>;
    using FFInfoMap = std::map<int, std::vector<FormFactorInfo>>;
    using FormFactorMap = std::map<FormFactorInfo::Type, std::complex<double>>;
    using ModelMap = std::unordered_map<NuclearMode, std::unique_ptr<NuclearModel>>;

    NuclearModel() = default;
    NuclearModel(const YAML::Node &, FormFactorBuilder &);
    NuclearModel(const NuclearModel &) = delete;
    NuclearModel(NuclearModel &&) = default;
    NuclearModel &operator=(const NuclearModel &) = delete;
    NuclearModel &operator=(NuclearModel &&) = default;
    virtual ~NuclearModel() = default;

    virtual NuclearMode Mode() const = 0;
    virtual std::string PhaseSpace() const = 0;
    virtual Currents CalcCurrents(const std::vector<FourVector> &, const std::vector<FourVector> &,
                                  const FourVector &, const FFInfoMap &) const = 0;
    virtual std::vector<ProcessInfo> AllowedStates(const ProcessInfo &) const = 0;
    virtual size_t NSpins() const = 0;
    virtual double InitialStateWeight(const std::vector<PID> &,
                                      const std::vector<FourVector> &) const = 0;

    virtual std::string GetName() const = 0;
    static std::string Name() { return "Nuclear Model"; }

  protected:
    FormFactor::Values EvalFormFactor(double q2) const { return m_form_factor->operator()(q2); }
    FormFactorMap CouplingsFF(const FormFactor::Values &,
                              const std::vector<FormFactorInfo> &) const;
    static YAML::Node LoadFormFactor(const YAML::Node &);

  private:
    std::unique_ptr<FormFactor> m_form_factor{nullptr};
};

template <typename Derived>
using RegistrableNuclearModel = Registrable<NuclearModel, Derived, const YAML::Node &>;
using NuclearModelFactory = Factory<NuclearModel, const YAML::Node &>;

NuclearModel::ModelMap LoadModels(const YAML::Node &);

class Coherent : public NuclearModel, RegistrableNuclearModel<Coherent> {
  public:
    Coherent(const YAML::Node &, const YAML::Node &, FormFactorBuilder &);

    NuclearMode Mode() const override { return NuclearMode::Coherent; }
    std::string PhaseSpace() const override { return Name(); }
    Currents CalcCurrents(const std::vector<FourVector> &, const std::vector<FourVector> &,
                          const FourVector &, const FFInfoMap &) const override;
    std::vector<ProcessInfo> AllowedStates(const ProcessInfo &) const override;
    size_t NSpins() const override { return 1; }
    double InitialStateWeight(const std::vector<PID> &,
                              const std::vector<FourVector> &) const override {
        return 1;
    }
    std::string GetName() const override { return Coherent::Name(); }

    // Required factory methods
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &);
    static std::string Name() { return "Coherent"; }

  private:
    PID nucleus_pid;
};

class QESpectral : public NuclearModel, RegistrableNuclearModel<QESpectral> {
  public:
    QESpectral(const YAML::Node &, const YAML::Node &, FormFactorBuilder &);

    NuclearMode Mode() const override { return NuclearMode::Quasielastic; }
    std::string PhaseSpace() const override { return Name(); }
    Currents CalcCurrents(const std::vector<FourVector> &, const std::vector<FourVector> &,
                          const FourVector &, const FFInfoMap &) const override;
    std::vector<ProcessInfo> AllowedStates(const ProcessInfo &) const override;
    size_t NSpins() const override { return 4; }
    double InitialStateWeight(const std::vector<PID> &,
                              const std::vector<FourVector> &) const override;
    std::string GetName() const override { return QESpectral::Name(); }

    // Required factory methods
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &);
    static std::string Name() { return "QESpectral"; }

  private:
    void CoulombGauge(VCurrent &, const FourVector &, double) const;
    void WeylGauge(VCurrent &, const FourVector &, double) const;
    void LandauGauge(VCurrent &, const FourVector &) const;
    Current HadronicCurrent(const std::array<Spinor, 2> &, const std::array<Spinor, 2> &,
                            const FourVector &, const FormFactorMap &) const;

    const WardGauge m_ward;
    SpectralFunction spectral_proton, spectral_neutron;
};

} // namespace achilles

#endif
