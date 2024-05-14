#ifndef NUCLEAR_MODEL_HH
#define NUCLEAR_MODEL_HH

#include "Achilles/Current.hh"
#include "Achilles/Event.hh"
#include "Achilles/Factory.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/SpectralFunction.hh"

#include "yaml-cpp/node/node.h"

#include <complex>
#include <stdexcept>

namespace achilles {

class PID;
class PSBuilder;
class Spinor;
class Process;

enum class NuclearMode : int {
    None = -1,
    Coherent = 1,
    Quasielastic = 2,
    MesonExchangeCurrent = 3,
    Interference_QE_MEC = 7,
    Resonance = 4,
    ShallowInelastic = 5,
    DeepInelastic = 6,
    Hyperon = 8
};

enum class NuclearFrame : int {
    Lab = 0,
    QZ,
    Custom,
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
    case NuclearMode::Hyperon:
        return "Hyperon";
    }

    throw std::runtime_error("Invalid NuclearMode");
}

inline int ToID(const NuclearMode mode) {
    if(mode == NuclearMode::None) throw std::runtime_error("Invalid nuclear mode");
    return static_cast<int>(mode) * 100;
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
    virtual std::string PhaseSpace(PID) const = 0;
    virtual Currents CalcCurrents(const std::vector<Particle> &, const std::vector<Particle> &, 
        const std::vector<Particle> &, const FourVector &, const FFInfoMap &) const = 0;
    virtual std::vector<ProcessInfo> AllowedStates(const ProcessInfo &) const;
    virtual size_t NSpins() const;
    virtual double InitialStateWeight(const std::vector<Particle> &, const std::vector<Particle> &, size_t, size_t) const = 0;

    virtual std::string GetName() const = 0;
    static std::string Name() { return "Nuclear Model"; }
    virtual std::string PSName() const = 0;
    virtual std::string InspireHEP() const = 0;
    virtual NuclearFrame Frame() const { return NuclearFrame::Lab; }
    void TransformFrame(Event &, const Process &, bool) const;
    void SetTransform();

  protected:
    FormFactor::Values EvalFormFactor(double q2) const { return m_form_factor->operator()(q2); }
    FormFactorMap CouplingsFF(const FormFactor::Values &,
                              const std::vector<FormFactorInfo> &) const;
    std::function<void(Event &, const Process &, bool)> transform;
    void TransformLab(Event &, const Process &, bool) const {}
    void TransformQZ(Event &, const Process &, bool);
    virtual void TransformCustom(Event &, const Process &, bool) const {}
    static YAML::Node LoadFormFactor(const YAML::Node &);
    static YAML::Node LoadModelParams(const YAML::Node &);

    void CoulombGauge(VCurrent &, const FourVector &, double) const;
    void WeylGauge(VCurrent &, const FourVector &, double) const;
    void LandauGauge(VCurrent &, const FourVector &) const;

  private:
    FourVector::RotMat rotation;
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
    std::string PhaseSpace(PID) const override;
    Currents CalcCurrents(const std::vector<Particle> &, const std::vector<Particle> &,
            const std::vector<Particle> &,const FourVector &, const FFInfoMap &) const override;
    std::vector<ProcessInfo> AllowedStates(const ProcessInfo &) const override;
    size_t NSpins() const override { return 1; }
    double InitialStateWeight(const std::vector<Particle> &, const std::vector<Particle> &, size_t, size_t) const override {
        return 1;
    }
    std::string GetName() const override { return Coherent::Name(); }
    std::string InspireHEP() const override { return ""; }
    std::string PSName() const override { return "Coherent"; }

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
    std::string PhaseSpace(PID) const override;
    Currents CalcCurrents(const std::vector<Particle> &, const std::vector<Particle> &,
         const std::vector<Particle> &,const FourVector &, const FFInfoMap &) const override;
    size_t NSpins() const override { return 4; }
    double InitialStateWeight(const std::vector<Particle> &, const std::vector<Particle> &, size_t, size_t) const override;
    std::string GetName() const override { return QESpectral::Name(); }
    std::string InspireHEP() const override { return ""; }
    std::string PSName() const override { return "OneBodySpectral"; }

    // Required factory methods
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &);
    static std::string Name() { return "QESpectral"; }

  private:
    Current HadronicCurrent(const std::array<Spinor, 2> &, const std::array<Spinor, 2> &,
                            const FourVector &, const FormFactorMap &) const;

    const WardGauge m_ward;
    SpectralFunction spectral_proton, spectral_neutron;
    // TODO: This is a code smell. Should figure out a better solution
    mutable bool is_hydrogen{false};
};

class HyperonSpectral : public NuclearModel, RegistrableNuclearModel<HyperonSpectral> {
  public:
    HyperonSpectral(const YAML::Node &, const YAML::Node &, FormFactorBuilder &);

    NuclearMode Mode() const override { return NuclearMode::Hyperon; }
    std::string PhaseSpace(PID) const override;
    Currents CalcCurrents(const std::vector<Particle> &, const std::vector<Particle> &,
         const std::vector<Particle> &,const FourVector &, const FFInfoMap &) const override;
    size_t NSpins() const override { return 4; }
    double InitialStateWeight(const std::vector<Particle> &, const std::vector<Particle> &, size_t, size_t) const override;
    std::string GetName() const override { return HyperonSpectral::Name(); }
    std::string InspireHEP() const override { return ""; }
    std::string PSName() const override { return "OneBodySpectral"; }

    // Required factory methods
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &);
    static std::string Name() { return "HyperonSpectral"; }

  private:
    Current HadronicCurrent(const std::array<Spinor, 2> &, const std::array<Spinor, 2> &,
                            const FourVector &, const FormFactorMap &, PID) const;

    const WardGauge m_ward;
    SpectralFunction spectral_proton, spectral_neutron;
    // TODO: This is a code smell. Should figure out a better solution
    mutable bool is_hydrogen{false};
};

} // namespace achilles

#endif
