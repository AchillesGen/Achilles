#pragma once

#include "Achilles/Event.hh"
#include "Achilles/Factory.hh"
#include "Achilles/Integrand.hh"
#include "Achilles/LeptonicCurrent.hh"

namespace YAML {
class Node;
}

namespace achilles {

class SherpaInterface;
class Process;

class XSecBackend {
  public:
    XSecBackend() {}
    virtual ~XSecBackend() = default;
    virtual double CrossSection(const Event &event, const Process &process) const = 0;
    virtual void SetOptions(const YAML::Node &) {}
    virtual void SetSherpa(SherpaInterface *) {}
    virtual void AddNuclearModel(std::unique_ptr<NuclearModel> model) {
        m_model = std::move(model);
    }
    virtual NuclearModel *GetNuclearModel() { return m_model.get(); }
    virtual void AddProcess(Process &process) = 0;
    virtual bool Validate() { return m_model != nullptr; }
    virtual void SetupChannels(const ProcessInfo &, std::shared_ptr<Beam>,
                               Integrand<FourVector> &) = 0;

  protected:
    double FluxFactor(const FourVector &, const FourVector &, const ProcessInfo &) const;
    double InitialStateFactor(size_t, size_t, const std::vector<FourVector> &,
                              const ProcessInfo &) const;
    double SpinAvg(const ProcessInfo &) const;

    std::unique_ptr<NuclearModel> m_model = nullptr;
};

template <typename Derived> using RegistrableBackend = Registrable<XSecBackend, Derived>;
using XSecBackendFactory = Factory<XSecBackend>;

class DefaultBackend : public XSecBackend, RegistrableBackend<DefaultBackend> {
  public:
    DefaultBackend();
    double CrossSection(const Event &event, const Process &process) const override;
    void AddProcess(Process &process) override;
    void SetupChannels(const ProcessInfo &, std::shared_ptr<Beam>,
                       Integrand<FourVector> &) override;

    // Required factory methods
    static std::unique_ptr<XSecBackend> Construct() { return std::make_unique<DefaultBackend>(); }
    static std::string Name() { return "Default"; }

  private:
    std::map<std::pair<PID, PID>, LeptonicCurrent> m_currents;
    FFDictionary form_factors;
};

#ifdef ACHILLES_SHERPA_INTERFACE
class BSMBackend : public XSecBackend, RegistrableBackend<BSMBackend> {
  public:
    BSMBackend();
    void SetSherpa(SherpaInterface *sherpa) override { p_sherpa = sherpa; }
    double CrossSection(const Event &event, const Process &process) const override;
    void AddProcess(Process &process) override;
    void SetupChannels(const ProcessInfo &, std::shared_ptr<Beam>,
                       Integrand<FourVector> &) override;

    // Required factory methods
    static std::unique_ptr<XSecBackend> Construct() { return std::make_unique<BSMBackend>(); }
    static std::string Name() { return "BSM"; }

  private:
    Currents CalcLeptonCurrents(const std::vector<FourVector> &, const ProcessInfo &) const;
    SherpaInterface *p_sherpa = nullptr;
};

class SherpaBackend : public XSecBackend, RegistrableBackend<SherpaBackend> {
  public:
    SherpaBackend();
    void SetSherpa(SherpaInterface *sherpa) override { p_sherpa = sherpa; }
    void SetOptions(const YAML::Node &) override;
    double CrossSection(const Event &event, const Process &process) const override;
    void AddProcess(Process &process) override;
    void SetupChannels(const ProcessInfo &, std::shared_ptr<Beam>,
                       Integrand<FourVector> &) override;

    // Required factory methods
    static std::unique_ptr<XSecBackend> Construct() { return std::make_unique<SherpaBackend>(); }
    static std::string Name() { return "Sherpa"; }

  private:
    SherpaInterface *p_sherpa = nullptr;
};
#endif

class XSecBuilder {
  public:
    XSecBuilder(const std::string &);
    XSecBuilder &AddOptions(const YAML::Node &node);
    XSecBuilder &AddProcess(Process &process);
    XSecBuilder &AddSherpa(SherpaInterface *sherpa);
    XSecBuilder &AddNuclearModel(std::unique_ptr<NuclearModel> model);
    std::unique_ptr<XSecBackend> build();

  private:
    std::unique_ptr<XSecBackend> m_backend = nullptr;
};

} // namespace achilles
