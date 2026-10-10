// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "Achilles/Factory.hh"

#include "yaml-cpp/node/node.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace achilles {

class Event;
class FormFactor;
class NuclearModel;
class Process;
class SpectralFunction;
class XSecBackend;

/// Everything a variation may need to evaluate the weight of an accepted event.
/// The event holds the lab-frame phase space point with its phase space weight,
/// i.e. exactly what was handed to XSecBackend::CrossSection for the nominal weight.
struct VariationContext {
    const Event &event;
    const Process &process;
    XSecBackend &backend;
    double nominal;
};

/// A single alternative weight. Implementations return w_var / w_nominal.
/// Variations must not consume random numbers, so that enabling them does not
/// change the generated event sample.
class Variation {
  public:
    Variation() = default;
    Variation(const Variation &) = delete;
    Variation &operator=(const Variation &) = delete;
    virtual ~Variation() = default;

    virtual double Ratio(const VariationContext &) = 0;
    /// Called once per nuclear model before generation to validate the variation
    virtual void Initialize(const NuclearModel &) {}
};

/// The generation settings themselves, i.e. a ratio of exactly one
class NominalVariation : public Variation {
  public:
    double Ratio(const VariationContext &) override { return 1; }
};

/// Variations that change an ingredient of the calculation, recompute the full
/// cross section with it, and restore the nominal state afterwards.
class RecomputeVariation : public Variation {
  public:
    double Ratio(const VariationContext &) override;

  protected:
    /// Return false if the variation does not change this backend (ratio = 1)
    virtual bool Apply(XSecBackend &) = 0;
    virtual void Restore(XSecBackend &) = 0;
};

/// Rebuild the nucleon form factors from the nominal configuration of each model with
/// edits applied: overrides merged in (parameters or functional forms), a complete
/// replacement configuration, and/or a single parameter set to a value.
class FormFactorVariation : public RecomputeVariation {
  public:
    struct Edit {
        YAML::Node overrides{}, replacement{};
        std::string parameter{};
        double value{};
    };
    explicit FormFactorVariation(Edit edit) : m_edit{std::move(edit)} {}
    void Initialize(const NuclearModel &) override;

    /// Merge overrides into a base form factor config
    static YAML::Node Merge(const YAML::Node &base, const YAML::Node &overrides);
    /// Access a parameter by a '/' separated path, where integers index into sequences
    static YAML::Node Lookup(const YAML::Node &config, const std::string &path);
    static YAML::Node SetParameter(const YAML::Node &config, const std::string &path, double value);
    /// The nominal config with the edit applied, or a null node if the edit changes nothing
    static YAML::Node EditedConfig(const YAML::Node &nominal, const Edit &edit);

  protected:
    bool Apply(XSecBackend &) override;
    void Restore(XSecBackend &) override;

  private:
    std::shared_ptr<FormFactor> Alternate(const NuclearModel &);

    Edit m_edit;
    // A null entry means the edit leaves this model's form factor unchanged
    std::map<const NuclearModel *, std::shared_ptr<FormFactor>> m_cache;
    std::shared_ptr<FormFactor> m_saved;
};

/// Replace the spectral function of single-nucleon knockout models. Since the
/// phase space is sampled independently of S(p, E) and the cross section is linear
/// in the initial state weight, the ratio is exactly S'(p, E) / S(p, E).
class SpectralFunctionVariation : public Variation {
  public:
    SpectralFunctionVariation(const std::string &proton, const std::string &neutron);
    double Ratio(const VariationContext &) override;

  private:
    std::shared_ptr<SpectralFunction> m_proton, m_neutron;
};

/// Uncertainty band from the central value and the member values of a group, for any
/// observable (total cross section, histogram bin, ...).
struct Band {
    double central{}, lower{}, upper{};
};

class Combination {
  public:
    explicit Combination(const YAML::Node &);
    Combination(const Combination &) = delete;
    Combination &operator=(const Combination &) = delete;
    virtual ~Combination() = default;

    Band Combine(double central, const std::vector<double> &members) const;
    virtual void Validate(const std::string &group, size_t nmembers) const;
    virtual std::string Type() const = 0;
    /// Factor applied to the deviations, e.g. to convert a 90% CL band to 68% CL
    double Scale() const { return m_scale; }

    static std::string Name() { return "Variation Combination"; }
    /// Accepts either a string (the type) or a map with 'Type' and optional 'Scale'
    static std::unique_ptr<Combination> Build(const YAML::Node &, const std::string &default_type);

  protected:
    /// Return the (lower, upper) deviations from central, both >= 0
    virtual std::pair<double, double> Deviations(double central,
                                                 const std::vector<double> &members) const = 0;

  private:
    double m_scale{1};
};

template <typename Derived>
using RegistrableCombination = Registrable<Combination, Derived, const YAML::Node &>;
using CombinationFactory = Factory<Combination, const YAML::Node &>;

/// Members are listed for reference only, no band is defined (e.g. parameter scans for fits)
class NoCombination : public Combination, RegistrableCombination<NoCombination> {
  public:
    using Combination::Combination;
    std::string Type() const override { return Name(); }
    void Validate(const std::string &, size_t) const override {}
    static std::unique_ptr<Combination> Construct(const YAML::Node &node) {
        return std::make_unique<NoCombination>(node);
    }
    static std::string Name() { return "None"; }

  protected:
    std::pair<double, double> Deviations(double, const std::vector<double> &) const override {
        return {0, 0};
    }
};

/// Band spanned by the central value and all members
class EnvelopeCombination : public Combination, RegistrableCombination<EnvelopeCombination> {
  public:
    using Combination::Combination;
    std::string Type() const override { return Name(); }
    static std::unique_ptr<Combination> Construct(const YAML::Node &node) {
        return std::make_unique<EnvelopeCombination>(node);
    }
    static std::string Name() { return "Envelope"; }

  protected:
    std::pair<double, double> Deviations(double, const std::vector<double> &) const override;
};

/// Symmetric Hessian: one member per eigenvector, delta = sqrt(sum_i (X_i - X_0)^2)
class SymmetricHessianCombination : public Combination,
                                    RegistrableCombination<SymmetricHessianCombination> {
  public:
    using Combination::Combination;
    std::string Type() const override { return Name(); }
    static std::unique_ptr<Combination> Construct(const YAML::Node &node) {
        return std::make_unique<SymmetricHessianCombination>(node);
    }
    static std::string Name() { return "SymmetricHessian"; }

  protected:
    std::pair<double, double> Deviations(double, const std::vector<double> &) const override;
};

/// Asymmetric Hessian: members are (+, -) pairs per eigenvector,
/// delta_+ = sqrt(sum_i max(X_i+ - X_0, X_i- - X_0, 0)^2) and analogously for delta_-
class HessianCombination : public Combination, RegistrableCombination<HessianCombination> {
  public:
    using Combination::Combination;
    std::string Type() const override { return Name(); }
    void Validate(const std::string &, size_t) const override;
    static std::unique_ptr<Combination> Construct(const YAML::Node &node) {
        return std::make_unique<HessianCombination>(node);
    }
    static std::string Name() { return "Hessian"; }

  protected:
    std::pair<double, double> Deviations(double, const std::vector<double> &) const override;
};

/// Monte Carlo replicas: symmetric band of one standard deviation of the members
class ReplicasCombination : public Combination, RegistrableCombination<ReplicasCombination> {
  public:
    using Combination::Combination;
    std::string Type() const override { return Name(); }
    void Validate(const std::string &, size_t) const override;
    static std::unique_ptr<Combination> Construct(const YAML::Node &node) {
        return std::make_unique<ReplicasCombination>(node);
    }
    static std::string Name() { return "Replicas"; }

  protected:
    std::pair<double, double> Deviations(double, const std::vector<double> &) const override;
};

/// Metadata describing a group, written to the event file header
struct VariationGroupInfo {
    std::string name, type, combination;
    double combination_scale{1};
    std::vector<std::string> weight_names; // central first, then the members
    std::vector<std::string> member_labels;
    std::vector<double> member_values; // parameter values, empty for alternatives
    std::string parameter{};
};

/// One entry of the 'Variations' block: a central value and its members.
/// Weights are named '<group>:central' and '<group>:<member label>'.
class VariationGroup {
  public:
    explicit VariationGroup(const YAML::Node &);
    VariationGroup(const VariationGroup &) = delete;
    VariationGroup &operator=(const VariationGroup &) = delete;
    virtual ~VariationGroup() = default;

    const std::string &GroupName() const { return m_name; }
    virtual std::string Type() const = 0;
    size_t NWeights() const { return 1 + m_members.size(); }
    std::vector<std::string> WeightNames() const;
    VariationGroupInfo Info() const;
    const Combination &GetCombination() const { return *m_combination; }

    void Initialize(const NuclearModel &);
    /// Ratios for the central value followed by each member
    void Evaluate(const VariationContext &, std::vector<double> &) const;

    static std::string Name() { return "Variation Group"; }

  protected:
    void SetCentral(std::unique_ptr<Variation> central) { m_central = std::move(central); }
    void AddMember(std::string label, std::unique_ptr<Variation> member, double value = 0);
    /// Read member values from 'Values', 'Scan: {Min, Max, Steps}' or 'MinMax'.
    /// Returns (label, value) pairs and sets the default combination accordingly
    std::vector<std::pair<std::string, double>> ParameterValues(const YAML::Node &);
    void SetupCombination(const YAML::Node &);
    std::string m_parameter{};

  private:
    std::string m_name;
    std::string m_default_combination{"Envelope"};
    std::unique_ptr<Variation> m_central;
    std::vector<std::string> m_labels;
    std::vector<double> m_values;
    std::vector<std::unique_ptr<Variation>> m_members;
    std::unique_ptr<Combination> m_combination;
};

template <typename Derived>
using RegistrableVariationGroup = Registrable<VariationGroup, Derived, const YAML::Node &>;
using VariationGroupFactory = Factory<VariationGroup, const YAML::Node &>;

/// Form factor parameter scans / min-max pairs, or alternative functional forms:
///
///   - Name: MA
///     Type: FormFactor
///     Parameter: AxialDipole/MA
///     Central: 1.0            # optional, defaults to the value used for generation
///     MinMax: [0.9, 1.1]      # or Values: [...] or Scan: {Min: 0.9, Max: 1.1, Steps: 5}
///   - Name: VectorFF
///     Type: FormFactor
///     Alternatives:
///       - {Name: BBBA, Overrides: {vector: BBBA}}
///       - {Name: Alt, File: alt_formfactors.yml}
class FormFactorGroup : public VariationGroup, RegistrableVariationGroup<FormFactorGroup> {
  public:
    explicit FormFactorGroup(const YAML::Node &);
    std::string Type() const override { return Name(); }
    static std::unique_ptr<VariationGroup> Construct(const YAML::Node &node) {
        return std::make_unique<FormFactorGroup>(node);
    }
    static std::string Name() { return "FormFactor"; }

  private:
    static FormFactorVariation::Edit ParseEdit(const YAML::Node &);
};

/// Alternative spectral functions for single-nucleon knockout models:
///
///   - Name: SF
///     Type: SpectralFunction
///     Central: {SpectralP: ..., SpectralN: ...}   # optional, defaults to the generation one
///     Alternatives:
///       - {Name: MF, SpectralP: ..._MF.data, SpectralN: ..._MF.data}
class SpectralFunctionGroup : public VariationGroup,
                              RegistrableVariationGroup<SpectralFunctionGroup> {
  public:
    explicit SpectralFunctionGroup(const YAML::Node &);
    std::string Type() const override { return Name(); }
    static std::unique_ptr<VariationGroup> Construct(const YAML::Node &node) {
        return std::make_unique<SpectralFunctionGroup>(node);
    }
    static std::string Name() { return "SpectralFunction"; }
};

/// Owns all requested variation groups and evaluates them for accepted events.
class VariationHandler {
  public:
    VariationHandler() = default;
    explicit VariationHandler(const YAML::Node &groups);

    void AddGroup(std::unique_ptr<VariationGroup> group);
    const std::vector<std::unique_ptr<VariationGroup>> &Groups() const { return m_groups; }
    std::vector<std::string> WeightNames() const;
    std::vector<VariationGroupInfo> Info() const;
    size_t NWeights() const { return m_nweights; }

    void Initialize(const NuclearModel &) const;
    /// Ratios w_var / w_nominal for every weight, in the order of WeightNames()
    std::vector<double> Evaluate(const VariationContext &) const;

  private:
    std::vector<std::unique_ptr<VariationGroup>> m_groups;
    size_t m_nweights{};
};

} // namespace achilles
