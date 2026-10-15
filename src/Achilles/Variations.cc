// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Variations.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Process.hh"
#include "Achilles/SpectralFunction.hh"
#include "Achilles/System.hh"
#include "Achilles/Utilities.hh"
#include "Achilles/XSecBackend.hh"

#include "yaml-cpp/yaml.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <regex>

using achilles::Combination;
using achilles::EnvelopeCombination;
using achilles::FFType;
using achilles::FormFactorGroup;
using achilles::FormFactorVariation;
using achilles::HessianCombination;
using achilles::RecomputeVariation;
using achilles::ReplicasCombination;
using achilles::SpectralFunctionGroup;
using achilles::SpectralFunctionVariation;
using achilles::SymmetricHessianCombination;
using achilles::VariationGroup;
using achilles::VariationHandler;
using achilles::VariationShorthand;

namespace {

std::string FormatValue(double value) {
    return fmt::format("{:.10g}", value);
}

std::vector<std::string> SplitPath(const std::string &path) {
    std::vector<std::string> tokens;
    achilles::tokenize(path, tokens, "/");
    return tokens;
}

bool IsIndex(const std::string &token) {
    return !token.empty() && std::all_of(token.begin(), token.end(), ::isdigit);
}

} // namespace

// ---------------------------------------------------------------------------
// Individual variations
// ---------------------------------------------------------------------------

double RecomputeVariation::Ratio(const VariationContext &ctx) {
    if(ctx.nominal == 0) return 0;
    if(!Apply(ctx.backend)) return 1;
    double weight = 0;
    try {
        weight = ctx.backend.CrossSection(ctx.event, ctx.process);
    } catch(...) {
        Restore(ctx.backend);
        throw;
    }
    Restore(ctx.backend);
    return weight / ctx.nominal;
}

YAML::Node FormFactorVariation::Merge(const YAML::Node &base, const YAML::Node &overrides) {
    YAML::Node result = YAML::Clone(base);
    for(const auto &entry : overrides) {
        const auto key = entry.first.as<std::string>();
        if(entry.second.IsMap() && result[key] && result[key].IsMap()) {
            result[key] = Merge(result[key], entry.second);
        } else {
            result[key] = YAML::Clone(entry.second);
        }
    }
    return result;
}

YAML::Node FormFactorVariation::Lookup(const YAML::Node &config, const std::string &path) {
    // NOTE: Must traverse a const node, since operator[] on a non-const node inserts keys
    const YAML::Node root = config;
    YAML::Node current = YAML::Clone(root);
    for(const auto &token : SplitPath(path)) {
        YAML::Node next;
        if(current.IsSequence() && IsIndex(token)) {
            const auto idx = std::stoul(token);
            if(idx >= current.size())
                throw std::runtime_error(
                    fmt::format("FormFactorVariation: Index {} out of range in '{}'", idx, path));
            next = current[idx];
        } else if(current.IsMap() && current[token]) {
            next = current[token];
        } else {
            throw std::runtime_error(fmt::format(
                "FormFactorVariation: Could not find '{}' of parameter '{}'", token, path));
        }
        current.reset(next);
    }
    return current;
}

YAML::Node FormFactorVariation::SetParameter(const YAML::Node &config, const std::string &path,
                                             double value) {
    if(!Lookup(config, path).IsScalar())
        throw std::runtime_error(
            fmt::format("FormFactorVariation: Parameter '{}' is not a single number", path));

    YAML::Node result = YAML::Clone(config);
    YAML::Node current = result;
    for(const auto &token : SplitPath(path)) {
        YAML::Node next = current.IsSequence() ? current[std::stoul(token)] : current[token];
        // NOTE: reset rebinds the handle, assignment would overwrite the referenced node
        current.reset(next);
    }
    current = value;
    return result;
}

YAML::Node FormFactorVariation::EditedConfig(const YAML::Node &nominal, const Edit &edit) {
    // NOTE: Assigning to a YAML::Node overwrites the node it refers to, so every step works on a
    // fresh clone to never modify the nominal configuration of the model
    const bool only_parameter = edit.replacement.IsNull() && edit.overrides.IsNull();
    const YAML::Node base = YAML::Clone(edit.replacement.IsNull() ? nominal : edit.replacement);
    const YAML::Node merged = edit.overrides.IsNull() ? base : Merge(base, edit.overrides);
    if(edit.parameters.empty()) return merged;

    YAML::Node result = YAML::Clone(merged);
    bool unchanged = only_parameter;
    for(const auto &parameter : edit.parameters) {
        const double nominal_value = Lookup(merged, parameter.path).as<double>();
        const double value = parameter.shift ? nominal_value + parameter.value : parameter.value;
        unchanged = unchanged && value == nominal_value;
        // NOTE: reset rebinds the handle, assignment would write through to the node it refers to
        result.reset(SetParameter(result, parameter.path, value));
    }
    // Setting parameters to the values used for generation leaves the model unchanged
    if(unchanged) return {};
    return result;
}

std::shared_ptr<achilles::FormFactor> FormFactorVariation::Alternate(const NuclearModel &model) {
    auto it = m_cache.find(&model);
    if(it != m_cache.end()) return it->second;

    std::shared_ptr<FormFactor> form_factor = nullptr;
    const auto config = EditedConfig(model.FormFactorConfig(), m_edit);
    if(!config.IsNull()) {
        FormFactorBuilder builder;
        form_factor = NuclearModel::BuildFormFactor(config, builder);
    }
    m_cache[&model] = form_factor;
    return form_factor;
}

void FormFactorVariation::Initialize(const NuclearModel &model) {
    if(model.FormFactorConfig().IsNull()) return;
    Alternate(model);
}

bool FormFactorVariation::Apply(XSecBackend &backend) {
    auto *model = backend.GetNuclearModel();
    if(!model || model->FormFactorConfig().IsNull()) return false;
    auto alternate = Alternate(*model);
    if(!alternate) return false;
    m_saved = model->SwapFormFactor(std::move(alternate));
    return true;
}

void FormFactorVariation::Restore(XSecBackend &backend) {
    backend.GetNuclearModel()->SwapFormFactor(std::move(m_saved));
}

SpectralFunctionVariation::SpectralFunctionVariation(const std::string &proton,
                                                     const std::string &neutron)
    : m_proton{std::make_shared<SpectralFunction>(proton)},
      m_neutron{std::make_shared<SpectralFunction>(neutron)} {}

double SpectralFunctionVariation::Ratio(const VariationContext &ctx) {
    const auto *model = ctx.backend.GetNuclearModel();
    Particle lepton_in;
    std::vector<Particle> hadron_in, lepton_out, hadron_out, spect;
    ctx.process.ExtractParticles(ctx.event, lepton_in, hadron_in, lepton_out, hadron_out, spect);

    // Only single nucleon knockout models are linear in a single S(p, E)
    if(hadron_in.size() != 1 || !spect.empty()) return 1;
    const auto pid = hadron_in[0].ID();
    if(pid != PID::proton() && pid != PID::neutron()) return 1;
    const auto *nominal = model->GetSpectralFunction(pid);
    if(!nominal) return 1;

    const double mom = hadron_in[0].Momentum().P();
    const double removal_energy = Constant::mN - hadron_in[0].E();
    const double denom = (*nominal)(mom, removal_energy);
    if(denom <= 0) return 0;
    const auto &alternate = pid == PID::proton() ? *m_proton : *m_neutron;
    return alternate(mom, removal_energy) / denom;
}

// ---------------------------------------------------------------------------
// Combinations
// ---------------------------------------------------------------------------

Combination::Combination(const YAML::Node &node) {
    if(node.IsMap() && node["Scale"]) m_scale = node["Scale"].as<double>();
}

std::unique_ptr<Combination> Combination::Build(const YAML::Node &node,
                                                const std::string &default_type) {
    std::string type = default_type;
    YAML::Node options;
    if(node && node.IsScalar()) {
        type = node.as<std::string>();
    } else if(node && node.IsMap()) {
        type = node["Type"].as<std::string>();
        options = node;
    }
    try {
        return CombinationFactory::Initialize(type, options);
    } catch(std::out_of_range &) {
        spdlog::error("Variations: Requested combination \"{}\", did you mean \"{}\"", type,
                      GetSuggestion(CombinationFactory::List(), type));
        throw;
    }
}

achilles::Band Combination::Combine(double central, const std::vector<double> &members) const {
    auto [lower, upper] = Deviations(central, members);
    return {central, central - m_scale * lower, central + m_scale * upper};
}

void Combination::Validate(const std::string &group, size_t nmembers) const {
    if(nmembers == 0)
        throw std::runtime_error(
            fmt::format("Variations: Group {} needs members to define a {} band", group, Type()));
}

std::pair<double, double> EnvelopeCombination::Deviations(double central,
                                                          const std::vector<double> &x) const {
    const auto [min, max] = std::minmax_element(x.begin(), x.end());
    return {std::max(0.0, central - *min), std::max(0.0, *max - central)};
}

std::pair<double, double>
SymmetricHessianCombination::Deviations(double central, const std::vector<double> &x) const {
    double sum2 = 0;
    for(const auto &value : x) sum2 += (value - central) * (value - central);
    return {std::sqrt(sum2), std::sqrt(sum2)};
}

void HessianCombination::Validate(const std::string &group, size_t nmembers) const {
    if(nmembers == 0 || nmembers % 2 != 0)
        throw std::runtime_error(fmt::format(
            "Variations: Group {} needs (+, -) pairs of members for a Hessian band, got {}", group,
            nmembers));
}

std::pair<double, double> HessianCombination::Deviations(double central,
                                                         const std::vector<double> &x) const {
    double up2 = 0, down2 = 0;
    for(size_t i = 0; i + 1 < x.size(); i += 2) {
        const double plus = x[i] - central, minus = x[i + 1] - central;
        const double up = std::max({plus, minus, 0.0});
        const double down = std::max({-plus, -minus, 0.0});
        up2 += up * up;
        down2 += down * down;
    }
    return {std::sqrt(down2), std::sqrt(up2)};
}

void ReplicasCombination::Validate(const std::string &group, size_t nmembers) const {
    if(nmembers < 2)
        throw std::runtime_error(fmt::format(
            "Variations: Group {} needs at least two replicas, got {}", group, nmembers));
}

std::pair<double, double> ReplicasCombination::Deviations(double,
                                                          const std::vector<double> &x) const {
    const double mean = std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
    double sum2 = 0;
    for(const auto &value : x) sum2 += (value - mean) * (value - mean);
    const double sigma = std::sqrt(sum2 / static_cast<double>(x.size() - 1));
    return {sigma, sigma};
}

// ---------------------------------------------------------------------------
// Groups
// ---------------------------------------------------------------------------

VariationGroup::VariationGroup(const YAML::Node &node) : m_name{node["Name"].as<std::string>()} {
    if(m_name.empty() || m_name.find_first_of(" \t:|") != std::string::npos)
        throw std::runtime_error(fmt::format(
            "Variations: Invalid group name '{}', must not contain whitespace, ':' or '|'",
            m_name));
    m_central = std::make_unique<NominalVariation>();
}

void VariationGroup::AddMember(std::string label, std::unique_ptr<Variation> member, double value) {
    if(label.empty() || label.find_first_of(" \t|") != std::string::npos || label == "central" ||
       std::find(m_labels.begin(), m_labels.end(), label) != m_labels.end())
        throw std::runtime_error(fmt::format(
            "Variations: Invalid or duplicate member name '{}' in group {}", label, m_name));
    m_labels.push_back(std::move(label));
    m_values.push_back(value);
    m_members.push_back(std::move(member));
}

std::vector<std::pair<std::string, double>>
VariationGroup::ParameterValues(const YAML::Node &node) {
    const int nspecs = (node["Values"] ? 1 : 0) + (node["Scan"] ? 1 : 0) + (node["MinMax"] ? 1 : 0);
    if(nspecs != 1)
        throw std::runtime_error(fmt::format(
            "Variations: Group {} needs exactly one of 'Values', 'Scan' or 'MinMax'", m_name));

    std::vector<std::pair<std::string, double>> result;
    if(node["MinMax"]) {
        const auto minmax = node["MinMax"].as<std::vector<double>>();
        if(minmax.size() != 2 || minmax[0] > minmax[1])
            throw std::runtime_error(
                fmt::format("Variations: Group {} expects 'MinMax: [min, max]'", m_name));
        result = {{"min", minmax[0]}, {"max", minmax[1]}};
        m_default_combination = EnvelopeCombination::Name();
    } else {
        std::vector<double> values;
        if(node["Values"]) {
            values = node["Values"].as<std::vector<double>>();
        } else {
            const auto min = node["Scan"]["Min"].as<double>();
            const auto max = node["Scan"]["Max"].as<double>();
            const auto steps = node["Scan"]["Steps"].as<size_t>();
            if(steps < 2 || min >= max)
                throw std::runtime_error(fmt::format(
                    "Variations: Group {} expects 'Scan: {{Min, Max, Steps >= 2}}'", m_name));
            for(size_t i = 0; i < steps; ++i)
                values.push_back(min + (max - min) * static_cast<double>(i) /
                                           static_cast<double>(steps - 1));
        }
        for(const auto &value : values) result.emplace_back(FormatValue(value), value);
        m_default_combination = NoCombination::Name();
    }
    return result;
}

void VariationGroup::SetupCombination(const YAML::Node &node, double default_scale) {
    YAML::Node combination;
    if(node["Combination"]) {
        combination = YAML::Clone(node["Combination"]);
    } else {
        combination["Type"] = m_default_combination;
        combination["Scale"] = default_scale;
    }
    m_combination = Combination::Build(combination, m_default_combination);
    m_combination->Validate(m_name, m_members.size());
}

std::vector<std::string> VariationGroup::WeightNames() const {
    std::vector<std::string> names{m_name + ":central"};
    for(const auto &label : m_labels) names.push_back(m_name + ":" + label);
    return names;
}

achilles::VariationGroupInfo VariationGroup::Info() const {
    VariationGroupInfo info;
    info.name = m_name;
    info.type = Type();
    info.combination = m_combination->Type();
    info.combination_scale = m_combination->Scale();
    info.weight_names = WeightNames();
    info.member_labels = m_labels;
    if(!m_parameter.empty() || !m_parameters.empty()) {
        info.parameter = m_parameter;
        info.parameters = m_parameters;
        info.member_values = m_values;
    }
    info.vectors = m_info_vectors;
    return info;
}

void VariationGroup::Initialize(const NuclearModel &model) {
    m_central->Initialize(model);
    for(auto &member : m_members) member->Initialize(model);
}

void VariationGroup::Evaluate(const VariationContext &ctx, std::vector<double> &ratios) const {
    ratios.push_back(m_central->Ratio(ctx));
    for(const auto &member : m_members) ratios.push_back(member->Ratio(ctx));
}

FormFactorVariation::Edit FormFactorGroup::ParseEdit(const YAML::Node &node) {
    FormFactorVariation::Edit edit;
    if(node["Overrides"]) edit.overrides = YAML::Clone(node["Overrides"]);
    if(node["File"])
        edit.replacement =
            YAML::LoadFile(Filesystem::FindFile(node["File"].as<std::string>(), "Variations"));
    if(node["Parameter"]) {
        edit.parameters.push_back(
            {node["Parameter"].as<std::string>(), node["Value"].as<double>(), false});
    }
    if(edit.overrides.IsNull() && edit.replacement.IsNull() && edit.parameters.empty())
        throw std::runtime_error(
            "Variations: A form factor variation needs 'Overrides', 'File' or 'Parameter'");
    return edit;
}

FormFactorGroup::FormFactorGroup(const YAML::Node &node) : VariationGroup(node) {
    if(node["Parameter"] && node["Alternatives"])
        throw std::runtime_error(fmt::format(
            "Variations: Group {} has both 'Parameter' and 'Alternatives'", GroupName()));

    if(node["Parameter"]) {
        m_parameter = node["Parameter"].as<std::string>();
        auto edit_for = [&](double value) {
            FormFactorVariation::Edit edit;
            edit.parameters.push_back({m_parameter, value, false});
            return std::make_unique<FormFactorVariation>(edit);
        };
        if(node["Central"]) SetCentral(edit_for(node["Central"].as<double>()));
        for(const auto &[label, value] : ParameterValues(node))
            AddMember(label, edit_for(value), value);
    } else if(node["Alternatives"]) {
        if(node["Central"])
            SetCentral(std::make_unique<FormFactorVariation>(ParseEdit(node["Central"])));
        for(const auto &alternative : node["Alternatives"]) {
            AddMember(alternative["Name"].as<std::string>(),
                      std::make_unique<FormFactorVariation>(ParseEdit(alternative)));
        }
    } else {
        throw std::runtime_error(fmt::format(
            "Variations: Group {} needs either 'Parameter' or 'Alternatives'", GroupName()));
    }
    SetupCombination(node);
}

SpectralFunctionGroup::SpectralFunctionGroup(const YAML::Node &node) : VariationGroup(node) {
    auto build = [](const YAML::Node &spec) {
        return std::make_unique<SpectralFunctionVariation>(spec["SpectralP"].as<std::string>(),
                                                           spec["SpectralN"].as<std::string>());
    };
    if(node["Central"]) SetCentral(build(node["Central"]));
    if(!node["Alternatives"])
        throw std::runtime_error(
            fmt::format("Variations: Group {} needs 'Alternatives'", GroupName()));
    for(const auto &alternative : node["Alternatives"])
        AddMember(alternative["Name"].as<std::string>(), build(alternative));
    SetupCombination(node);
}

// ---------------------------------------------------------------------------
// Flux variations
// ---------------------------------------------------------------------------

double achilles::FluxVariation::Scale(double energy) const {
    if(energy < m_edges.front() || energy >= m_edges.back()) return 1;
    const auto bin = std::upper_bound(m_edges.begin(), m_edges.end(), energy) - m_edges.begin() - 1;
    return m_scales[static_cast<size_t>(bin)];
}

double achilles::FluxVariation::Normalization(const Beam &beam, int pid) {
    const auto key = std::make_pair(&beam, pid);
    auto it = m_norms.find(key);
    if(it != m_norms.end()) return it->second;

    // Flux weighted average of the scale factors over the bins
    static constexpr size_t npoints = 512;
    double total = 0, scaled = 0;
    for(size_t ibin = 0; ibin + 1 < m_edges.size(); ++ibin) {
        const double lo = m_edges[ibin], hi = m_edges[ibin + 1];
        const double step = (hi - lo) / static_cast<double>(npoints);
        double integral = 0;
        for(size_t i = 0; i <= npoints; ++i) {
            const double energy = lo + step * static_cast<double>(i);
            const double weight = (i == 0 || i == npoints) ? 0.5 : 1.0;
            integral += weight * beam.EvaluateFlux(PID(pid), {energy, 0, 0, energy});
        }
        integral *= step;
        total += integral;
        scaled += m_scales[ibin] * integral;
    }
    const double norm = total > 0 ? scaled / total : 1;
    m_norms[key] = norm;
    return norm;
}

double achilles::FluxVariation::Ratio(const VariationContext &ctx) {
    const auto pid = static_cast<int>(ctx.process.Info().m_leptonic.first.AsInt());
    if(m_pid != 0 && pid != m_pid) return 1;
    if(!ctx.beam) throw std::runtime_error("FluxVariation: No beam available");
    // A fixed beam energy has no flux shape, so flux averaged cross sections do not change
    const bool monochromatic = ctx.beam->at(PID(pid))->Type() == "Monochromatic";
    const double scale = Scale(ctx.event.Momentum()[0].E());
    if(!m_normalize) return scale;
    if(monochromatic) return 1;
    return scale / Normalization(*ctx.beam, pid);
}

// ---------------------------------------------------------------------------
// Covariance groups
// ---------------------------------------------------------------------------

achilles::EigenSystem achilles::SymmetricEigen(Matrix a) {
    const size_t n = a.size();
    Matrix v(n, std::vector<double>(n, 0));
    for(size_t i = 0; i < n; ++i) v[i][i] = 1;

    // Cyclic Jacobi rotations, robust for the small matrices used here
    for(size_t sweep = 0; sweep < 100; ++sweep) {
        double off = 0, diag = 0;
        for(size_t p = 0; p < n; ++p) {
            diag += a[p][p] * a[p][p];
            for(size_t q = p + 1; q < n; ++q) off += a[p][q] * a[p][q];
        }
        if(off <= 1e-30 * std::max(diag, 1e-300)) break;
        for(size_t p = 0; p < n; ++p) {
            for(size_t q = p + 1; q < n; ++q) {
                if(a[p][q] == 0) continue;
                const double theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
                const double t =
                    (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1));
                const double c = 1 / std::sqrt(t * t + 1), s = t * c;
                for(size_t k = 0; k < n; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for(size_t k = 0; k < n; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for(size_t k = 0; k < n; ++k) {
                    const double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }

    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t i, size_t j) { return a[i][i] > a[j][j]; });
    EigenSystem result;
    for(const auto k : order) {
        result.values.push_back(a[k][k]);
        std::vector<double> vec(n);
        size_t imax = 0;
        for(size_t i = 0; i < n; ++i) {
            vec[i] = v[i][k];
            if(std::abs(vec[i]) > std::abs(vec[imax])) imax = i;
        }
        if(vec[imax] < 0)
            for(auto &x : vec) x = -x;
        result.vectors.push_back(vec);
    }
    return result;
}

namespace {

achilles::Matrix ReadMatrix(const YAML::Node &node, size_t n, const std::string &name) {
    auto matrix = node.as<achilles::Matrix>();
    if(matrix.size() != n)
        throw std::runtime_error(fmt::format("Variations: {} must be {}x{}", name, n, n));
    for(const auto &row : matrix)
        if(row.size() != n)
            throw std::runtime_error(fmt::format("Variations: {} must be {}x{}", name, n, n));
    for(size_t i = 0; i < n; ++i)
        for(size_t j = i + 1; j < n; ++j)
            if(std::abs(matrix[i][j] - matrix[j][i]) >
               1e-9 * std::max({std::abs(matrix[i][j]), std::abs(matrix[j][i]), 1e-300}))
                throw std::runtime_error(fmt::format("Variations: {} is not symmetric", name));
    return matrix;
}

} // namespace

achilles::Matrix achilles::ParseCovariance(const YAML::Node &node, size_t n) {
    YAML::Node spec = YAML::Clone(node);
    if(node["File"]) {
        const auto file =
            YAML::LoadFile(Filesystem::FindFile(node["File"].as<std::string>(), "Variations"));
        for(const auto &entry : file) {
            const auto key = entry.first.as<std::string>();
            if(!spec[key]) spec[key] = entry.second;
        }
    }

    const bool has_cov = spec["Covariance"].IsDefined() && !spec["Covariance"].IsNull();
    const bool has_corr = spec["Correlation"].IsDefined() && !spec["Correlation"].IsNull();
    const bool has_unc = spec["Uncertainties"].IsDefined() && !spec["Uncertainties"].IsNull();
    const bool has_hess = spec["Hessian"].IsDefined() && !spec["Hessian"].IsNull();
    if(static_cast<int>(has_cov) + static_cast<int>(has_unc || has_corr) +
           static_cast<int>(has_hess) !=
       1)
        throw std::runtime_error("Variations: Give exactly one of 'Covariance', "
                                 "'Uncertainties' (with optional 'Correlation') or 'Hessian'");

    if(has_cov) return ReadMatrix(spec["Covariance"], n, "Covariance");

    if(has_hess) {
        // Invert through the eigen decomposition, which shares the eigenvectors
        const auto hessian = SymmetricEigen(ReadMatrix(spec["Hessian"], n, "Hessian"));
        Matrix cov(n, std::vector<double>(n, 0));
        for(size_t k = 0; k < n; ++k) {
            if(hessian.values[k] <= 1e-12 * std::abs(hessian.values.front()))
                throw std::runtime_error("Variations: Hessian is not positive definite");
            for(size_t i = 0; i < n; ++i)
                for(size_t j = 0; j < n; ++j)
                    cov[i][j] += hessian.vectors[k][i] * hessian.vectors[k][j] / hessian.values[k];
        }
        return cov;
    }

    if(!has_unc) throw std::runtime_error("Variations: 'Correlation' needs 'Uncertainties'");
    const auto sigma = spec["Uncertainties"].as<std::vector<double>>();
    if(sigma.size() != n)
        throw std::runtime_error(fmt::format("Variations: Expected {} uncertainties", n));
    Matrix corr(n, std::vector<double>(n, 0));
    if(has_corr) {
        corr = ReadMatrix(spec["Correlation"], n, "Correlation");
    } else {
        for(size_t i = 0; i < n; ++i) corr[i][i] = 1;
    }
    Matrix cov(n, std::vector<double>(n, 0));
    for(size_t i = 0; i < n; ++i) {
        if(std::abs(corr[i][i] - 1) > 1e-9)
            throw std::runtime_error("Variations: Correlation diagonal must be 1");
        for(size_t j = 0; j < n; ++j) {
            if(std::abs(corr[i][j]) > 1 + 1e-9)
                throw std::runtime_error("Variations: Correlations must be within [-1, 1]");
            cov[i][j] = corr[i][j] * sigma[i] * sigma[j];
        }
    }
    return cov;
}

achilles::EigenShifts achilles::CovarianceEigenShifts(const Matrix &covariance, double tolerance,
                                                      size_t keep) {
    const auto eigen = SymmetricEigen(covariance);
    if(eigen.values.empty() || eigen.values.front() <= 0)
        throw std::runtime_error("Variations: Covariance has no positive eigenvalues");
    const double largest = eigen.values.front();
    if(eigen.values.back() < -1e-9 * largest)
        throw std::runtime_error(fmt::format(
            "Variations: Covariance is not positive semi-definite (eigenvalue {} vs largest {})",
            eigen.values.back(), largest));

    EigenShifts shifts;
    for(size_t k = 0; k < eigen.values.size(); ++k) {
        if(eigen.values[k] <= tolerance * largest) break;
        if(keep > 0 && shifts.sigmas.size() == keep) break;
        shifts.sigmas.push_back(std::sqrt(eigen.values[k]));
        shifts.directions.push_back(eigen.vectors[k]);
    }
    return shifts;
}

std::unique_ptr<achilles::Variation> achilles::CovarianceTarget::Central() const {
    return std::make_unique<NominalVariation>();
}

achilles::FormFactorTarget::FormFactorTarget(const YAML::Node &node) {
    if(node.IsSequence()) {
        m_parameters = node.as<std::vector<std::string>>();
    } else {
        m_parameters = node["Parameters"].as<std::vector<std::string>>();
        if(node["Central"]) m_central = node["Central"].as<std::vector<double>>();
    }
    if(m_parameters.empty())
        throw std::runtime_error("Variations: FormFactor covariance needs parameters");
    if(!m_central.empty() && m_central.size() != m_parameters.size())
        throw std::runtime_error(
            "Variations: FormFactor covariance needs one central value per parameter");
}

std::unique_ptr<achilles::Variation> achilles::FormFactorTarget::Central() const {
    if(m_central.empty()) return std::make_unique<NominalVariation>();
    return Shifted(std::vector<double>(m_parameters.size(), 0));
}

std::unique_ptr<achilles::Variation>
achilles::FormFactorTarget::Shifted(const std::vector<double> &delta) const {
    FormFactorVariation::Edit edit;
    for(size_t i = 0; i < m_parameters.size(); ++i) {
        if(m_central.empty()) {
            edit.parameters.push_back({m_parameters[i], delta[i], true});
        } else {
            edit.parameters.push_back({m_parameters[i], m_central[i] + delta[i], false});
        }
    }
    return std::make_unique<FormFactorVariation>(edit);
}

achilles::FluxTarget::FluxTarget(const YAML::Node &node) {
    m_edges = node["Bins"].as<std::vector<double>>();
    if(m_edges.size() < 2 || !std::is_sorted(m_edges.begin(), m_edges.end()) ||
       std::adjacent_find(m_edges.begin(), m_edges.end()) != m_edges.end())
        throw std::runtime_error("Variations: Flux 'Bins' must be increasing bin edges");
    if(node["PID"]) m_pid = node["PID"].as<int>();
    if(node["Normalize"]) m_normalize = node["Normalize"].as<bool>();
}

std::vector<std::string> achilles::FluxTarget::Parameters() const {
    std::vector<std::string> names;
    for(size_t i = 0; i + 1 < m_edges.size(); ++i)
        names.push_back(
            fmt::format("Flux[{},{})", FormatValue(m_edges[i]), FormatValue(m_edges[i + 1])));
    return names;
}

std::unique_ptr<achilles::Variation>
achilles::FluxTarget::Shifted(const std::vector<double> &delta) const {
    std::vector<double> scales(delta.size());
    for(size_t i = 0; i < delta.size(); ++i) scales[i] = 1 + delta[i];
    return std::make_unique<FluxVariation>(m_edges, scales, m_pid, m_normalize);
}

achilles::CovarianceGroup::CovarianceGroup(const YAML::Node &node) : VariationGroup(node) {
    const YAML::Node target_node = node["Target"];
    if(!target_node || !target_node.IsMap() || target_node.size() != 1)
        throw std::runtime_error(
            fmt::format("Variations: Group {} needs a 'Target' with exactly one of: {}",
                        GroupName(), fmt::join(CovarianceTargetFactory::List(), ", ")));
    const auto target_name = target_node.begin()->first.as<std::string>();
    std::unique_ptr<CovarianceTarget> target;
    try {
        target = CovarianceTargetFactory::Initialize(target_name, target_node.begin()->second);
    } catch(std::out_of_range &) {
        throw std::runtime_error(
            fmt::format("Variations: Unknown covariance target '{}', did you mean '{}'?",
                        target_name, GetSuggestion(CovarianceTargetFactory::List(), target_name)));
    }

    const double nsigma = node["NSigma"] ? node["NSigma"].as<double>() : 1.0;
    const double tolerance = node["Tolerance"] ? node["Tolerance"].as<double>() : 1e-10;
    const size_t keep = node["Keep"] ? node["Keep"].as<size_t>() : 0;
    if(nsigma <= 0) throw std::runtime_error("Variations: NSigma must be positive");

    const auto covariance = ParseCovariance(node, target->Dimension());
    const auto shifts = CovarianceEigenShifts(covariance, tolerance, keep);
    spdlog::info("Variations: Group {} keeps {} of {} eigen directions", GroupName(),
                 shifts.sigmas.size(), target->Dimension());

    m_parameters = target->Parameters();
    SetCentral(target->Central());
    std::vector<double> flat_vectors;
    for(size_t k = 0; k < shifts.sigmas.size(); ++k) {
        std::vector<double> delta(shifts.directions[k].size());
        for(size_t i = 0; i < delta.size(); ++i)
            delta[i] = nsigma * shifts.sigmas[k] * shifts.directions[k][i];
        std::vector<double> minus(delta.size());
        std::transform(delta.begin(), delta.end(), minus.begin(), [](double x) { return -x; });
        AddMember(fmt::format("ev{}+", k + 1), target->Shifted(delta), nsigma * shifts.sigmas[k]);
        AddMember(fmt::format("ev{}-", k + 1), target->Shifted(minus), -nsigma * shifts.sigmas[k]);
        flat_vectors.insert(flat_vectors.end(), shifts.directions[k].begin(),
                            shifts.directions[k].end());
    }
    m_info_vectors = {{"Sigmas", shifts.sigmas}, {"Eigenvectors", flat_vectors}};

    // The members are NSigma away, so rescale the Hessian band to one standard deviation
    SetDefaultCombination(HessianCombination::Name());
    SetupCombination(node, 1 / nsigma);
}

// ---------------------------------------------------------------------------
// Shorthand notation
// ---------------------------------------------------------------------------

namespace {

bool IsNumber(const YAML::Node &node) {
    double value{};
    return node.IsScalar() && YAML::convert<double>::decode(node, value);
}

/// Group names may not contain whitespace, ':' or '|', so derive them from parameter paths
std::string GroupNameFromPath(const std::string &path) {
    std::string name = path;
    std::replace(name.begin(), name.end(), '/', '.');
    std::replace_if(name.begin(), name.end(), [](char c) { return std::isspace(c) != 0; }, '_');
    return name;
}

/// Round to 15 significant digits so that e.g. 1.05 - 0.1 gives the same double as 0.95
double RoundDecimal(double value) {
    return std::stod(fmt::format("{:.15g}", value));
}

/// Parse "C +- d", "C ± d", "C +u -d" or "C -d +u" into (central, down, up)
std::array<double, 3> ParseUncertainty(const std::string &value) {
    static const std::string number = R"(([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?))";
    static const std::regex symmetric("^\\s*" + number + "\\s*(?:\\+-|\\+/-|±)\\s*" + number +
                                      "\\s*$");
    static const std::regex up_down("^\\s*" + number + "\\s*\\+\\s*" + number + "\\s*-\\s*" +
                                    number + "\\s*$");
    static const std::regex down_up("^\\s*" + number + "\\s*-\\s*" + number + "\\s*\\+\\s*" +
                                    number + "\\s*$");
    std::smatch match;
    if(std::regex_match(value, match, symmetric)) {
        const double central = std::stod(match[1]), delta = std::stod(match[2]);
        return {central, RoundDecimal(central - std::abs(delta)),
                RoundDecimal(central + std::abs(delta))};
    }
    if(std::regex_match(value, match, up_down)) {
        const double central = std::stod(match[1]);
        return {central, RoundDecimal(central - std::stod(match[3])),
                RoundDecimal(central + std::stod(match[2]))};
    }
    if(std::regex_match(value, match, down_up)) {
        const double central = std::stod(match[1]);
        return {central, RoundDecimal(central - std::stod(match[2])),
                RoundDecimal(central + std::stod(match[3]))};
    }
    throw std::runtime_error(
        fmt::format("Variations: Could not parse '{}', expected a number, a list, a range or "
                    "'central +- delta' / 'central +up -down'",
                    value));
}

/// Translate a Sherpa style range {Min, Max, Number | Step} into the explicit Scan
YAML::Node ScanFromRange(const std::string &key, const YAML::Node &range) {
    const auto min = range["Min"].as<double>();
    const auto max = range["Max"].as<double>();
    YAML::Node scan;
    scan["Min"] = min;
    scan["Max"] = max;
    if(range["Number"] && !range["Step"]) {
        scan["Steps"] = range["Number"].as<size_t>();
    } else if(range["Step"] && !range["Number"]) {
        const auto step = range["Step"].as<double>();
        const double intervals = (max - min) / step;
        const auto nintervals = std::llround(intervals);
        if(step <= 0 || nintervals < 1 ||
           std::abs(intervals - static_cast<double>(nintervals)) > 1e-9 * std::max(1.0, intervals))
            throw std::runtime_error(fmt::format(
                "Variations: Step of {} must divide the range of {} into equal intervals", key,
                key));
        scan["Steps"] = nintervals + 1;
    } else {
        throw std::runtime_error(
            fmt::format("Variations: Range of {} needs exactly one of 'Number' or 'Step'", key));
    }
    return scan;
}

std::string SpectralFunctionFile(const std::string &file) {
    return file.find('/') == std::string::npos ? "data/Spectral_Functions/" + file : file;
}

YAML::Node SpectralFunctionPair(const std::string &key, const YAML::Node &value) {
    if(!value.IsSequence() || value.size() != 2)
        throw std::runtime_error(fmt::format(
            "Variations: Spectral function {} expects [proton file, neutron file]", key));
    YAML::Node spec;
    spec["SpectralP"] = SpectralFunctionFile(value[0].as<std::string>());
    spec["SpectralN"] = SpectralFunctionFile(value[1].as<std::string>());
    return spec;
}

[[maybe_unused]] const bool registered_form_factors =
    VariationShorthand::Register("FormFactors", achilles::FormFactorGroup::ExpandShorthand);
[[maybe_unused]] const bool registered_spectral_functions = VariationShorthand::Register(
    "SpectralFunctions", achilles::SpectralFunctionGroup::ExpandShorthand);
[[maybe_unused]] const bool registered_covariance =
    VariationShorthand::Register("Covariance", achilles::CovarianceGroup::ExpandShorthand);

} // namespace

std::map<std::string, VariationShorthand::Expander> &VariationShorthand::Registry() {
    static std::map<std::string, Expander> registry;
    return registry;
}

bool VariationShorthand::Register(const std::string &key, Expander expander) {
    return Registry().emplace(key, std::move(expander)).second;
}

std::vector<std::string> VariationShorthand::Keys() {
    std::vector<std::string> keys;
    for(const auto &entry : Registry()) keys.push_back(entry.first);
    return keys;
}

YAML::Node VariationShorthand::Expand(const YAML::Node &shorthand) {
    YAML::Node groups(YAML::NodeType::Sequence);
    for(const auto &entry : shorthand) {
        const auto key = entry.first.as<std::string>();
        auto it = Registry().find(key);
        if(it == Registry().end()) {
            throw std::runtime_error(fmt::format(
                "Variations: Unknown variation kind '{}', did you mean '{}'? Known kinds: {}", key,
                GetSuggestion(Keys(), key), fmt::join(Keys(), ", ")));
        }
        if(!entry.second.IsMap())
            throw std::runtime_error(fmt::format("Variations: '{}' expects a map", key));
        it->second(entry.second, groups);
    }
    return groups;
}

void FormFactorGroup::ExpandShorthand(const YAML::Node &entries, YAML::Node &groups) {
    static const std::vector<std::string> slots{
        FFTypeToString(FFType::vector),         FFTypeToString(FFType::axial),
        FFTypeToString(FFType::coherent),       FFTypeToString(FFType::resonancevector),
        FFTypeToString(FFType::resonanceaxial), FFTypeToString(FFType::mecvector),
        FFTypeToString(FFType::mecaxial),       FFTypeToString(FFType::hyperon)};

    for(const auto &entry : entries) {
        const auto key = entry.first.as<std::string>();
        const YAML::Node value = entry.second;
        YAML::Node group;
        group["Type"] = Name();

        if(key.find('/') == std::string::npos) {
            // Functional form alternatives for one of the form factor slots
            if(std::find(slots.begin(), slots.end(), key) == slots.end())
                throw std::runtime_error(fmt::format(
                    "Variations: '{}' is neither a 'Block/parameter' path nor one of the form "
                    "factor types: {}",
                    key, fmt::join(slots, ", ")));
            group["Name"] = key;
            YAML::Node names(YAML::NodeType::Sequence);
            if(value.IsSequence()) {
                names = value;
            } else {
                names.push_back(value);
            }
            for(const auto &name : names) {
                YAML::Node alternative;
                alternative["Name"] = name.as<std::string>();
                alternative["Overrides"][key] = name.as<std::string>();
                group["Alternatives"].push_back(alternative);
            }
        } else {
            group["Name"] = GroupNameFromPath(key);
            group["Parameter"] = key;
            if(value.IsSequence()) {
                group["Values"] = value;
            } else if(value.IsMap()) {
                for(const auto &option : value) {
                    const auto name = option.first.as<std::string>();
                    if(name == "Min" || name == "Max" || name == "Number" || name == "Step")
                        continue;
                    if(name != "Values" && name != "MinMax" && name != "Central" &&
                       name != "Combination")
                        throw std::runtime_error(fmt::format(
                            "Variations: Unknown option '{}' for parameter {}", name, key));
                    group[name] = option.second;
                }
                if(value["Min"] || value["Max"]) group["Scan"] = ScanFromRange(key, value);
            } else if(IsNumber(value)) {
                group["Values"].push_back(value);
            } else {
                const auto [central, down, up] = ParseUncertainty(value.as<std::string>());
                group["Central"] = central;
                group["MinMax"].push_back(down);
                group["MinMax"].push_back(up);
            }
        }
        groups.push_back(group);
    }
}

void achilles::CovarianceGroup::ExpandShorthand(const YAML::Node &entries, YAML::Node &groups) {
    const auto targets = CovarianceTargetFactory::List();
    for(const auto &entry : entries) {
        const auto name = entry.first.as<std::string>();
        if(!entry.second.IsMap())
            throw std::runtime_error(fmt::format("Variations: Covariance {} expects a map", name));
        YAML::Node group;
        group["Name"] = name;
        group["Type"] = Name();
        for(const auto &option : entry.second) {
            const auto key = option.first.as<std::string>();
            if(std::find(targets.begin(), targets.end(), key) != targets.end()) {
                if(group["Target"])
                    throw std::runtime_error(
                        fmt::format("Variations: Covariance {} has more than one target", name));
                group["Target"][key] = option.second;
            } else {
                group[key] = option.second;
            }
        }
        groups.push_back(group);
    }
}

void SpectralFunctionGroup::ExpandShorthand(const YAML::Node &entries, YAML::Node &groups) {
    YAML::Node group;
    group["Name"] = "SF";
    group["Type"] = Name();
    for(const auto &entry : entries) {
        const auto key = entry.first.as<std::string>();
        if(key == "Combination") {
            group["Combination"] = entry.second;
        } else if(key == "Central") {
            group["Central"] = SpectralFunctionPair(key, entry.second);
        } else {
            auto alternative = SpectralFunctionPair(key, entry.second);
            alternative["Name"] = key;
            group["Alternatives"].push_back(alternative);
        }
    }
    groups.push_back(group);
}

// ---------------------------------------------------------------------------
// Handler
// ---------------------------------------------------------------------------

VariationHandler::VariationHandler(const YAML::Node &variations) {
    const YAML::Node groups =
        variations.IsMap() ? VariationShorthand::Expand(variations) : variations;
    if(variations.IsMap())
        spdlog::debug("Variations: Expanded shorthand to\n{}", YAML::Dump(groups));
    for(const auto &node : groups) {
        const auto type = node["Type"].as<std::string>();
        std::unique_ptr<VariationGroup> group;
        try {
            group = VariationGroupFactory::Initialize(type, node);
        } catch(std::out_of_range &) {
            spdlog::error("Variations: Requested variation type \"{}\", did you mean \"{}\"", type,
                          GetSuggestion(VariationGroupFactory::List(), type));
            throw;
        }
        AddGroup(std::move(group));
    }
}

void VariationHandler::AddGroup(std::unique_ptr<VariationGroup> group) {
    for(const auto &other : m_groups)
        if(other->GroupName() == group->GroupName())
            throw std::runtime_error(
                fmt::format("Variations: Duplicate group name {}", group->GroupName()));
    spdlog::info("Variations: Added {} group \"{}\" with {} weights ({} band)", group->Type(),
                 group->GroupName(), group->NWeights(), group->GetCombination().Type());
    m_nweights += group->NWeights();
    m_groups.push_back(std::move(group));
}

std::vector<std::string> VariationHandler::WeightNames() const {
    std::vector<std::string> names;
    for(const auto &group : m_groups) {
        auto group_names = group->WeightNames();
        names.insert(names.end(), group_names.begin(), group_names.end());
    }
    return names;
}

std::vector<achilles::VariationGroupInfo> VariationHandler::Info() const {
    std::vector<VariationGroupInfo> info;
    for(const auto &group : m_groups) info.push_back(group->Info());
    return info;
}

void VariationHandler::Initialize(const NuclearModel &model) const {
    for(const auto &group : m_groups) group->Initialize(model);
}

std::vector<double> VariationHandler::Evaluate(const VariationContext &ctx) const {
    std::vector<double> ratios;
    ratios.reserve(m_nweights);
    for(const auto &group : m_groups) group->Evaluate(ctx, ratios);
    return ratios;
}
