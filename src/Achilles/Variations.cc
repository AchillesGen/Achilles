// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Variations.hh"
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
    if(edit.parameter.empty()) return merged;

    // Setting a parameter to the value used for generation leaves the model unchanged
    if(only_parameter && Lookup(merged, edit.parameter).as<double>() == edit.value) return {};
    return SetParameter(merged, edit.parameter, edit.value);
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

void VariationGroup::SetupCombination(const YAML::Node &node) {
    m_combination = Combination::Build(node["Combination"], m_default_combination);
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
    if(!m_parameter.empty()) {
        info.parameter = m_parameter;
        info.member_values = m_values;
    }
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
        edit.parameter = node["Parameter"].as<std::string>();
        edit.value = node["Value"].as<double>();
    }
    if(edit.overrides.IsNull() && edit.replacement.IsNull() && edit.parameter.empty())
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
            edit.parameter = m_parameter;
            edit.value = value;
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
