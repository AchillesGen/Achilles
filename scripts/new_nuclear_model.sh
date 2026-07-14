#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
#
# new_nuclear_model.sh — Scaffold a new Achilles nuclear-model plugin.
#
# Prompts for a model name and writes a compilable C++ skeleton that:
#   * derives from achilles::NuclearModel and self-registers with the nuclear
#     model factory (via RegistrableNuclearModel<>),
#   * implements every pure-virtual method of NuclearModel with a TODO stub,
#   * exports the Achilles version symbol the PluginManager checks on load.
#
# Compile + install the result with the companion script:
#   scripts/install_nuclear_model.sh <ModelName>NuclearModel.cc
#
# Usage:
#   scripts/new_nuclear_model.sh [-n NAME] [-o OUTPUT_DIR]
#
#   -n NAME   Model name (a valid C++ identifier). Prompted for if omitted.
#   -o DIR    Directory to write the source file into (default: current directory).
#
set -euo pipefail

NAME=""
OUTPUT_DIR="$(pwd)"

while getopts ":n:o:h" opt; do
    case "$opt" in
        n) NAME="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :)  echo "Option -$OPTARG requires an argument" >&2; exit 1 ;;
    esac
done

log() { printf '\033[1;34m[new_nuclear_model]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[new_nuclear_model] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

# Ask for the model name if it was not supplied on the command line.
if [[ -z "$NAME" ]]; then
    printf 'Name for the new nuclear model (C++ identifier, e.g. MySpectral): ' >&2
    read -r NAME
fi

# A model name becomes a C++ class name and part of a library file name, so it
# must be a valid identifier.
[[ "$NAME" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] \
    || die "'$NAME' is not a valid C++ identifier (letters, digits, underscore; not starting with a digit)."

mkdir -p "$OUTPUT_DIR"
DEST="$OUTPUT_DIR/${NAME}NuclearModel.cc"
[[ -e "$DEST" ]] && die "Refusing to overwrite existing file: $DEST"

cat > "$DEST" <<EOF
// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// ${NAME} — an Achilles nuclear-model plugin scaffolded by new_nuclear_model.sh.
//
// Build and install it with:
//     scripts/install_nuclear_model.sh ${NAME}NuclearModel.cc
//
// Fill in the TODOs below, then select the model from a run card with:
//     NuclearModels:
//       - Model: ${NAME}
#include "Achilles/NuclearModel.hh"
#include "Achilles/Version.hh"

#include "spdlog/spdlog.h"

// Exports the ExpectedVersion symbol the PluginManager validates before loading.
EXPORT_ACHILLES_VERSION()

namespace achilles {

// Inheriting from RegistrableNuclearModel<${NAME}> self-registers the model with
// the nuclear-model factory at static-initialisation time; no explicit
// registration call is needed.
class ${NAME} : public NuclearModel, RegistrableNuclearModel<${NAME}> {
  public:
    // The factory hands the raw run-card node; the second argument is the parsed
    // form-factor configuration produced by Construct() below.
    ${NAME}(const YAML::Node &config, const YAML::Node &form_factor,
            FormFactorBuilder &builder = FormFactorBuilder::Instance())
        : NuclearModel(form_factor, builder) {
        // TODO: read any model-specific parameters, e.g.
        //   auto value = config["NuclearModel"]["MyParameter"].as<double>();
        (void)config;
    }

    // --- Required NuclearModel interface ------------------------------------ //

    // Which interaction mode this model provides (Quasielastic, Coherent, ...).
    NuclearMode Mode() const override { return NuclearMode::Quasielastic; } // TODO

    // Name of a registered phase space compatible with this mode (see
    // 'achilles --display-ps'). Throw achilles::InvalidChannel for unsupported
    // nuclei.
    std::string PhaseSpace(PID) const override { return PSName(); } // TODO

    // Heart of the model: build the hadronic current(s) for the transferred
    // four-momentum qVec, keyed by the form-factor group in 'ff'.
    Currents CalcCurrents(const std::vector<Particle> &had_in,
                          const std::vector<Particle> &had_out,
                          const std::vector<Particle> &spectators, const FourVector &qVec,
                          const FFInfoMap &ff) const override {
        (void)had_in;
        (void)had_out;
        (void)spectators;
        (void)qVec;
        (void)ff;
        // TODO: for each entry in 'ff', evaluate the form factors with
        //   EvalFormFactor(-qVec.M2()) and CouplingsFF(...), then fill a Current.
        throw std::runtime_error("${NAME}: CalcCurrents is not implemented yet");
    }

    // Number of spin combinations summed over (1 for spin-0 coherent, 4 for a
    // spin-1/2 nucleon in and out, ...).
    size_t NSpins() const override { return 4; } // TODO

    // Phase-space weight of the sampled initial nuclear state.
    double InitialStateWeight(const std::vector<Particle> &had_in,
                              const std::vector<Particle> &had_out, size_t nprotons,
                              size_t nneutrons) const override {
        (void)had_in;
        (void)had_out;
        (void)nprotons;
        (void)nneutrons;
        return 1.0; // TODO
    }

    // Metadata.
    std::string GetName() const override { return ${NAME}::Name(); }
    std::string PSName() const override { return "OneBodySpectral"; } // TODO
    std::string InspireHEP() const override { return ""; } // TODO: BibTeX/INSPIRE key

    // --- Required factory hooks --------------------------------------------- //

    // Called by the factory to build an instance from a run-card node. The base
    // class helper LoadFormFactor() resolves the form-factor file referenced by
    // config["NuclearModel"]["FormFactorFile"].
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &config) {
        auto form_factor = LoadFormFactor(config);
        return std::make_unique<${NAME}>(config, form_factor);
    }

    // The name used to select this model from a run card.
    static std::string Name() { return "${NAME}"; }
};

} // namespace achilles
EOF

log "Created nuclear-model template: $DEST"
log "Next: scripts/install_nuclear_model.sh \"$DEST\""
echo "$DEST"
