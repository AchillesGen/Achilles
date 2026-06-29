# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""Physics validation for Achilles run card configurations.

Two-layer validation:
  1. Rules.yml — structural rules shared with the C++ engine
  2. Physics rules — domain-specific checks (MCP-only)
"""

import importlib.resources
from pathlib import Path

import yaml

from .defaults import (
    NUCLEUS_SPECTRAL_FUNCTIONS,
    SINGLE_NUCLEON_TARGETS,
    VALID_CUT_TYPES,
)
from .run_card import IncludeTag


class ValidationError(Exception):
    """Raised when a run card configuration is invalid."""

    pass


def validate_config(config: dict) -> None:
    """Validate a complete run card config dict.

    Applies Rules.yml constraints, then physics-specific checks.
    Raises ValidationError with a descriptive message.
    """
    _apply_rules_yml(config)
    _validate_beam(config)
    _validate_nucleus_model_compatibility(config)
    _validate_cascade(config)
    _validate_hard_cuts(config)
    _validate_parameter_ranges(config)


# ---------------------------------------------------------------------------
# Layer 1: Rules.yml validation (ported from C++ SettingsValidator)
# ---------------------------------------------------------------------------


def _load_rules_yml() -> list[dict]:
    """Load Rules.yml from the Achilles data directory."""
    try:
        data_dir = Path(importlib.resources.files("achilles")) / "data"
    except (TypeError, AttributeError):
        data_dir = Path(__file__).resolve().parent.parent.parent.parent / "data"
    rules_path = data_dir / "Rules.yml"
    # For editable installs, data/ may be at the repo root
    if not rules_path.exists():
        repo_root = Path(__file__).resolve().parent.parent.parent.parent
        rules_path = repo_root / "data" / "Rules.yml"
    if not rules_path.exists():
        return []
    with open(rules_path) as f:
        doc = yaml.safe_load(f)
    return doc.get("Rules", [])


def _resolve_path(config, path: str):
    """Navigate a config dict using slash-separated paths like 'Main/NEvents'.

    Returns the value at the path, or raises KeyError if not found.
    Supports wildcard '*' matching for sequences/maps.
    """
    segments = [s for s in path.split("/") if s]
    return _resolve_segments(config, segments)


def _resolve_segments(node, segments: list[str]):
    """Recursively resolve path segments against a node."""
    if not segments:
        return node

    segment = segments[0]
    rest = segments[1:]

    if segment == "*":
        # Wildcard: iterate over list items or dict values
        results = []
        if isinstance(node, list):
            for item in node:
                try:
                    results.append(_resolve_segments(item, rest))
                except (KeyError, TypeError, IndexError):
                    pass
        elif isinstance(node, dict):
            for val in node.values():
                try:
                    results.append(_resolve_segments(val, rest))
                except (KeyError, TypeError, IndexError):
                    pass
        if not results:
            raise KeyError(f"Wildcard '*' matched nothing")
        return results

    # Skip IncludeTag nodes — they can't be traversed
    if isinstance(node, IncludeTag):
        raise KeyError(f"Cannot traverse IncludeTag for '{segment}'")

    if isinstance(node, dict):
        if segment in node:
            return _resolve_segments(node[segment], rest)
        raise KeyError(f"Key '{segment}' not found")

    if isinstance(node, list):
        # Try numeric index
        if segment.isdigit():
            idx = int(segment)
            if idx < len(node):
                return _resolve_segments(node[idx], rest)
        raise KeyError(f"Cannot index list with '{segment}'")

    raise KeyError(f"Cannot traverse {type(node).__name__} with '{segment}'")


def _path_exists(config, path: str) -> bool:
    """Check if a path exists in the config dict."""
    try:
        val = _resolve_path(config, path)
        return val is not None
    except (KeyError, TypeError, IndexError):
        return False


def _apply_rules_yml(config: dict) -> None:
    """Load data/Rules.yml and apply all rules to the config."""
    rules = _load_rules_yml()
    for rule in rules:
        rule_type = rule.get("Type", "")
        if rule_type == "RequiredFields":
            _apply_required_fields(config, rule)
        elif rule_type == "RangeConstraint":
            _apply_range_constraint(config, rule)
        elif rule_type == "ConditionalInteractionOption":
            # Cascade interactions use IncludeTag in the config dict,
            # so we can't inspect them here. Skip — the C++ engine
            # validates these after resolving includes.
            pass


def _apply_required_fields(config: dict, rule: dict) -> None:
    """Check that all required fields exist in the config."""
    fields = rule.get("Fields", [])
    for field in fields:
        if not _path_exists(config, field):
            raise ValidationError(
                f"Required field '{field}' is missing."
            )


def _apply_range_constraint(config: dict, rule: dict) -> None:
    """Apply a scalar range constraint from Rules.yml."""
    if "Path" not in rule:
        # List-based range constraints reference Cascade/Interactions
        # which is an IncludeTag — skip for MCP validation.
        return

    path = rule["Path"]
    optional = rule.get("Optional", False)

    try:
        value = _resolve_path(config, path)
    except (KeyError, TypeError, IndexError):
        if optional:
            return
        raise ValidationError(
            f"Required field '{path}' is missing."
        )

    if isinstance(value, IncludeTag):
        return  # Can't validate include references

    try:
        value = float(value)
    except (ValueError, TypeError):
        return  # Non-numeric values can't be range-checked

    min_val = rule.get("Min")
    max_val = rule.get("Max")
    error_msg = rule.get("ErrorMessage", f"Value at '{path}' is out of range.")

    if min_val is not None and value < float(min_val):
        raise ValidationError(error_msg)
    if max_val is not None and value > float(max_val):
        raise ValidationError(error_msg)


# ---------------------------------------------------------------------------
# Layer 2: Physics-specific validation
# ---------------------------------------------------------------------------


def _get_nucleus_name(config: dict) -> str | None:
    """Extract the nucleus name from the config.

    Handles both IncludeTag references and inline nucleus configs.
    """
    nuclei = config.get("Nuclei", [])
    if not nuclei:
        return None

    nucleus = nuclei[0]
    nuc_val = nucleus.get("Nucleus") if isinstance(nucleus, dict) else None
    if nuc_val is None:
        return None

    if isinstance(nuc_val, IncludeTag):
        # e.g. IncludeTag("data/default/12C.yml") -> "12C"
        stem = Path(nuc_val.path).stem
        return stem

    if isinstance(nuc_val, dict):
        return nuc_val.get("Name")

    return None


def _get_nuclear_models(config: dict) -> list[str]:
    """Get list of nuclear model types from the config."""
    models = config.get("NuclearModels", [])
    result = []
    for entry in models:
        if isinstance(entry, dict):
            nm = entry.get("NuclearModel", {})
            if isinstance(nm, dict):
                model = nm.get("Model")
                if model:
                    result.append(model)
    return result


def _validate_beam(config: dict) -> None:
    """Validate beam energy settings."""
    beams = config.get("Beams", [])
    if not beams:
        return

    for beam_entry in beams:
        beam = beam_entry.get("Beam", {}) if isinstance(beam_entry, dict) else {}
        params = beam.get("Beam Params", {})
        if not isinstance(params, dict):
            continue

        beam_type = params.get("Type", "")

        if beam_type == "Monochromatic":
            energy = params.get("Energy")
            if energy is not None and energy <= 0:
                raise ValidationError("Beam energy must be positive.")

        elif beam_type == "FlatFlux":
            min_e = params.get("MinEnergy")
            max_e = params.get("MaxEnergy")
            if min_e is not None and min_e <= 0:
                raise ValidationError(
                    "FlatFlux minimum energy must be positive."
                )
            if min_e is not None and max_e is not None and max_e <= min_e:
                raise ValidationError(
                    "FlatFlux maximum energy must be greater than minimum energy."
                )

        elif beam_type == "Spectrum":
            histogram = params.get("Histogram")
            hepdata = params.get("HepData")
            if histogram:
                ext = Path(histogram).suffix
                if ext not in (".dat", ".csv"):
                    raise ValidationError(
                        f"Spectrum histogram file must be a .dat or .csv file, "
                        f"got '{ext}'."
                    )
            elif hepdata:
                ext = Path(hepdata).suffix
                if ext not in (".yaml", ".yml"):
                    raise ValidationError(
                        f"Spectrum HepData file must be a .yaml file, "
                        f"got '{ext}'."
                    )


def _validate_nucleus_model_compatibility(config: dict) -> None:
    """Check that QESpectral is only used with nuclei that have spectral data."""
    nucleus = _get_nucleus_name(config)
    if nucleus is None:
        return

    models = _get_nuclear_models(config)

    for model in models:
        if model == "QESpectral" and nucleus in NUCLEUS_SPECTRAL_FUNCTIONS:
            if NUCLEUS_SPECTRAL_FUNCTIONS[nucleus] is None:
                raise ValidationError(
                    f"QESpectral requires spectral function data not available "
                    f"for {nucleus}. Use FortranModel instead."
                )


def _validate_cascade(config: dict) -> None:
    """Block cascade on single-nucleon targets."""
    cascade = config.get("Cascade", {})
    if not isinstance(cascade, dict):
        return

    run_cascade = cascade.get("Run", False)
    if not run_cascade:
        return

    nucleus = _get_nucleus_name(config)
    if nucleus in SINGLE_NUCLEON_TARGETS:
        raise ValidationError(
            f"Cannot run cascade on single-nucleon target {nucleus}."
        )


def _validate_hard_cuts(config: dict) -> None:
    """Validate hard cut definitions."""
    cuts = config.get("HardCuts", [])
    if not cuts:
        return

    for i, cut in enumerate(cuts):
        if not isinstance(cut, dict):
            raise ValidationError(f"Hard cut {i} must be a dictionary.")

        cut_type = cut.get("Type")
        if cut_type is None:
            raise ValidationError(
                f"Hard cut {i} is missing required field 'Type'."
            )
        if cut_type not in VALID_CUT_TYPES:
            raise ValidationError(
                f"Hard cut {i} has invalid Type '{cut_type}'. "
                f"Must be one of: {', '.join(sorted(VALID_CUT_TYPES))}."
            )

        if "PIDs" not in cut:
            raise ValidationError(
                f"Hard cut {i} is missing required field 'PIDs'."
            )

        cut_range = cut.get("range")
        if cut_range is None:
            raise ValidationError(
                f"Hard cut {i} is missing required field 'range'."
            )
        if not isinstance(cut_range, list) or len(cut_range) != 2:
            raise ValidationError(
                f"Hard cut {i} 'range' must be a 2-element list [min, max]."
            )
        if cut_range[0] >= cut_range[1]:
            raise ValidationError(
                f"Hard cut {i} range[0] must be less than range[1]."
            )

        # Type-specific range validation
        if cut_type == "AngleTheta":
            if cut_range[0] < 0 or cut_range[1] > 180:
                raise ValidationError(
                    "AngleTheta range must be within [0, 180] degrees."
                )
        elif cut_type in ("Energy", "Momentum", "TransverseMomentum"):
            if cut_range[0] < 0:
                raise ValidationError(
                    f"{cut_type} range values must be non-negative."
                )


def _validate_parameter_ranges(config: dict) -> None:
    """Validate parameter ranges under Options/."""
    options = config.get("Options", {})
    if not isinstance(options, dict):
        return

    # Options/Initialize/Accuracy
    init = options.get("Initialize", {})
    if isinstance(init, dict):
        accuracy = init.get("Accuracy")
        if accuracy is not None:
            if accuracy <= 0 or accuracy >= 1:
                raise ValidationError(
                    "Integration accuracy must be between 0 and 1."
                )

    # Options/Unweighting/percentile
    unweight = options.get("Unweighting", {})
    if isinstance(unweight, dict):
        percentile = unweight.get("percentile")
        if percentile is not None:
            if percentile <= 0 or percentile > 100:
                raise ValidationError(
                    "Unweighting percentile must be between 0 and 100."
                )
