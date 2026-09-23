# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""Achilles MCP server — exposes neutrino event generation via MCP tools and resources."""

import importlib.resources
import json
import os
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from .defaults import (
    AVAILABLE_NUCLEI,
    BEAM_TYPES,
    CASCADE_PROBABILITY_TYPES,
    FORM_FACTOR_AXIAL,
    FORM_FACTOR_COHERENT,
    FORM_FACTOR_VECTOR,
    FORTRAN_MODEL_NAMES,
    NUCLEAR_MODELS,
    OUTPUT_FORMATS,
)
from .run_card import RunCardBuilder
from .runner import format_result, run_generation
from .validation import ValidationError, validate_config

mcp = FastMCP("achilles-neutrino-generator")

# Server-side run card builder state
_builder = RunCardBuilder()

# Resolve data directories from the installed achilles package
_ACHILLES_ROOT = Path(importlib.resources.files("achilles"))
_DATA_DIR = _ACHILLES_ROOT / "data"
_EXAMPLES_DIR = _ACHILLES_ROOT.parent / "examples"  # examples at install root
_FLUX_DIR = _ACHILLES_ROOT / "flux"


# ---------------------------------------------------------------------------
# High-level tool: one-shot event generation
# ---------------------------------------------------------------------------


@mcp.tool()
async def generate_events(
    beam_type: str,
    nucleus: str = "12C",
    nevents: int = 10000,
    beam_energy: float | None = None,
    beam_params_type: str = "Monochromatic",
    histogram: str | None = None,
    hepdata: str | None = None,
    min_energy: float | None = None,
    max_energy: float | None = None,
    nuclear_model: str = "QESpectral",
    output_format: str = "NuHepMC",
    run_cascade: bool = False,
    hard_cuts: list[dict] | None = None,
    seed: int = 12345678,
    accuracy: float = 0.01,
    output_dir: str | None = None,
) -> str:
    """Generate neutrino/electron-nucleus scattering events with sensible defaults.

    Args:
        beam_type: Probe particle — "electron", "nue", "nuebar", "numu", or "numubar".
        nucleus: Target nucleus — "12C", "16O", "40Ar", "1H", or "1N".
        nevents: Number of events to generate.
        beam_energy: Beam energy in MeV (required for Monochromatic beams).
        beam_params_type: Beam profile — "Monochromatic", "Spectrum", or "FlatFlux".
        histogram: Path to a .dat flux histogram file (for Spectrum beams).
        hepdata: Path to a HepData .yaml flux file (for Spectrum beams).
        min_energy: Minimum energy in MeV (required for FlatFlux beams).
        max_energy: Maximum energy in MeV (required for FlatFlux beams).
        nuclear_model: Nuclear interaction model — "QESpectral" or "FortranModel".
        output_format: Output file format — "NuHepMC" or "HepMC3".
        run_cascade: Whether to run intranuclear cascade.
        hard_cuts: Optional kinematic cuts, e.g. [{"Type":"AngleTheta","PIDs":11,"range":[10,90]}].
        seed: Random number seed.
        accuracy: Integration accuracy target.
        output_dir: Directory for output files (defaults to a temp directory).

    Returns:
        A summary of the generation run including output file location and metadata.
    """
    builder = RunCardBuilder()
    builder.set_beam(
        beam_type, energy=beam_energy, beam_params_type=beam_params_type,
        histogram=histogram, hepdata=hepdata,
        min_energy=min_energy, max_energy=max_energy,
    )
    builder.set_nucleus(nucleus)

    # Build a descriptive output filename
    if beam_params_type == "Monochromatic":
        tag = f"{int(beam_energy)}"
    elif histogram:
        tag = Path(histogram).stem
    elif hepdata:
        tag = Path(hepdata).stem
    else:
        tag = f"{int(min_energy)}_{int(max_energy)}"
    builder.set_main(
        nevents=nevents,
        output_format=output_format,
        output_name=f"achilles_{beam_type}_{nucleus}_{tag}.hepmc",
    )
    builder.set_nuclear_model(nuclear_model)
    builder.set_cascade(run=run_cascade)
    builder.set_initialize(seed=seed, accuracy=accuracy)
    if hard_cuts:
        builder.set_hard_cuts(hard_cuts)

    try:
        validate_config(builder.get_config())
    except ValidationError as e:
        return f"Validation error: {e}"

    run_card_yaml = builder.to_yaml()
    result = await run_generation(run_card_yaml, output_dir=output_dir)
    return format_result(result)


# ---------------------------------------------------------------------------
# Low-level tools: run card builder
# ---------------------------------------------------------------------------


@mcp.tool()
async def runcard_reset() -> str:
    """Reset the run card builder to default values. Returns the full default config as YAML."""
    _builder.reset()
    return _builder.to_yaml()


@mcp.tool()
async def runcard_get() -> str:
    """Get the current run card configuration as YAML."""
    return _builder.to_yaml()


@mcp.tool()
async def runcard_set_main(
    nevents: int | None = None,
    hard_cuts_enabled: bool | None = None,
    event_cuts_enabled: bool | None = None,
    do_rotate: bool | None = None,
    run_decays: bool | None = None,
    output_format: str | None = None,
    output_name: str | None = None,
    output_zipped: bool | None = None,
) -> str:
    """Configure the Main section of the run card. Only provided fields are updated.

    Args:
        nevents: Number of events to generate.
        hard_cuts_enabled: Enable hard kinematic cuts.
        event_cuts_enabled: Enable event-level cuts.
        do_rotate: Rotate events to lab frame.
        run_decays: Run particle decays.
        output_format: "NuHepMC" or "HepMC3".
        output_name: Output filename.
        output_zipped: Compress output file.
    """
    _builder.set_main(
        nevents=nevents,
        hard_cuts_enabled=hard_cuts_enabled,
        event_cuts_enabled=event_cuts_enabled,
        do_rotate=do_rotate,
        run_decays=run_decays,
        output_format=output_format,
        output_name=output_name,
        output_zipped=output_zipped,
    )
    return _builder.to_yaml()


@mcp.tool()
async def runcard_set_beam(
    beam_type: str,
    beam_params_type: str = "Monochromatic",
    energy: float | None = None,
    histogram: str | None = None,
    hepdata: str | None = None,
    min_energy: float | None = None,
    max_energy: float | None = None,
) -> str:
    """Set the beam configuration. Also updates Processes/Leptons to match the beam type.

    Args:
        beam_type: "electron", "nue", "nuebar", "numu", or "numubar".
        beam_params_type: "Monochromatic", "Spectrum", or "FlatFlux".
        energy: Beam energy in MeV (required for Monochromatic).
        histogram: Path to a .dat flux histogram file (for Spectrum beams).
            Available files in the flux/ directory include miniboone.dat,
            minerva_numu_fhc.dat, T2K_nu.dat, etc.
        hepdata: Path to a HepData .yaml flux file (for Spectrum beams).
            Example: flux/microboone_flux_numu.yaml.
        min_energy: Minimum energy in MeV (required for FlatFlux).
        max_energy: Maximum energy in MeV (required for FlatFlux).
    """
    _builder.set_beam(
        beam_type, energy=energy, beam_params_type=beam_params_type,
        histogram=histogram, hepdata=hepdata,
        min_energy=min_energy, max_energy=max_energy,
    )
    return _builder.to_yaml()


@mcp.tool()
async def runcard_set_nucleus(
    name: str,
    use_defaults: bool = True,
) -> str:
    """Set the target nucleus.

    Args:
        name: Nucleus name — "12C", "16O", "40Ar", "1H", or "1N".
        use_defaults: If True, use the shipped default configuration for this nucleus.
    """
    _builder.set_nucleus(name, use_defaults)
    return _builder.to_yaml()


@mcp.tool()
async def runcard_set_nuclear_model(
    model: str,
    fortran_model_name: str | None = None,
    form_factor_file: str = "FormFactors.yml",
    spectral_p: str | None = None,
    spectral_n: str | None = None,
    config_file: str | None = None,
    ward: str = "None",
) -> str:
    """Set the nuclear interaction model, replacing any existing models.

    Args:
        model: "QESpectral" or "FortranModel".
        fortran_model_name: Required if model is "FortranModel" — e.g. "RES_Spectral_Func".
        form_factor_file: Path to form factor YAML file.
        spectral_p: Path to proton spectral function data file.
        spectral_n: Path to neutron spectral function data file.
        config_file: Path to nuclear configuration data file.
        ward: Ward identity enforcement — "None" or other options.
    """
    _builder.set_nuclear_model(
        model, fortran_model_name, form_factor_file,
        spectral_p, spectral_n, config_file, ward,
    )
    return _builder.to_yaml()


@mcp.tool()
async def runcard_add_nuclear_model(
    model: str,
    fortran_model_name: str | None = None,
    form_factor_file: str = "FormFactors.yml",
    spectral_p: str | None = None,
    spectral_n: str | None = None,
    config_file: str | None = None,
    ward: str = "None",
) -> str:
    """Append an additional nuclear model to the run card.

    Args:
        model: "QESpectral" or "FortranModel".
        fortran_model_name: Required if model is "FortranModel".
        form_factor_file: Path to form factor YAML file.
        spectral_p: Path to proton spectral function data file.
        spectral_n: Path to neutron spectral function data file.
        config_file: Path to nuclear configuration data file.
        ward: Ward identity enforcement.
    """
    _builder.add_nuclear_model(
        model, fortran_model_name, form_factor_file,
        spectral_p, spectral_n, config_file, ward,
    )
    return _builder.to_yaml()


@mcp.tool()
async def runcard_set_cascade(
    run: bool,
    step: float = 0.04,
    probability: str = "Cylinder",
    in_medium: str = "None",
    potential_prop: bool = False,
    algorithm: str = "Base",
) -> str:
    """Configure the intranuclear cascade (final-state interactions).

    Args:
        run: Whether to run the cascade.
        step: Cascade step size in fm.
        probability: Interaction probability model — "Cylinder" or "Gaussian".
        in_medium: In-medium corrections — "None" or other options.
        potential_prop: Use potential propagation.
        algorithm: Cascade algorithm — "Base".
    """
    _builder.set_cascade(
        run, step, probability, in_medium, potential_prop, algorithm,
    )
    return _builder.to_yaml()


@mcp.tool()
async def runcard_set_hard_cuts(cuts: list[dict]) -> str:
    """Set hard kinematic cuts. Pass an empty list to clear all cuts.

    Args:
        cuts: List of cut definitions. Each dict should have "Type", "PIDs", and "range" keys.
              Example: [{"Type": "AngleTheta", "PIDs": 11, "range": [10, 90]}]
    """
    _builder.set_hard_cuts(cuts)
    return _builder.to_yaml()


@mcp.tool()
async def runcard_set_initialize(
    seed: int | None = None,
    accuracy: float | None = None,
) -> str:
    """Set integration initialization parameters.

    Args:
        seed: Random number seed for reproducibility.
        accuracy: Target integration accuracy (e.g. 0.01 for 1%).
    """
    _builder.set_initialize(seed, accuracy)
    return _builder.to_yaml()


@mcp.tool()
async def runcard_set_unweighting(
    name: str = "Percentile",
    percentile: float = 99,
) -> str:
    """Configure the event unweighting strategy.

    Args:
        name: Unweighting algorithm name.
        percentile: Percentile threshold for unweighting.
    """
    _builder.set_unweighting(name, percentile)
    return _builder.to_yaml()


@mcp.tool()
async def runcard_run(output_dir: str | None = None) -> str:
    """Execute event generation using the current run card builder configuration.

    Args:
        output_dir: Directory for output files. Defaults to a temporary directory.

    Returns:
        A summary of the generation run including output file location and metadata.
    """
    try:
        validate_config(_builder.get_config())
    except ValidationError as e:
        return f"Validation error: {e}"

    run_card_yaml = _builder.to_yaml()
    result = await run_generation(run_card_yaml, output_dir=output_dir)
    return format_result(result)


@mcp.tool()
async def runcard_validate() -> str:
    """Validate the current run card for physics consistency.

    Checks required fields, parameter ranges, nucleus/model compatibility,
    cascade settings, and cut definitions without running the generator.
    """
    try:
        validate_config(_builder.get_config())
    except ValidationError as e:
        return f"Validation error: {e}"
    return "Configuration is valid."


# ---------------------------------------------------------------------------
# MCP Resources
# ---------------------------------------------------------------------------


@mcp.resource("achilles://version")
async def get_version() -> str:
    """Get the Achilles version string."""
    try:
        from achilles import __version__
        return __version__
    except ImportError:
        return "Achilles not installed — install with: pip install 'Achilles[mcp]'"


@mcp.resource("achilles://beam-types")
async def list_beam_types() -> str:
    """List supported beam/probe types with their PDG IDs and process configurations."""
    return json.dumps(BEAM_TYPES, indent=2)


@mcp.resource("achilles://nuclei")
async def list_nuclei() -> str:
    """List all available nucleus configurations."""
    results = {}
    for name in AVAILABLE_NUCLEI:
        config_path = _DATA_DIR / "default" / f"{name}.yml"
        if config_path.exists():
            results[name] = config_path.read_text()
        else:
            results[name] = f"Config file not found: {config_path}"
    return json.dumps(results, indent=2)


@mcp.resource("achilles://nuclei/{name}")
async def get_nucleus(name: str) -> str:
    """Get the default configuration for a specific nucleus."""
    config_path = _DATA_DIR / "default" / f"{name}.yml"
    if config_path.exists():
        return config_path.read_text()
    return f"Nucleus '{name}' not found. Available: {AVAILABLE_NUCLEI}"


@mcp.resource("achilles://form-factors")
async def list_form_factors() -> str:
    """List available form factor models and their parameters."""
    ff_path = _ACHILLES_ROOT / "FormFactors.yml"
    if ff_path.exists():
        return ff_path.read_text()
    ff_path = _DATA_DIR / "default" / "FormFactors.yml"
    if ff_path.exists():
        return ff_path.read_text()
    return json.dumps({
        "vector": FORM_FACTOR_VECTOR,
        "axial": FORM_FACTOR_AXIAL,
        "coherent": FORM_FACTOR_COHERENT,
    }, indent=2)


@mcp.resource("achilles://nuclear-models")
async def list_nuclear_models() -> str:
    """List available nuclear model types and Fortran model names."""
    return json.dumps({
        "models": NUCLEAR_MODELS,
        "fortran_model_names": FORTRAN_MODEL_NAMES,
        "output_formats": OUTPUT_FORMATS,
        "cascade_probability_types": CASCADE_PROBABILITY_TYPES,
    }, indent=2)


@mcp.resource("achilles://flux-files")
async def list_flux_files() -> str:
    """List available flux files that can be used with Spectrum beams.

    Use these paths with the histogram or hepdata parameters in set_beam.
    .dat files are used with the histogram parameter.
    .yaml files are used with the hepdata parameter.
    """
    if not _FLUX_DIR.exists():
        return "No flux directory found."
    flux_files = {}
    for p in sorted(_FLUX_DIR.iterdir()):
        if p.suffix in (".dat", ".yaml", ".csv"):
            flux_files[p.name] = f"flux/{p.name}"
    return json.dumps(flux_files, indent=2)


@mcp.resource("achilles://examples")
async def list_examples() -> str:
    """List all available example run cards."""
    if not _EXAMPLES_DIR.exists():
        return "No examples directory found."
    examples = sorted(p.name for p in _EXAMPLES_DIR.glob("*.yml"))
    return json.dumps(examples, indent=2)


@mcp.resource("achilles://examples/{name}")
async def get_example(name: str) -> str:
    """Get the content of a specific example run card."""
    example_path = _EXAMPLES_DIR / name
    if not example_path.exists():
        # Try adding .yml extension
        example_path = _EXAMPLES_DIR / f"{name}.yml"
    if example_path.exists():
        return example_path.read_text()
    return f"Example '{name}' not found. Use achilles://examples to list available examples."


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
