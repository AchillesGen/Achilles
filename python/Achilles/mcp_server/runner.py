"""Subprocess-based wrapper for running Achilles event generation."""

import asyncio
import os
import sys
import tempfile
from pathlib import Path

import yaml


class _IncludeLoader(yaml.SafeLoader):
    """YAML loader that handles !include tags by returning the path as a string."""
    pass


_IncludeLoader.add_constructor(
    "!include",
    lambda loader, node: f"!include:{loader.construct_scalar(node)}",
)


async def run_generation(run_card_yaml: str, output_dir: str | None = None) -> dict:
    """Write a run card to a temp directory and run Achilles in a subprocess.

    Uses subprocess isolation so C++ global state (logger, paths) is fresh
    each invocation.

    Returns a dict with run metadata: output file path, size, return code, etc.
    """
    work_dir = output_dir or tempfile.mkdtemp(prefix="achilles_mcp_")
    os.makedirs(work_dir, exist_ok=True)

    run_card_path = os.path.join(work_dir, "run.yml")
    with open(run_card_path, "w") as f:
        f.write(run_card_yaml)

    # Parse output file name from the run card (using include-aware loader)
    config = yaml.load(run_card_yaml, Loader=_IncludeLoader)
    output_name = config.get("Main", {}).get("Output", {}).get("Name", "achilles.hepmc")
    output_format = config.get("Main", {}).get("Output", {}).get("Format", "NuHepMC")
    nevents = config.get("Main", {}).get("NEvents", 0)

    script = (
        "import os, sys\n"
        f"os.chdir({work_dir!r})\n"
        "from Achilles import initialize_logging, generate_events\n"
        "initialize_logging(verbosity=3, log_verbosity=2, logfile='achilles_mcp.log')\n"
        f"generate_events(runcard={run_card_path!r}, batch_mode=True)\n"
    )

    proc = await asyncio.create_subprocess_exec(
        sys.executable, "-c", script,
        cwd=work_dir,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await proc.communicate()

    output_path = os.path.join(work_dir, output_name)
    output_exists = os.path.isfile(output_path)
    output_size = os.path.getsize(output_path) if output_exists else 0

    log_path = os.path.join(work_dir, "achilles_mcp.log")
    log_tail = ""
    if os.path.isfile(log_path):
        with open(log_path) as lf:
            lines = lf.readlines()
            log_tail = "".join(lines[-30:])

    result = {
        "return_code": proc.returncode,
        "work_dir": work_dir,
        "run_card_path": run_card_path,
        "output_file": output_path if output_exists else None,
        "output_format": output_format,
        "output_size_bytes": output_size,
        "nevents_requested": nevents,
        "log_tail": log_tail,
    }

    if proc.returncode != 0:
        result["stderr"] = stderr.decode(errors="replace")[-2000:]

    return result


def format_result(result: dict) -> str:
    """Format a run result dict into a human-readable summary."""
    lines = []
    if result["return_code"] == 0:
        lines.append("Event generation completed successfully.")
    else:
        lines.append(
            f"Event generation failed (exit code {result['return_code']})."
        )

    lines.append(f"Working directory: {result['work_dir']}")
    lines.append(f"Run card: {result['run_card_path']}")
    lines.append(f"Events requested: {result['nevents_requested']}")

    if result["output_file"]:
        size_mb = result["output_size_bytes"] / (1024 * 1024)
        lines.append(
            f"Output: {result['output_file']} "
            f"({result['output_format']}, {size_mb:.2f} MB)"
        )
    else:
        lines.append("Output: No output file produced.")

    if result.get("stderr"):
        lines.append(f"\nStderr (last 2000 chars):\n{result['stderr']}")

    if result["log_tail"]:
        lines.append(f"\nLog (last 30 lines):\n{result['log_tail']}")

    return "\n".join(lines)
