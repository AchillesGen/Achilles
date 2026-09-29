#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Decide what a physval workflow run does; prints GITHUB_OUTPUT lines.

  push       compare, if a commit *subject* carries !physval or !physval(<scope>);
             the scope is comma-separated setup names and/or dry-run
  schedule   baseline, skipped when the stored one is already current
  dispatch   whatever the inputs say
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys

import yaml

# A change under these paths since the stored baseline's commit makes it stale.
PHYSICS_PATHS = ["src", "include", "data", "examples", "flux", "validation/physval",
                 "CMakeLists.txt", "CMake"]
EVENTS = {"compare": 500_000, "baseline": 2_000_000}


def _git(*args) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], capture_output=True, text=True)


def parse_markers(messages) -> tuple:
    """(marked, dry_run, setups) from commit messages, reading subject lines only."""
    marked, dry, setups = False, False, []
    for msg in messages:
        subject = msg.split("\n", 1)[0]
        if "!physval" not in subject:
            continue
        marked = True
        for scope in re.findall(r"!physval\(([^)]*)\)", subject):
            for tok in (t.strip() for t in scope.split(",")):
                if tok == "dry-run":
                    dry = True
                elif tok and tok not in setups:
                    setups.append(tok)
    return marked, dry, setups


def baseline_is_current(branch: str, key: str, setups) -> bool:
    """True if every setup has a baseline for ``key`` from a commit with no physics
    change since."""
    if _git("fetch", "--depth", "1", "origin", branch).returncode:
        return False
    shas = set()
    for s in setups:
        shown = _git("show", f"FETCH_HEAD:baselines/{key}/{s}.json")
        if shown.returncode:
            return False
        shas.add(json.loads(shown.stdout).get("main_sha", ""))
    for sha in shas:
        if (_git("fetch", "--depth", "1", "origin", sha).returncode
                or _git("diff", "--quiet", sha, "HEAD", "--", *PHYSICS_PATHS).returncode):
            return False
    return True


def main() -> int:
    env = os.environ
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "measurements.yml")) as fh:
        names = [e["name"] for e in yaml.safe_load(fh)["experiments"]]

    event = env.get("EVENT", "")
    run, mode, dry, only = True, "compare", False, []
    if event == "workflow_dispatch":
        mode = env.get("IN_MODE") or "compare"
        dry = env.get("IN_DRY", "").lower() == "true"
        only = [t.strip() for t in env.get("IN_ONLY", "").split(",") if t.strip()]
    elif event == "schedule":
        mode = "baseline"
    else:
        run, dry, only = parse_markers(json.loads(env.get("COMMIT_MESSAGES") or "[]"))

    if mode == "baseline" and not dry and env.get("REF") != "refs/heads/main":
        sys.exit("::error::a real baseline can only be made from main")

    unknown = [o for o in only if o not in names]
    if unknown:
        sys.exit(f"::error::unknown physval setups: {', '.join(unknown)}. "
                 f"Valid: {', '.join(names)}")
    if only:
        names = [n for n in names if n in only]

    if run and event == "schedule" and baseline_is_current(
            env["BASELINE_BRANCH"], env["KEY"], names):
        print("baseline is current; nothing to do", file=sys.stderr)
        run = False

    events = env.get("IN_EVENTS") or str(EVENTS[mode])
    print(f"run={'true' if run else 'false'}")
    print(f"mode={mode}")
    print(f"dry_run={'true' if dry else 'false'}")
    print(f"only={','.join(only)}")
    print(f"events={events}")
    print("matrix=" + json.dumps({"experiment": names}))
    print(f"run={run} mode={mode} dry_run={dry} events={events} "
          f"setups={','.join(only) or 'all'}", file=sys.stderr)
    return 0


def _selftest() -> int:
    checks = {
        "subject marker": parse_markers(["fix: x !physval"]) == (True, False, []),
        "body is ignored": parse_markers(["fix: x\n\n!physval(A)"]) == (False, False, []),
        "scoped": parse_markers(["x !physval(A, dry-run)", "y !physval(B,A)"])
                  == (True, True, ["A", "B"]),
    }
    for name, ok in checks.items():
        print(f"[{'ok' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(_selftest() if "--selftest" in sys.argv else main())
