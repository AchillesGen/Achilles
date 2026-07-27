<!--
SPDX-FileCopyrightText: 2018-2026 Achilles Developers
SPDX-License-Identifier: GPL-3.0-or-later
-->

## Build performance

CMake wires up `ccache` automatically when it is installed, but its defaults are a poor
fit for a project this size. Recommended one-time setup:

```bash
ccache -M 20G                          # the 5G default evicts constantly
ccache -o compiler_check=content       # the mtime default invalidates the whole cache
                                       # whenever the toolchain is touched
ccache -o sloppiness=locale,time_macros
```

Check your hit rate with `ccache -s`. Note that `ccache` does not appear in
`compile_commands.json` — CMake strips compiler launchers from the compile database.

Two build options trade debuggability and tooling for speed:

- `ACHILLES_DEBUG_INFO_LEVEL` (default `1`) compiles `RelWithDebInfo` with `-g1` rather
  than `-g`. Line tables are kept, so backtraces, `perf` and `addr2line` still work, but
  full type and local-variable DWARF is dropped. For debugger stepping, use
  `-DCMAKE_BUILD_TYPE=Debug` or `-DACHILLES_DEBUG_INFO_LEVEL=2`.
- `ACHILLES_USE_FAST_LINKER` (default `ON`) uses `mold` or `lld` when either is
  installed, falling back to the system linker.
