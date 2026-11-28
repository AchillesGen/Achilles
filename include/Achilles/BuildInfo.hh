// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef BUILD_INFO_HH
#define BUILD_INFO_HH

#include <string>
#include <vector>

namespace achilles {

// One line of the build information, e.g. {"Compiler", "GCC 13.2.0"}. A row with an empty
// name continues the row above it (for example a list of files).
struct InfoRow {
    std::string name, value;
};

// How this copy of Achilles was built: version, build type, compilers and features
std::vector<InfoRow> BuildRows();
// Where it is running: operating system and the paths searched for data files
std::vector<InfoRow> SystemRows();

// Plain text with one "name: value" line per row, for pasting into bug reports
std::string InfoText(const std::vector<InfoRow> &rows);
// The same rows under a title, printed to the terminal and the log file
void PrintInfo(const std::string &title, const std::vector<InfoRow> &rows, bool debug = false);

} // namespace achilles

#endif
