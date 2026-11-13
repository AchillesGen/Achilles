// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef CONTRIBUTORS_HH
#define CONTRIBUTORS_HH

#include <string>
#include <vector>

namespace achilles {

// People credited in the start-up banner.
struct Contributors {
    std::vector<std::string> authors, former_authors, undergraduate_contributors;

    // Throws std::runtime_error if the YAML is not in the expected format
    static Contributors Parse(const std::string &yaml);
    static const Contributors &Builtin();
    static const std::string &BuiltinYAML();
};

} // namespace achilles

#endif
