// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef GIT_INFO_HH
#define GIT_INFO_HH

#include <string>
#include <vector>

namespace achilles {

// The git checkout this copy of Achilles was built from, recorded at build time by
// CMake/WriteGitInfo.cmake. Builds from a Python package or a release tarball get it from
// the information recorded when the package was made; see "source".
struct GitInfo {
    bool available;
    std::string commit;   // full hash of HEAD
    std::string branch;   // "HEAD" for a detached checkout
    std::string describe; // git describe --tags --always --dirty
    std::string date;     // commit date, ISO 8601
    std::string subject;  // first line of the commit message
    std::string remote;   // URL of origin, without credentials
    bool dirty;           // tracked files had uncommitted changes
    bool shallow;         // shallow clone, so describe may lack tags
    // SHA-256 of "git diff HEAD": identifies the uncommitted changes without including them
    std::string diff_hash;
    std::vector<std::string> changed_files;
    // "git checkout", "package snapshot", "release archive" or "none"
    std::string source;
};

// Defined in the generated GitInfo.cc of the achilles_git_info library
const GitInfo &BuildGitInfo();

} // namespace achilles

#endif
