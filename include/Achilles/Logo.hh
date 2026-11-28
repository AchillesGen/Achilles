// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "Achilles/BuildInfo.hh"
#include "Achilles/GitInfo.hh"
#include "Achilles/Splash.hh"
#include "Achilles/Version.hh"
#include "fmt/format.h"
#include "spdlog/spdlog.h"

#include <regex>
#include <string>
#include <vector>

inline void Splash() {
    achilles::PrintSplash(ACHILLES_VERSION);
}

inline std::string BranchName(const achilles::GitInfo &git) {
    return git.branch.empty() || git.branch == "HEAD" ? "detached HEAD" : git.branch;
}

inline std::string FilesChanged(size_t nfiles) {
    return fmt::format("{} file{} changed", nfiles, nfiles == 1 ? "" : "s");
}

// Without a tag, git describe falls back to the abbreviated commit
inline bool DescribeHasTag(const achilles::GitInfo &git) {
    return !git.describe.empty() && git.commit.rfind(git.describe.substr(0, 7), 0) != 0;
}

// The release tag the build was made from, empty when it is not on a tag
inline std::string GitTag(const achilles::GitInfo &git) {
    if(git.source == "release archive") return git.branch;
    if(!DescribeHasTag(git)) return "";
    // "tag[-<commits since>-g<commit>][-dirty]"
    static const std::regex describe{R"((.+?)(?:-(\d+)-g[0-9a-f]+)?(?:-dirty)?)"};
    std::smatch match;
    std::regex_match(git.describe, match, describe);
    if(!match[2].matched) return match[1].str();
    return fmt::format("{} (+{} commits)", match[1].str(), match[2].str());
}

inline std::vector<achilles::InfoRow> GitRows() {
    const auto &git = achilles::BuildGitInfo();
    if(!git.available)
        return {{"Git", "no git information (not built from a checkout, package or release)"}};

    std::vector<achilles::InfoRow> rows{{"Commit", git.commit}};
    // A release archive has no branch, only the tags on its commit
    if(git.source == "release archive")
        rows.push_back({"Tags", git.branch.empty() ? "none" : git.branch});
    else
        rows.push_back({"Branch", BranchName(git)});
    if(DescribeHasTag(git)) rows.push_back({"Describe", git.describe});
    rows.push_back({"Commit date", git.date});
    rows.push_back({"Commit subject", git.subject});
    if(!git.remote.empty()) rows.push_back({"Remote", git.remote});
    if(git.source != "git checkout") rows.push_back({"Git info from", git.source});
    if(git.shallow) rows.push_back({"Clone", "shallow (tags and history may be missing)"});

    if(!git.dirty) {
        rows.push_back({"Local changes", "none"});
        return rows;
    }
    const size_t nfiles = git.changed_files.size();
    rows.push_back({"Local changes", fmt::format("{}, diff sha256 {}", FilesChanged(nfiles),
                                                 git.diff_hash.substr(0, 16))});
    // Enough to see what was touched without flooding the report
    constexpr size_t max_files = 10;
    for(size_t i = 0; i < nfiles && i < max_files; ++i) rows.push_back({"", git.changed_files[i]});
    if(nfiles > max_files) rows.push_back({"", fmt::format("... and {} more", nfiles - max_files)});
    return rows;
}

// Everything worth including in a bug report
inline std::vector<achilles::InfoRow> BuildReportRows() {
    auto rows = achilles::BuildRows();
    const auto git = GitRows();
    // Keep the version first, then the source tree, then the build and the system
    rows.insert(rows.begin() + 1, git.begin(), git.end());
    const auto system = achilles::SystemRows();
    rows.insert(rows.end(), system.begin(), system.end());
    return rows;
}

inline std::vector<achilles::InfoRow> KeyBuildRows() {
    const auto &git = achilles::BuildGitInfo();
    if(!git.available) return {{"Version", ACHILLES_VERSION}, {"Git", "no git information"}};

    const std::string tag = GitTag(git);
    const std::string branch = git.source == "release archive" ? "release" : BranchName(git);
    std::vector<achilles::InfoRow> rows{
        {"Version", ACHILLES_VERSION},
        {"Git tag", tag.empty() ? "none" : tag},
        {"Commit", fmt::format("{} ({})", git.commit.substr(0, 12), branch)}};
    if(git.source != "git checkout") rows.push_back({"Git info from", git.source});
    // The files are listed with -v and by achilles --version-full
    rows.push_back({"Local changes", git.dirty ? FilesChanged(git.changed_files.size()) : "none"});
    return rows;
}

// One line for --version
inline std::string ShortVersion(const std::string &program) {
    const auto &git = achilles::BuildGitInfo();
    if(!git.available) return fmt::format("{} {}", program, ACHILLES_VERSION);
    const std::string ref = !git.describe.empty() ? git.describe : git.commit.substr(0, 12);
    return fmt::format("{} {} ({})", program, ACHILLES_VERSION, ref);
}

// Plain text for --version-full, meant to be pasted into bug reports as is
inline std::string BuildReport(const std::string &program) {
    return fmt::format("{} {}\n\n{}\nPlease include this output when reporting a problem.", program,
                       ACHILLES_VERSION, achilles::InfoText(BuildReportRows()));
}

// The full build information only with -v (and in a log file that records debug output)
inline void GitInformation() {
    achilles::PrintInfo("Build (run with -v or achilles --version-full for details)",
                        KeyBuildRows());
    achilles::PrintInfo("Build information", BuildReportRows(), true);

    const auto &git = achilles::BuildGitInfo();
    if(git.available && git.dirty)
        spdlog::warn("This build includes uncommitted changes to the Achilles source; "
                     "results may differ from commit {}",
                     git.commit.substr(0, 12));
}
