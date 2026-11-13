// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "Achilles/Logging.hh"
#include "Achilles/Splash.hh"
#include "Achilles/Version.hh"
#include "fmt/core.h"
#include "git.h"

inline void Splash() {
    achilles::PrintSplash(ACHILLES_VERSION);
}

inline void GitInformation() {
    // std::string msg;
    if(git::IsPopulated()) {
        spdlog::info("Achilles git information");
        spdlog::info("    Commit: {}", git::CommitSHA1());
        spdlog::info("    Branch: {}", git::Branch());
        spdlog::info("    Local Changes: {}", git::AnyUncommittedChanges() ? "Yes" : "No");
    } else {
        spdlog::warn("This is not a git repository version of Achilles");
    }
}
