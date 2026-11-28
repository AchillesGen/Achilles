// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Logging.hh"
#include "Achilles/Logo.hh"
#include "Achilles/RunCascade.hh"
#include "Achilles/System.hh"
#include "Achilles/Version.hh"

#include "docopt.h"
#include "fmt/color.h"

static const std::string USAGE =
    R"(
    Usage:
      achilles-cascade [<input>] [-v | -vv] [-l | -ll]
      achilles-cascade (-h | --help)
      achilles-cascade --version
      achilles-cascade --version-full

    Options:
      -v[v]            Increase verbosity level.
      -l[l]            Increase log verbosity 
                       (Note: Log verbosity is never lower than total level)
      -h --help        Show this screen.
      --version        Show the version and the git tag or commit.
      --version-full   Show the version and all build details (include these in bug reports).
)";

int main(int argc, char *argv[]) {
    // The data search path in the --version-full report needs the paths set up first
    achilles::PathVariables::Instance().InitializeStandalone();
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc},
                       true,                              // show help if requested
                       ShortVersion("achilles-cascade")); // version string

    if(args["--version-full"].asBool()) {
        fmt::print("{}\n", BuildReport("achilles-cascade"));
        return 0;
    }

    Splash();
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::red),
               "NOTE: You are running in cascade test mode. Hard interactions "
               "are not included\n");

    std::string runcard = "cascade.yml";
    if(args["<input>"].isString()) runcard = args["<input>"].asString();

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    auto log_verbosity = std::min(verbosity, static_cast<int>(2 - args["-l"].asLong()));
    CreateLogger(verbosity, log_verbosity, 5);
    GitInformation();

    achilles::CascadeTest::RunCascade(runcard);

    return 0;
}
