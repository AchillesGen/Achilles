// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Logging.hh"
#include "Achilles/Logo.hh"
#include "Achilles/Precomputed.hh"
#include "Achilles/Version.hh"

#include "docopt.h"
#include "fmt/color.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

static const std::string USAGE =
    R"(
    Usage:
      achilles-precomputed [<input>] [-v | -vv]
      achilles-precomputed (-h | --help)
      achilles-precomputed --version

    Options:
      -v[v]            Increase verbosity level.
      -h --help        Show this screen.
      --version        Show version.
)";

int main(int argc, char *argv[]) {
    Splash();
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::red),
               "NOTE: You are running in precomputed mode. Hard "
               "interactions are not included\n");
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc},
                       true,                                          // show help if requested
                       fmt::format("Achilles {}", ACHILLES_VERSION)); // version string

    std::string runcard = "precomputed.yml";
    if(args["<input>"].isString()) runcard = args["<input>"].asString();

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    CreateLogger(verbosity, verbosity, 5, "achilles-precomputed.log");

    achilles::Precomputed::RunCascade runner(runcard);
    runner.RunAll();

    return 0;
}
