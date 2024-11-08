#include "Achilles/Logging.hh"
#include "Achilles/Precomputed.hh"
#include "Achilles/Version.hh"

#include "docopt.h"
#include "fmt/color.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

void Splash() {
    fmt::print(R"splash(
+=====================================================================+
|                                                                     | 
|    .d8b.   .o88b. db   db d888888b db      db      d88888b .d8888.  |
|   d8' `8b d8P  Y8 88   88   `88'   88      88      88'     88'  YP  |
|   88ooo88 8P      88ooo88    88    88      88      88ooooo `8bo.    |
|   88~~~88 8b      88~~~88    88    88      88      88~~~~~  `Y8b.   |
|   88   88 Y8b  d8 88   88   .88.   88booo. 88booo. 88.     db  8D   |
|   YP   YP  `Y88P' YP   YP Y888888P Y88888P Y88888P Y88888P `8888Y'  |
|                                                                     | 
+---------------------------------------------------------------------+
|                                                                     |
|                 .d88888888888888888888888888b.                      |
|            .d8888888888888888888888888888888888b.                   |
|         .d88888888888888888888888888888888888888888b.               |
|      88888888888888888888888888888888888888888888888888b.           |
|       `Y8888888888888888888888888888888888888888888888888b.         |
|            `Y88888888888888888888888888888888888888888888888        |
|              `Y88888888888888888888888888888888888888888888888      |
|                 `Y888888888888888888888888888P  Y88888888888888     |
|                    `Y88  8888  8888  88Y'        Y8888888888888     |
|                   .d888888888888888888888b.       Y888888888888.    |
|                .d88888888888888888888888888b.      8888888888888    |
|              .d8888888888888888888888888888888b   88888888888888    |
|            .d8888888888888888888888888888888888   8888888888888D    |
|            d88888888888888888888888888888888888   8888888888888     |
|           d888888888888888888888888888888888888   Y888888888888     |
|          d888888P'   d888888888888888888888888P       88888888      |
|          8888P     .d888888888888888888888888P         888888       |
|          888P    .d8888888888888888888888888P          8888Y        |
|          888b  .d888888888888888888888888888           8888P        |
|          Y88  d888888888888888888888888888888           888         |
|          `8' d8888888888888888888888888888888            88         |
|              8888888888888888888888888888888P             8         |
|              8888888888888888P     Y8888888P                        |
|             d8888888888888888       Y88888P                         |
|             8888888888888888P        Y888P                          |
|            d8888888888888P            Y8P                           |
|            888888888P                  8                            |
|           d888P                                                     |
|                                                                     |
+---------------------------------------------------------------------+
|                                                                     |
|    Version: {:56}|
|    Authors: Joshua Isaacson, William Jay, Alessandro Lovato,        | 
|        Pedro A. Machado, Noemi Rocco                                | 
|                                                                     |
|    Undergraduate Student Contributions:                             |
|        Diego Lopez Gutierrez, Sherry Wang, Russell Farnsworth       |
|                                                                     |
+=====================================================================+
)splash",
               ACHILLES_VERSION);

    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::red),
               "NOTE: You are running in precomputed mode. Hard "
               "interactions are not included\n");
}

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
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc},
                       true,                                          // show help if requested
                       fmt::format("Achilles {}", ACHILLES_VERSION)); // version string

    std::string runcard = "precomputed.yml";
    if(args["<input>"].isString()) runcard = args["<input>"].asString();

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    CreateLogger(verbosity, 5);

    achilles::Precomputed::RunCascade runner(runcard);
    runner.RunAll();

    return 0;
}
