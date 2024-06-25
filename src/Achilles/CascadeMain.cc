#include "Achilles/Logging.hh"
#include "Achilles/Logo.hh"
#include "Achilles/RunCascade.hh"
#include "Achilles/Version.hh"

#include "docopt.h"
#include "fmt/color.h"

static const std::string USAGE =
    R"(
    Usage:
      achilles-cascade [<input>] [-v | -vv]
      achilles-cascade (-h | --help)
      achilles-cascade --version

    Options:
      -v[v]            Increase verbosity level.
      -h --help        Show this screen.
      --version        Show version.
)";

int main(int argc, char *argv[]) {
    Splash();
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::red),
               "NOTE: You are running in cascade test mode. Hard interactions "
               "are not included\n");
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc},
                       true,                                          // show help if requested
                       fmt::format("Achilles {}", ACHILLES_VERSION)); // version string

    std::string runcard = "cascade.yml";
    if(args["<input>"].isString()) runcard = args["<input>"].asString();

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    CreateLogger(verbosity, 5);
    GitInformation();

    achilles::CascadeTest::RunCascade(runcard);

    return 0;
}
