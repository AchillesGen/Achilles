// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/AchillesInit.hh"
#include "Achilles/Cuts.hh"
#include "Achilles/EventGen.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Logging.hh"
#include "Achilles/Logo.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/OneParticleCuts.hh"
#include "Achilles/Particle.hh"
#include "Achilles/System.hh"
#include "Achilles/TwoParticleCuts.hh"
#include "Achilles/Version.hh"
#include "Achilles/fortran/FNuclearModel.hh"
#include "fmt/core.h"
#ifdef ACHILLES_SHERPA_INTERFACE
#include "Plugins/Sherpa/Channels.hh"
#include "Plugins/Sherpa/SherpaInterface.hh"
#endif

#include "docopt.h"

#include <algorithm>

// docopt only parses sections whose title ends in "options:"
static const std::string USAGE =
    R"(
    Achilles: A CHIcago Land Lepton Event Simulator

    Usage:
      achilles [<input>] [-v | -vv | -q] [-l | -ll] [-b] [--ascii] [--logfile=<logfile_name>] [-s | --sherpa=<sherpa>...]
      achilles --display-cuts
      achilles --display-ps
      achilles --display-ff
      achilles --display-int-models
      achilles --display-nuc-models
      achilles (-h | --help)
      achilles --version
      achilles --version-full

    <input> is the run card to use (default: run.yml in the current directory).

    Run options:
      -q --quiet                            Show only the progress line, warnings and the final
                                              summary (the log file keeps everything)
      -v[v]                                 Show more: -v for debug output, -vv for trace output
      -l[l]                                 Write more to the log file only: -l for debug,
                                              -ll for trace (never less than the terminal)
      -b                                    Batch mode: no live progress line, for log files
                                              and clusters
      --ascii                               Draw tables and progress bars with plain ASCII
                                              (also set by the ACHILLES_ASCII environment variable)
      --logfile=<logfile_name>              Write the log to this file instead of achilles.log
      -s <sherpa> --sherpa=<sherpa>         Pass an option to Sherpa (can be repeated)

    Information options:
      --display-cuts                        List the cuts that can be used in a run card
      --display-ps                          List the phase space mappings
      --display-ff                          List the form factors
      --display-int-models                  List the cascade interaction models
      --display-nuc-models                  List the nuclear models
      --version                             Show the version and the git tag or commit
      --version-full                        Show the version and all build details
                                              (include these in bug reports)
      -h --help                             Show this screen

    Examples:
      achilles                              Run ./run.yml (a default one is copied here if missing)
      achilles my_run.yml -b                Run in batch mode, e.g. on a cluster
      achilles my_run.yml -q                Show only the progress and the final summary
      achilles --display-nuc-models         List the nuclear models for the run card
      achilles --version-full               Build details to paste into a bug report
)";

int main(int argc, char *argv[]) {
    achilles::PathVariables::Instance().InitializeStandalone();
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc},
                       true,                      // show help if requested
                       ShortVersion("achilles")); // version string

    if(args["--version-full"].asBool()) {
        fmt::print("{}\n", BuildReport("achilles"));
        return 0;
    }

    achilles::console::Configure(args["-b"].asBool(), args["--ascii"].asBool());

    // Quiet shows only warnings on the terminal; the log file stays at info
    const bool quiet = args["--quiet"].asBool();
    int verbosity = quiet ? 3 : static_cast<int>(2 - args["-v"].asLong());
    int log_verbosity = std::min({verbosity, 2, static_cast<int>(2 - args["-l"].asLong())});

    std::string logFilePath = "achilles.log";
    if(args["--logfile"].isString()) logFilePath = args["--logfile"].asString();

    if(args["--display-cuts"].asBool()) {
        achilles::CutFactory<achilles::OneParticleCut>::Display();
        achilles::CutFactory<achilles::TwoParticleCut>::Display();
        return 0;
    }

    if(args["--display-ps"].asBool()) {
        achilles::Factory<achilles::HadronicBeamMapper, const achilles::ProcessInfo &,
                          size_t>::Display();
        achilles::Factory<achilles::FinalStateMapper, std::vector<double>>::Display();
#ifdef ACHILLES_SHERPA_INTERFACE
        achilles::Factory<PHASIC::Channels, std::vector<double>>::Display();
#endif // ACHILLES_SHERPA_INTERFACE
        return 0;
    }

    if(args["--display-ff"].asBool()) {
        achilles::FormFactorFactory::Display();
        return 0;
    }

    if(args["--display-int-models"].asBool()) {
        achilles::InteractionFactory::Display();
        return 0;
    }

    achilles::FortranModel::RegisterModels();
    if(args["--display-nuc-models"].asBool()) {
        achilles::FortranModel::DisplayModels();
        return 0;
    }

    if(!quiet) Splash();

    std::string runcard = "run.yml";
    if(args["<input>"].isString()) runcard = args["<input>"].asString();

    std::vector<std::string> shargs;
    if(args["--sherpa"].isStringList()) shargs = args["--sherpa"].asStringList();

    bool batchMode = args["-b"].asBool();

    achilles::InitializeLogging(verbosity, log_verbosity, logFilePath);
    achilles::GenerateEvents(runcard, shargs, batchMode);

    return 0;
}
