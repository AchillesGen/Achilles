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
#include "git.h"
#ifdef ACHILLES_SHERPA_INTERFACE
#include "Plugins/Sherpa/Channels.hh"
#include "Plugins/Sherpa/SherpaInterface.hh"
#endif

#include "docopt.h"

static const std::string USAGE =
    R"(
    Usage:
      achilles [<input>] [-v | -vv] [-l | -ll] [-b] [--logfile=<logfile_name>] [-s | --sherpa=<sherpa>...]
      achilles --display-cuts
      achilles --display-ps
      achilles --display-ff
      achilles --display-int-models
      achilles --display-nuc-models
      achilles (-h | --help)
      achilles --version

    Options:
      -v[v]                                 Increase verbosity level.
      -l[l]                                 Increase log verbosity
                                              (Note: Log verbosity is never lower than total level)
      -b                                    Batch Mode (makes output more log-friendly)
      -h --help                             Show this screen.
      --version                             Show version.
      --logfile=<logfile_name>              Change the logging output destination
      -s <sherpa> --sherpa=<sherpa>         Define Sherpa option.
      --display-cuts                        Display the available cuts
      --display-ps                          Display the available phase spaces
      --display-ff                          Display the available form factors
      --display-int-models                  Display the available cascade interaction models
      --display-nuc-models                  Display the available nuclear interaction models
)";

int main(int argc, char *argv[]) {
    Splash();
    achilles::PathVariables::Instance().InitializeStandalone();
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc},
                       true,                                          // show help if requested
                       fmt::format("achilles {}", ACHILLES_VERSION)); // version string

    int verbosity = static_cast<int>(2 - args["-v"].asLong());
    int log_verbosity = std::min(verbosity, static_cast<int>(2 - args["-l"].asLong()));

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
        achilles::NuclearModelFactory::Display();
        achilles::FortranModel::DisplayModels();
        return 0;
    }

    std::string runcard = "run.yml";
    if(args["<input>"].isString()) runcard = args["<input>"].asString();

    std::vector<std::string> shargs;
    if(args["--sherpa"].isStringList()) shargs = args["--sherpa"].asStringList();

    bool batchMode = args["-b"].asBool();

    achilles::InitializeLogging(verbosity, log_verbosity, logFilePath);
    achilles::GenerateEvents(runcard, shargs, batchMode);

    return 0;
}
