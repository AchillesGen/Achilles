#include "Achilles/EventGen.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Logging.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Particle.hh"
#include "Achilles/System.hh"
#include "Achilles/Logo.hh"
#include "Achilles/Version.hh"
#include "Achilles/fortran/FNuclearModel.hh"
#include "git.h"
#include "plugins/Manager/PluginManager.hh"
#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/Channels.hh"
#include "plugins/Sherpa/SherpaInterface.hh"
#endif

#include "docopt.h"

#include <dlfcn.h>

using namespace achilles::SystemVariables;
using namespace achilles::PathVariables;

static const std::string USAGE =
    R"(
    Usage:
      achilles [<input>] [-v | -vv] [-s | --sherpa=<sherpa>...]
      achilles --display-cuts
      achilles --display-ps
      achilles --display-ff
      achilles --display-int-models
      achilles --display-nuc-models
      achilles (-h | --help)
      achilles --version

    Options:
      -v[v]                                 Increase verbosity level.
      -h --help                             Show this screen.
      --version                             Show version.
      -s <sherpa> --sherpa=<sherpa>         Define Sherpa option.
      --display-cuts                        Display the available cuts
      --display-ps                          Display the available phase spaces
      --display-ff                          Display the available form factors
      --display-int-models                  Display the available cascade interaction models
      --display-nuc-models                  Display the available nuclear interaction models
)";

void GenerateEvents(const std::string &runcard, const std::vector<std::string> &shargs) {
    achilles::EventGen generator(runcard, shargs);
    generator.Initialize();
    generator.GenerateEvents();
}

int main(int argc, char *argv[]) {
    Splash();
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE, {argv + 1, argv + argc},
                       true,                                          // show help if requested
                       fmt::format("achilles {}", ACHILLES_VERSION)); // version string

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    CreateLogger(verbosity, 5);
    GitInformation();
    achilles::Plugin::Manager plugin_manager;

    if(args["--display-cuts"].asBool()) {
        achilles::CutFactory<achilles::OneParticleCut>::DisplayCuts();
        achilles::CutFactory<achilles::TwoParticleCut>::DisplayCuts();
        return 0;
    }

    if(args["--display-ps"].asBool()) {
        achilles::PSFactory<achilles::HadronicBeamMapper, size_t>::DisplayPhaseSpaces();
        achilles::PSFactory<achilles::FinalStateMapper, std::vector<double>>::DisplayPhaseSpaces();
#ifdef ACHILLES_SHERPA_INTERFACE
        achilles::PSFactory<PHASIC::Channels, std::vector<double>>::DisplayPhaseSpaces();
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

    GenerateEvents(runcard, shargs);

    return 0;
}
